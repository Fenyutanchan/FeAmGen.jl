
#################################################################
"""
    degree( expr::Basic, x::Basic; max_n::Int64=100 )::Int64

A little extension for SymEngine. This function `degree` calculate the degree of polynomial `expr` of `x`.
For now we do not have toolkit to constraint `expr` as polynomial, so we only choose a maximum degree in the loop.
"""
function degree( expr::Basic, x::Basic; max_n::Int64=100 )::Int64
#################################################################

  diff_expr = expr
  for ii in 1:(max_n+1)

    diff_expr = diff(diff_expr,x)
    if diff_expr == 0
      return ii-1
    end # if

  end # for ii

  error( "Are you sure ",expr," is a polynomial of ",x,"? Or max_n = ",max_n," is too small?\n" )

end # function degree






###################################################
"""
    convert_to_array( ::Val{:Symbol}, mom::Basic )::Array{NamedTuple{(:num, :ki),Tuple{Basic,Basic}}}

This specific convertion is applied only on single symbol momentum, e.g. `k1`.

# Examples
```julia-repl
julia> using SymEngine, FeAmGen

julia> @vars k1
(k1,)

julia> FeAmGen.convert_to_array(k1)
1-element Array{NamedTuple{(:num, :ki),Tuple{Basic,Basic}},1}:
 (num = 1, ki = k1)
```
"""
function convert_to_array( ::Val{:Symbol}, mom::Basic )::Array{NamedTuple{(:num, :ki),Tuple{Basic,Basic}}}
###################################################

  return NamedTuple{(:num,:ki),Tuple{Basic,Basic}}[ ( num=one(Basic), ki=mom ) ] 

end # function convert_to_array



###################################################
"""
    convert_to_array( ::Val{:Mul}, mom::Basic )::Array{NamedTuple{(:num, :ki),Tuple{Basic,Basic}}}

This specific convertion is applied only on non-unit coefficient with single symbol momentum, e.g. `2*k1`.

# Examples
```julia-repl
julia> using SymEngine, FeAmGen

julia> @vars k1
(k1,)

julia> FeAmGen.convert_to_array(2*k1)
1-element Array{NamedTuple{(:num, :ki),Tuple{Basic,Basic}},1}:
 (num = 2, ki = k1)
```
"""
function convert_to_array( ::Val{:Mul}, mom::Basic )::Array{NamedTuple{(:num, :ki),Tuple{Basic,Basic}}}
###################################################

  arg_list = get_args(mom)
  @assert length(arg_list) == 2
  first_arg_class = SymEngine.get_symengine_class(arg_list[1])
  coeff, symbol = first_arg_class == :Integer ? arg_list : reverse(arg_list)

  return NamedTuple{(:num,:ki),Tuple{Basic,Basic}}[ ( num=coeff, ki=symbol ) ] 

end # function convert_to_array



###################################################
"""
    convert_to_array( ::Val{:Add}, mom::Basic )::Array{NamedTuple{(:num, :ki),Tuple{Basic,Basic}}}

This specific convertion is applied only on the summation of non-unit coefficient with single symbol momentum, e.g. `2*k1+3*k2+k3`.

# Examples
```julia-repl
julia> using SymEngine, FeAmGen

julia> @vars k1, k2, k3
(k1, k2, k3)

julia> FeAmGen.convert_to_array(2*k1+3*k2+k3)
3-element Array{NamedTuple{(:num, :ki),Tuple{Basic,Basic}},1}:
 (num = 1, ki = k3)
 (num = 2, ki = k1)
 (num = 3, ki = k2)
```
"""
function convert_to_array( ::Val{:Add}, mom::Basic )::Array{NamedTuple{(:num, :ki),Tuple{Basic,Basic}}}
###################################################

  arg_list = get_args(mom)
  n_solo = length(arg_list)
  @assert n_solo >= 2

  return map( arg_list ) do solo_mom

    solo_arg_list = get_args(solo_mom)
    @assert length(solo_arg_list) == 2 || SymEngine.get_symengine_class(solo_mom) == :Symbol

    if length(solo_arg_list) == 2

      first_arg_class = SymEngine.get_symengine_class(solo_arg_list[1])
      coeff, symbol = first_arg_class == :Integer ? solo_arg_list : reverse(solo_arg_list)
      return ( num=coeff, ki=symbol ) # return for map

    else # solo_mom is :Symbol

      return ( num=one(Basic), ki=solo_mom ) # return for map

    end # if

  end # do solo_mom


end # function convert_to_array







###################################################
"""
    convert_to_array( mom::Basic )::Array{NamedTuple{(:num, :ki),Tuple{Basic,Basic}}}

This is the interface to different specific convertion functions according to the class of `mom`.
Now we assume there only there different classes of momentum `mom`. For example, 
1. `k1`
2. `2*k1`
3. `2*k1+3*k2+k3`
"""
function convert_to_array( mom::Basic )::Array{NamedTuple{(:num, :ki),Tuple{Basic,Basic}}}
###################################################

  mom_class = SymEngine.get_symengine_class(mom)
   
  return convert_to_array( Val(mom_class), mom )

end # function convert_to_array






###################################################
"""
    make_SP( mom1::Basic, mom2::Basic)::Basic

This function linearly expands the scalar product of `mom1` and `mom2`. 

# Examples
```julia-repl
julia> using SymEngine, FeAmGen

julia> @vars k1, k2, k3
(k1, k2, k3)

julia> FeAmGen.make_SP( 2*k1+3*k2+k3, k1+k2 )
2*SP(k1, k1) + 2*SP(k1, k2) + 3*SP(k2, k1) + 3*SP(k2, k2) + SP(k3, k1) + SP(k3, k2)
```
"""
function make_SP( mom1::Basic, mom2::Basic )::Basic
###################################################

  @funs SP

  mom1_array = convert_to_array( mom1 )
  mom2_array = convert_to_array( mom2 )

  result_SP = zero(Basic)
  for pair1 in mom1_array, pair2 in mom2_array
    result_SP += pair1[:num] * pair2[:num] * SP( sort( Basic[pair1[:ki],pair2[:ki]], by=string )... )
  end # for pair1, pair2

  return result_SP

end # function make_SP





###################################################
"""
This funciton is deprecated and replaced by the above function `make_SP`.
The benchmark test to show the improvement of this replacement is in the following.

```julia
@vars k1, k2, k3, k4, q1, q2
test_mom = 2*k1+3*k2+k3+5*k4+6*q1+7*q2

println( "make_SP: " )
@btime aa = FeAmGen.make_SP_new( test_mom, test_mom )

println( "make_SP_old: " )
@btime aa = FeAmGen.make_SP( test_mom, test_mom )
```

The result is
```
make_SP: 
  92.063 μs (518 allocations: 10.96 KiB)
make_SP_old: 
  137.740 μs (1527 allocations: 55.01 KiB)
```
"""
function make_SP_old( mom1::Basic, mom2::Basic )::Basic
###################################################

  @funs SP

  class1 = SymEngine.get_symengine_class(mom1)
  class2 = SymEngine.get_symengine_class(mom2)

  if class1 == :Symbol && class2 == :Symbol

    return SP( sort( Basic[mom1,mom2], by=string )... )

  elseif class1 == :Symbol && class2 == :Mul

    arg_list = get_args(mom2)
    @assert length(arg_list) == 2
    first_arg_class = SymEngine.get_symengine_class(arg_list[1])
    coeff2, symbol2 = first_arg_class == :Integer ? arg_list : reverse(arg_list)
    return coeff2*SP( sort( Basic[mom1,symbol2], by=string )... )

  elseif class1 == :Mul && class2 == :Symbol

    return make_SP( mom2, mom1 )

  elseif class1 == :Mul && class2 == :Mul

    arg1_list = get_args(mom1)
    @assert length(arg1_list) == 2
    first_arg1_class = SymEngine.get_symengine_class(arg1_list[1])
    coeff1, symbol1 = first_arg1_class == :Integer ? arg1_list : reverse(arg1_list)

    arg2_list = get_args(mom2)
    @assert length(arg2_list) == 2
    first_arg2_class = SymEngine.get_symengine_class(arg2_list[1])
    coeff2, symbol2 = first_arg2_class == :Integer ? arg2_list : reverse(arg2_list)

    return coeff1*coeff2*SP( sort( Basic[symbol1,symbol2], by=string )... )

  elseif class1 == :Add 

    arg_list = get_args(mom1)
    return sum( map( m_ -> make_SP(m_,mom2), arg_list ) )

  elseif class2 == :Add

    return make_SP( mom2, mom1 )

  else

    error( "Exception found: $(mom1), $(mom2)" )

  end # if

  return SP(mom1,mom2)

end # function make_SP



###########################################
gen_mma_str( ::Val{:Symbol}, expr::Basic )::String = string(expr)
gen_mma_str( ::Val{:Integer}, expr::Basic )::String = string(expr)
gen_mma_str( ::Val{:Rational}, expr::Basic )::String = string(expr)
gen_mma_str( ::Val{:Complex}, expr::Basic )::String = string(expr)
###########################################


###########################################
function gen_mma_str( ::Val{:Pow}, expr::Basic )::String
###########################################
  arglist = get_args(expr)
  return "Power[ $(gen_mma_str(arglist[1])), $(arglist[2]) ]"
end # function gen_mma_str

###########################################
function gen_mma_str( ::Val{:FunctionSymbol}, expr::Basic )::String
###########################################
  name = replace( get_name(expr), "Trace" => "DiracTrace" )
  return "$(name)[ $(join( map( gen_mma_str, get_args(expr) ), "," )) ]"
end # function gen_mma_str

###########################################
function gen_mma_str( ::Val{:Mul}, expr::Basic )::String
###########################################
  return "Times[ $(join( (sort∘map)( gen_mma_str, get_args(expr) ), "," )) ]"
end # function gen_mma_str

###########################################
function gen_mma_str( ::Val{:Add}, expr::Basic )::String
###########################################
  return "Plus[ $(join( (sort∘map)( gen_mma_str, get_args(expr) ), "," )) ]"
end # function gen_mma_str


###########################################
"""
    gen_mma_str( expr::Basic )::String

This is a generic interface for different classes of the `expr`.
And it will generate the mathematica code for the expression `expr`.
"""
function gen_mma_str( expr::Basic )::String
###########################################

  expr_class = SymEngine.get_symengine_class(expr)
  return gen_mma_str( Val(expr_class), expr )

end # function gen_mma_str



###########################################
"""
    gen_mma_str_old( expr::Basic )::String

This is deprecated function for gen_mma_str.
Mostly the new version has similar efficiency than this old one. 

```julia
@vars k1, k2, k3, k4, q1, q2
test_mom = 2*k1+3*k2+k3+5*k4+6*q1+7*q2
expr = expand(test_mom^3+test_mom)

println( "FeAmGen.gen_mma_str(expr)" )
@btime FeAmGen.gen_mma_str(expr)
println( "FeAmGen.gen_mma_str_old(expr)" )
@btime FeAmGen.gen_mma_str_old(expr)
```

FeAmGen.gen_mma_str(expr)
  1.344 ms (2554 allocations: 115.23 KiB)
FeAmGen.gen_mma_str_old(expr)
  1.330 ms (2554 allocations: 115.23 KiB)

"""
function gen_mma_str_old( expr::Basic )::String
###########################################

  expr_class = SymEngine.get_symengine_class(expr)

  if expr_class == :Add
    return "Plus[ $(join( (sort∘map)( gen_mma_str, get_args(expr) ), "," )) ]"
  elseif expr_class == :Mul
    return "Times[ $(join( (sort∘map)( gen_mma_str, get_args(expr) ), "," )) ]"
  elseif expr_class == :FunctionSymbol
    name = replace( get_name(expr), "Trace" => "DiracTrace" )
    return "$(name)[ $(join( map( gen_mma_str, get_args(expr) ), "," )) ]"
  elseif expr_class == :Pow
    arglist = get_args(expr)
    return "Power[ $(gen_mma_str(arglist[1])), $(arglist[2]) ]"
  elseif expr_class in [:Symbol,:Integer,:Rational,:Complex]
    return string(expr)
  else
    error("Excpetion: $(expr)")
  end # if

  return string(expr)

end # function gen_mma_str





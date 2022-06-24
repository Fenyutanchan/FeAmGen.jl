
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
gen_sorted_str( ::Val{:Symbol}, expr::Basic )::String = string(expr)
gen_sorted_str( ::Val{:Integer}, expr::Basic )::String = string(expr)
gen_sorted_str( ::Val{:Rational}, expr::Basic )::String = string(expr)
gen_sorted_str( ::Val{:Complex}, expr::Basic )::String = string(expr)
###########################################


###########################################
function gen_sorted_str( ::Val{:Pow}, expr::Basic )::String
###########################################
  arglist = get_args(expr)
  return ":Pow( $(gen_sorted_str(arglist[1])), $(arglist[2]) )"
end # function gen_sorted_str

###########################################
function gen_sorted_str( ::Val{:FunctionSymbol}, expr::Basic )::String
###########################################
  name = replace( get_name(expr), "Trace" => "DiracTrace" )
  return "$(name)( $(join( map( gen_sorted_str, get_args(expr) ), "," )) )"
end # function gen_sorted_str

###########################################
function gen_sorted_str( ::Val{:Mul}, expr::Basic )::String
###########################################
  return ":Mul( $(join( (sort∘map)( gen_sorted_str, get_args(expr) ), "," )) )"
end # function gen_sorted_str

###########################################
function gen_sorted_str( ::Val{:Add}, expr::Basic )::String
###########################################
  return ":Add( $(join( (sort∘map)( gen_sorted_str, get_args(expr) ), "," )) )"
end # function gen_sorted_str


###########################################
"""
    gen_sorted_str( expr::Basic )::String

This is a generic interface for different classes of the `expr`.
And it will generate the sorted string format for the expression `expr`.
"""
function gen_sorted_str( expr::Basic )::String
###########################################

  expr_class = SymEngine.get_symengine_class(expr)
  return gen_sorted_str( Val(expr_class), expr )

end # function gen_sorted_str




#####################################################
function get_add_vector( expr::Basic )::Vector{Basic}
#####################################################

  return SymEngine.get_symengine_class( expr ) == :Add ? get_args(expr) : Basic[ expr ]

end # function get_add_vector



#########################################
function drop_coeff( expr::Basic )::Basic
#########################################

  return SymEngine.get_symengine_class( expr ) == :Mul ?
         (prod∘filter)( x_-> SymEngine.get_symengine_class(x_) != :Integer, get_args(expr) ) : expr

end # function drop_coeff


#########################################
function generate_SPcombo( 
    rank_str::String, 
    independent_mom_list::Array{Basic}
)::Array{Basic}
#########################################

  # e.g. generate_SP([2,2],[k1,k2,K3])
  # only for two-loop level

  @vars q1,q2
  @funs SP

  ###### r1 is the rank of q1, r2 is the rank of q2
  r1 = count( "q1", rank_str )
  r2 = count( "q2", rank_str )
  num_mom = length(independent_mom_list)

  # term_q1_k = SP(q1,k1)+SP(q1,k2)+...
  term_q1_k = 0
  for i in 1:num_mom
    term_q1_k = SP(q1,independent_mom_list[i])+term_q1_k
  end # for

  # term_q2_k = SP(q2,k1)+SP(q2,k2)+...
  term_q2_k = 0
  for i in 1:num_mom
    term_q2_k = SP(q2,independent_mom_list[i])+term_q2_k
  end # for

  #q1_q2_exp is the exponent of SP(q1,q2)
  #q1_q1_exp is the exponent of SP(q1,q1)
  #q2_q2_exp is the exponent of SP(q2,q2)
  #q1_remain = total q1 rank - q1 rank in SP(q1,q2)
  #q2_remain = total q2 rank - q2 rank in SP(q1,q2)
  #q1_k_exp is the exponent of (SP(q1,k1)+SP(q1,k2)+...)
  #q2_k_exp is the exponent of (SP(q2,k1)+SP(q2,k2)+...)
  #total_term is the summation of all possible terms
  #expand the total_term and drop the coefficients

  total_term = 0
  for q1_q2_exp in 0:min(r1,r2)
    q1_remain = r1-q1_q2_exp
    q2_remain = r2-q1_q2_exp
    max_q1_q1_exp = floor(Int,q1_remain/2)
    max_q2_q2_exp = floor(Int,q2_remain/2)
    for q1_q1_exp in 0:max_q1_q1_exp
      q1_k_exp = q1_remain - 2*q1_q1_exp
      for q2_q2_exp in 0:max_q2_q2_exp
        q2_k_exp = q2_remain - 2*q2_q2_exp
        total_term += term_q1_k^q1_k_exp*term_q2_k^q2_k_exp*
                      SP(q1,q1)^q1_q1_exp*SP(q2,q2)^q2_q2_exp*SP(q1,q2)^q1_q2_exp
      end # for q2_q2_exp
    end # for q1_q1_exp
  end # for q1_q2_exp

  SP_list = (get_add_vector∘expand)(total_term)
  SP_list = map(drop_coeff,SP_list)

  return SP_list

end # function generate_SPcombo






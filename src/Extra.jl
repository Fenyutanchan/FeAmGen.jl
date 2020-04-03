
#################################################################
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

@testset "degree" begin
  @vars x, y
  poly = 1+2*x+y*3*x^3
  @test degree( poly, x ) == 3
  @test degree( poly, x, max_n=3 ) == 3
end # @testset




###################################################
function make_SP( mom1::Basic, mom2::Basic )::Basic
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


@testset "SP" begin
  @vars k1, k2, k3, k4
  @funs SP
  @test make_SP(k1,k2) == SP(k1,k2)
  @test make_SP(k1,-k2) == -SP(k1,k2)
  @test make_SP(-k1,-k2) == SP(k1,k2)
  @test make_SP(k1-k3,k2+k4) == SP(k1, k2) + SP(k1, k4) - SP(k2, k3) - SP(k3, k4)
end # @testset


###########################################
function gen_mma_str( expr::Basic )::String
###########################################

  expr_class = SymEngine.get_symengine_class(expr)

  if expr_class == :Add
    return "Plus[ $(join( (sort∘map)( gen_mma_str, get_args(expr) ), "," )) ]"
  elseif expr_class == :Mul
    return "Times[ $(join( (sort∘map)( gen_mma_str, get_args(expr) ), "," )) ]"
  elseif expr_class == :FunctionSymbol
    return "$(get_name(expr))[ $(join( map( gen_mma_str, get_args(expr) ), "," )) ]"
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






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
  exit()
end # function degree

@testset "degree" begin
  @vars x, y
  poly = 1+2*x+y*3*x^3
  @test degree( poly, x ) == 3
  @test degree( poly, x, max_n=3 ) == 3
end # @testset



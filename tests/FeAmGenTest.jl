using SymEngine, FeAmGen, Test, BenchmarkTools

@testset "degree" begin
  @vars x, y
  poly = 1+2*x+y*3*x^3
  @test FeAmGen.degree( poly, x ) == 3
  @test FeAmGen.degree( poly, x, max_n=3 ) == 3
end # @testset



@testset "make_SP" begin
  @funs SP
  @vars k1, k2, k3, k4, q1, q2
  @test FeAmGen.make_SP(k1,k2) == SP(k1,k2)
  @test FeAmGen.make_SP(k1,-k2) == -SP(k1,k2)
  @test FeAmGen.make_SP(-k1,-k2) == SP(k1,k2)
  @test FeAmGen.make_SP(k1-k3,k2+k4) == SP(k1, k2) + SP(k1, k4) - SP(k3, k2) - SP(k3, k4)

  @test FeAmGen.make_SP( k1, k2 ) == SP(k1, k2)
  @test FeAmGen.make_SP( k1, k2+2*k3 ) == SP(k1, k2) + 2*SP(k1, k3)
  test_mom = 2*k1+3*k2+k3+5*k4+6*q1+7*q2
  @test FeAmGen.make_SP( test_mom, test_mom ) == 4*SP(k1, k1) + 6*SP(k1, k2) + 2*SP(k1, k3) + 10*SP(k1, k4) + 12*SP(k1, q1) + 14*SP(k1, q2) + 6*SP(k2, k1) + 9*SP(k2, k2) + 3*SP(k2, k3) + 15*SP(k2, k4) + 18*SP(k2, q1) + 21*SP(k2, q2) + 2*SP(k3, k1) + 3*SP(k3, k2) + SP(k3, k3) + 5*SP(k3, k4) + 6*SP(k3, q1) + 7*SP(k3, q2) + 10*SP(k4, k1) + 15*SP(k4, k2) + 5*SP(k4, k3) + 25*SP(k4, k4) + 30*SP(k4, q1) + 35*SP(k4, q2) + 12*SP(q1, k1) + 18*SP(q1, k2) + 6*SP(q1, k3) + 30*SP(q1, k4) + 36*SP(q1, q1) + 42*SP(q1, q2) + 14*SP(q2, k1) + 21*SP(q2, k2) + 7*SP(q2, k3) + 35*SP(q2, k4) + 42*SP(q2, q1) + 49*SP(q2, q2)
end # @testset



@vars k1, k2, k3, k4, q1, q2
test_mom = 2*k1+3*k2+k3+5*k4+6*q1+7*q2
expr = expand(test_mom^3+test_mom)

println( "FeAmGen.gen_mma_str(expr)" )
@btime FeAmGen.gen_mma_str(expr) 
println( "FeAmGen.gen_mma_str_old(expr)" )
@btime FeAmGen.gen_mma_str_old(expr) 

exit()

@vars k1, k2, k3, k4, q1, q2
aa = FeAmGen.convert_to_array( k1 )
println( aa )
aa = FeAmGen.convert_to_array( 2*k1 )
println( aa )
aa = FeAmGen.convert_to_array( 2*k1+3*k2+k3 )
println( aa )

#test_mom = 2*k1+3*k2+k3+5*k4+6*q1+7*q2
#aa = FeAmGen.make_SP_new( test_mom, test_mom )
#println( aa )

#exit()

test_mom = 2*k1+3*k2+k3+5*k4+6*q1+7*q2

println( "make_SP_new: " )
@btime aa = FeAmGen.make_SP_new( test_mom, test_mom )

println( "make_SP: " )
@btime aa = FeAmGen.make_SP( test_mom, test_mom )


@testset "tensor_product" begin
  target_in = [ ["a","b","c"], ["d","e","f"] ]
  target_out = ["a,d","a,e","a,f","b,d","b,e","b,f","c,d","c,e","c,f"]
  @test tensor_product( target_in... ) == target_out
  target_in = [ ["a","b"], ["c","d"], ["e","f"] ]
  target_out = ["a,c,e","a,c,f","a,d,e","a,d,f","b,c,e","b,c,f","b,d,e","b,d,f"]
  @test tensor_product( target_in... ) == target_out
end # @testset


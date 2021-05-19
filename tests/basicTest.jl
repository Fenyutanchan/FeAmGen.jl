using SymEngine, FeAmGen, Test, BenchmarkTools, YAML, JLD, Pipe, Logging, Dates

io = open("basicTest.log", "w+")
logger = SimpleLogger(io)
global_logger(logger)

@info "basicTest starts @ $(now())"

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
  @test FeAmGen.make_SP(k1-k3,k2+k4) == SP(k1, k2) + SP(k1, k4) - SP(k2, k3) - SP(k3, k4)

  @test FeAmGen.make_SP( k1, k2 ) == SP(k1, k2)
  @test FeAmGen.make_SP( k1, k2+2*k3 ) == SP(k1, k2) + 2*SP(k1, k3)
  test_mom = 2*k1+3*k2+k3+5*k4+6*q1+7*q2
  @test FeAmGen.make_SP( test_mom, test_mom ) == 4*SP(k1, k1) + 12*SP(k1, k2) + 4*SP(k1, k3) + 20*SP(k1, k4) + 24*SP(k1, q1) + 28*SP(k1, q2) + 9*SP(k2, k2) + 6*SP(k2, k3) + 30*SP(k2, k4) + 36*SP(k2, q1) + 42*SP(k2, q2) + SP(k3, k3) + 10*SP(k3, k4) + 12*SP(k3, q1) + 14*SP(k3, q2) + 25*SP(k4, k4) + 60*SP(k4, q1) + 70*SP(k4, q2) + 36*SP(q1, q1) + 84*SP(q1, q2) + 49*SP(q2, q2)
end # @testset



@testset "gen_mma_str" begin
  @vars k1, k2, k3, k4, q1, q2
  test_mom = 2*k1+3*k2+k3+5*k4+6*q1+7*q2
  expr = expand(test_mom^3+test_mom)
  @test FeAmGen.gen_mma_str(expr) == "Plus[ Power[ k3, 3 ],Times[ 108,Power[ q1, 2 ],k3 ],Times[ 108,k2,k3,q1 ],Times[ 12,Power[ k1, 2 ],k3 ],Times[ 125,Power[ k4, 3 ] ],Times[ 126,k2,k3,q2 ],Times[ 1260,k4,q1,q2 ],Times[ 135,Power[ k2, 2 ],k4 ],Times[ 147,Power[ q2, 2 ],k3 ],Times[ 15,Power[ k3, 2 ],k4 ],Times[ 150,Power[ k4, 2 ],k1 ],Times[ 162,Power[ k2, 2 ],q1 ],Times[ 18,Power[ k3, 2 ],q1 ],Times[ 180,k1,k2,k4 ],Times[ 180,k3,k4,q1 ],Times[ 189,Power[ k2, 2 ],q2 ],Times[ 2,k1 ],Times[ 21,Power[ k3, 2 ],q2 ],Times[ 210,k3,k4,q2 ],Times[ 216,Power[ q1, 2 ],k1 ],Times[ 216,Power[ q1, 3 ] ],Times[ 216,k1,k2,q1 ],Times[ 225,Power[ k4, 2 ],k2 ],Times[ 252,k1,k2,q2 ],Times[ 252,k3,q1,q2 ],Times[ 27,Power[ k2, 2 ],k3 ],Times[ 27,Power[ k2, 3 ] ],Times[ 294,Power[ q2, 2 ],k1 ],Times[ 3,k2 ],Times[ 324,Power[ q1, 2 ],k2 ],Times[ 343,Power[ q2, 3 ] ],Times[ 36,Power[ k1, 2 ],k2 ],Times[ 36,k1,k2,k3 ],Times[ 360,k1,k4,q1 ],Times[ 420,k1,k4,q2 ],Times[ 441,Power[ q2, 2 ],k2 ],Times[ 450,Power[ k4, 2 ],q1 ],Times[ 5,k4 ],Times[ 504,k1,q1,q2 ],Times[ 525,Power[ k4, 2 ],q2 ],Times[ 54,Power[ k2, 2 ],k1 ],Times[ 540,Power[ q1, 2 ],k4 ],Times[ 540,k2,k4,q1 ],Times[ 6,Power[ k3, 2 ],k1 ],Times[ 6,q1 ],Times[ 60,Power[ k1, 2 ],k4 ],Times[ 60,k1,k3,k4 ],Times[ 630,k2,k4,q2 ],Times[ 7,q2 ],Times[ 72,Power[ k1, 2 ],q1 ],Times[ 72,k1,k3,q1 ],Times[ 735,Power[ q2, 2 ],k4 ],Times[ 75,Power[ k4, 2 ],k3 ],Times[ 756,Power[ q1, 2 ],q2 ],Times[ 756,k2,q1,q2 ],Times[ 8,Power[ k1, 3 ] ],Times[ 84,Power[ k1, 2 ],q2 ],Times[ 84,k1,k3,q2 ],Times[ 882,Power[ q2, 2 ],q1 ],Times[ 9,Power[ k3, 2 ],k2 ],Times[ 90,k2,k3,k4 ],k3 ]" 
end # @testset


@testset "tensor_product" begin
  target_in = Array{String}[ String["a","b","c"], String["d","e","f"] ]
  target_out = String["a,d","a,e","a,f","b,d","b,e","b,f","c,d","c,e","c,f"]
  @test FeAmGen.tensor_product( target_in... ) == target_out
  target_in = Array{String}[ String["a","b"], String["c","d"], String["e","f"] ]
  target_out = String["a,c,e","a,c,f","a,d,e","a,d,f","b,c,e","b,c,f","b,d,e","b,d,f"]
  @test FeAmGen.tensor_product( target_in... ) == target_out
end # @testset


@testset "expand_parton" begin
  @test FeAmGen.expand_parton( String["u","parton"], String["d","g"], String["u","d","g"] ) == String["u,u,d,g", "u,d,d,g", "u,g,d,g"]
end # @testset


@testset "sort_proc_str" begin
  @test FeAmGen.sort_proc_str( "u,d,g,d,u", 2 ) == "d,u,d,g,u"
end # @testset


@testset "get_list_quoted_str" begin
  @test FeAmGen.get_list_quoted_str( String["a","b","c"] ) == "[ \"a\",\"b\",\"c\" ]"
end # @testset

@info "basicTest ends @ $(now())"

close(io)


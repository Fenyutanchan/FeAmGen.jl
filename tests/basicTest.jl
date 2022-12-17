using SymEngine, FeAmGen, Test, Dates

@info "basicTest starts @ $(now())"

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


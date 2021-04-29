using SymEngine, FeAmGen, Test, BenchmarkTools, JLD, Pipe

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








#----------------------------------------------------------------------------
# gg->ttbar 0-loop, 1-loop tests
#----------------------------------------------------------------------------
generic_ggttbar_seed_proc_yaml_str( ; nloop::Int64 = 2::Int64 ) = """
# input file for calculation details
# model related information

# model name
model_name: "sm"

# use unitary_gauge for internal vector boson
unitary_gauge: false

# content of "quark-parton", the "parton" will also contain gluon. Only for seed program.
#const partons = Array{String,1}( [ "g", "u", "d", "ubar", "dbar", "s", "c", "b", "sbar", "cbar", "bbar" ] )   
# for single top @ NNLO
#const partons = Array{String,1}( [ "g", "u", "d", "ubar", "dbar", "b", "bbar" ] )   
# for Higgs+Jet @ NLO
partons: [ "g", "u", "ubar", "d", "dbar", "b", "bbar" ] 

# only for seed program
AllowLeptonNumberViolation: false
AllowQuarkGenerationViolation: false

# process information
DropTadpole: true              # drop tadpole?
DropWFcorrection: true         # drop WFcorrection?

# number of loops
n_loop: $(nloop)
# order of QCD counter-term vertices
QCDCT_order: 0   

# order of QCD coupling gs in the amplitude
Amp_QCD_order: $(2+2*nloop) 
# order of QED coupling ee in the amplitude
Amp_QED_order: 0  

# min eps power in the amplitude
Amp_Min_Eps_Xpt: $(-2*nloop)
# max eps power in the amplitude
Amp_Max_Eps_Xpt: 0

# incoming and outgoing information
incoming: [ "g", "g" ]          # incoming particles
outgoing: [ "t", "tbar" ]               # outgoing particles 
"""

for nloop in [0,1]

  open( "ggttbar_seed_proc_$(nloop)Loop.yaml", "w" ) do infile
    write( infile, generic_ggttbar_seed_proc_yaml_str(nloop=nloop) )
  end 

  digest_seed_proc( "ggttbar_seed_proc_$(nloop)Loop.yaml" )

  generate_amp( "g_g_TO_t_tbar_$(nloop)Loop/g_g_TO_t_tbar.yaml" )

end # for nloop

@testset "gg->ttbar" for nloop in [0,1]

  n_diagram = @pipe readdir( "g_g_TO_t_tbar_$(nloop)Loop_amplitudes" ) |> filter( name->name[end-3:end]==".jld", _ ) |> length

  @testset "$(nloop)-loop diagrams" for diagram_index in 1:n_diagram

    content_dict = load( "g_g_TO_t_tbar_$(nloop)Loop_amplitudes/amplitude_diagram$(diagram_index).jld" )
    content_dict_bench = load( "g_g_TO_t_tbar_$(nloop)Loop_amplitudes_benchmark/amplitude_diagram$(diagram_index).jld" )

    @test content_dict == content_dict_bench 

    visual_bench_file = open( "g_g_TO_t_tbar_$(nloop)Loop_visuals_benchmark/visual_diagram$(diagram_index).tex" ) 
    visual_list_bench = readlines( visual_bench_file )
    close( visual_bench_file )
    
    visual_file = open( "g_g_TO_t_tbar_$(nloop)Loop_visuals/visual_diagram$(diagram_index).tex" )
    visual_list = readlines( visual_file )
    close( visual_file )

    @test visual_list == visual_list_bench 
  end # testset for diagram_index

end # testset







#----------------------------------------------------------------------------
# ee->HZ 0-loop, 1-loop tests
#----------------------------------------------------------------------------
generic_eeHZ_seed_proc_yaml_str( ; nloop::Int64 = 2::Int64 ) = """
# input file for calculation details
# model related information

# model name
model_name: "sm"

# use unitary_gauge for internal vector boson
unitary_gauge: false

# content of "quark-parton", the "parton" will also contain gluon. Only for seed program.
#const partons = Array{String,1}( [ "g", "u", "d", "ubar", "dbar", "s", "c", "b", "sbar", "cbar", "bbar" ] )   
# for single top @ NNLO
#const partons = Array{String,1}( [ "g", "u", "d", "ubar", "dbar", "b", "bbar" ] )   
# for Higgs+Jet @ NLO
partons: [ "g", "u", "ubar", "d", "dbar", "b", "bbar" ] 

# only for seed program
AllowLeptonNumberViolation: false
AllowQuarkGenerationViolation: false

# process information
DropTadpole: true              # drop tadpole?
DropWFcorrection: true         # drop WFcorrection?

# number of loops
n_loop: $(nloop)
# order of QCD counter-term vertices
QCDCT_order: 0   

# order of QCD coupling gs in the amplitude
Amp_QCD_order: 0 
# order of QED coupling ee in the amplitude
Amp_QED_order: $(2+2*nloop)

# min eps power in the amplitude
Amp_Min_Eps_Xpt: $(-2*nloop)
# max eps power in the amplitude
Amp_Max_Eps_Xpt: 0

# incoming and outgoing information
incoming: [ "eplus", "eminus" ]          # incoming particles
outgoing: [ "H", "Z" ]               # outgoing particles 
"""

for nloop in [0,1]

  open( "eeHZ_seed_proc_$(nloop)Loop.yaml", "w" ) do infile
    write( infile, generic_eeHZ_seed_proc_yaml_str(nloop=nloop) )
  end 

  digest_seed_proc( "eeHZ_seed_proc_$(nloop)Loop.yaml" )

  generate_amp( "eplus_eminus_TO_H_Z_$(nloop)Loop/eminus_eplus_TO_H_Z.yaml" )

end # for nloop

@testset "ee->HZ" for nloop in [0,1]

  n_diagram = @pipe readdir( "eminus_eplus_TO_H_Z_$(nloop)Loop_amplitudes" ) |> filter( name->name[end-3:end]==".jld", _ ) |> length

  @testset "$(nloop)-loop diagrams" for diagram_index in 1:n_diagram

    content_dict = load( "eminus_eplus_TO_H_Z_$(nloop)Loop_amplitudes/amplitude_diagram$(diagram_index).jld" )
    content_dict_bench = load( "eminus_eplus_TO_H_Z_$(nloop)Loop_amplitudes_benchmark/amplitude_diagram$(diagram_index).jld" )

    @test content_dict == content_dict_bench 

    visual_bench_file = open( "eminus_eplus_TO_H_Z_$(nloop)Loop_visuals_benchmark/visual_diagram$(diagram_index).tex" ) 
    visual_list_bench = readlines( visual_bench_file )
    close( visual_bench_file )
    
    visual_file = open( "eminus_eplus_TO_H_Z_$(nloop)Loop_visuals/visual_diagram$(diagram_index).tex" )
    visual_list = readlines( visual_file )
    close( visual_file )

    @test visual_list == visual_list_bench 
  end # testset for diagram_index

end # testset








#----------------------------------------------------------------------------
# Drell-Yan 0-loop, 1-loop, 2-loop tests
#----------------------------------------------------------------------------
generic_dy_seed_proc_yaml_str( ; nloop::Int64 = 2::Int64 ) = """
# input file for calculation details
# model related information

# model name
model_name: "sm_CKMdiag_Haa"

# use unitary_gauge for internal vector boson
unitary_gauge: true

# content of "quark-parton", the "parton" will also contain gluon. Only for seed program.
#const partons = Array{String,1}( [ "g", "u", "d", "ubar", "dbar", "s", "c", "b", "sbar", "cbar", "bbar" ] )   
# for single top @ NNLO
#const partons = Array{String,1}( [ "g", "u", "d", "ubar", "dbar", "b", "bbar" ] )   
# for Higgs+Jet @ NLO
partons: [ "g", "u", "ubar", "d", "dbar", "b", "bbar" ] 

# only for seed program
AllowLeptonNumberViolation: false
AllowQuarkGenerationViolation: false

# process information
DropTadpole: true              # drop tadpole?
DropWFcorrection: true         # drop WFcorrection?

# number of loops
n_loop: $nloop    
# order of QCD counter-term vertices
QCDCT_order: 0   

# order of QCD coupling gs in the amplitude
Amp_QCD_order: $(2*nloop) 
# order of QED coupling ee in the amplitude
Amp_QED_order: 1 

# min eps power in the amplitude
Amp_Min_Eps_Xpt: $(-2*nloop)
# max eps power in the amplitude
Amp_Max_Eps_Xpt: 0

# incoming and outgoing information
incoming: [ "dbar", "u" ]          # incoming particles
outgoing: [ "Wplus" ]               # outgoing particles 
"""



for nloop in [0,1,2]

  open( "dy_seed_proc_$(nloop)Loop.yaml", "w" ) do infile
    write( infile, generic_dy_seed_proc_yaml_str(nloop=nloop) )
  end 

  digest_seed_proc( "dy_seed_proc_$(nloop)Loop.yaml" )

  generate_amp( "dbar_u_TO_Wplus_$(nloop)Loop/dbar_u_TO_Wplus.yaml" )

end # for nloop

@testset "Drell-Yan" for nloop in [0,1,2]

  n_diagram = @pipe readdir( "dbar_u_TO_Wplus_$(nloop)Loop_amplitudes" ) |> filter( name->name[end-3:end]==".jld", _ ) |> length

  @testset "$(nloop)-loop diagrams" for diagram_index in 1:n_diagram

    content_dict = load( "dbar_u_TO_Wplus_$(nloop)Loop_amplitudes/amplitude_diagram$(diagram_index).jld" )
    content_dict_bench = load( "dbar_u_TO_Wplus_$(nloop)Loop_amplitudes_benchmark/amplitude_diagram$(diagram_index).jld" )

    @test content_dict == content_dict_bench 

    visual_bench_file = open( "dbar_u_TO_Wplus_$(nloop)Loop_visuals_benchmark/visual_diagram$(diagram_index).tex" ) 
    visual_list_bench = readlines( visual_bench_file )
    close( visual_bench_file )
    
    visual_file = open( "dbar_u_TO_Wplus_$(nloop)Loop_visuals/visual_diagram$(diagram_index).tex" )
    visual_list = readlines( visual_file )
    close( visual_file )

    @test visual_list == visual_list_bench 
  end # testset for diagram_index

end # testset






#----------------------------------------------------------------------------
# single-top gb->tW 0-loop, 1-loop, 2-loop tests
#----------------------------------------------------------------------------
generic_gbtw_seed_proc_yaml_str( ; nloop::Int64 = 2::Int64 ) = """
# input file for calculation details
# model related information

# model name
model_name: "sm_CKMdiag_Haa"

# use unitary_gauge for internal vector boson
unitary_gauge: true

# content of "quark-parton", the "parton" will also contain gluon. Only for seed program.
#const partons = Array{String,1}( [ "g", "u", "d", "ubar", "dbar", "s", "c", "b", "sbar", "cbar", "bbar" ] )   
# for single top @ NNLO
#const partons = Array{String,1}( [ "g", "u", "d", "ubar", "dbar", "b", "bbar" ] )   
# for Higgs+Jet @ NLO
partons: [ "g", "u", "ubar", "d", "dbar", "b", "bbar" ] 

# only for seed program
AllowLeptonNumberViolation: false
AllowQuarkGenerationViolation: false

# process information
DropTadpole: true              # drop tadpole?
DropWFcorrection: true         # drop WFcorrection?

# number of loops
n_loop: $(nloop)    
# order of QCD counter-term vertices
QCDCT_order: 0   

# order of QCD coupling gs in the amplitude
Amp_QCD_order: $(1+2*nloop) 
# order of QED coupling ee in the amplitude
Amp_QED_order: 1 

# min eps power in the amplitude
Amp_Min_Eps_Xpt: $(-2*nloop)
# max eps power in the amplitude
Amp_Max_Eps_Xpt: 0

# incoming and outgoing information
incoming: [ "g", "b" ]          # incoming particles
outgoing: [ "t", "Wminus" ]               # outgoing particles 
"""

for nloop in [0,1,2]

  open( "gbtw_seed_proc_$(nloop)Loop.yaml", "w" ) do infile
    write( infile, generic_gbtw_seed_proc_yaml_str(nloop=nloop) )
  end 

  digest_seed_proc( "gbtw_seed_proc_$(nloop)Loop.yaml" )

  generate_amp( "g_b_TO_t_Wminus_$(nloop)Loop/b_g_TO_Wminus_t.yaml" )

end # for nloop

@testset "gb->tW" for nloop in [0,1,2]

  n_diagram = @pipe readdir( "b_g_TO_Wminus_t_$(nloop)Loop_amplitudes" ) |> filter( name->name[end-3:end]==".jld", _ ) |> length

  @testset "$(nloop)-loop diagrams" for diagram_index in 1:n_diagram

    content_dict = load( "b_g_TO_Wminus_t_$(nloop)Loop_amplitudes/amplitude_diagram$(diagram_index).jld" )
    content_dict_bench = load( "b_g_TO_Wminus_t_$(nloop)Loop_amplitudes_benchmark/amplitude_diagram$(diagram_index).jld" )

    @test content_dict == content_dict_bench 

    visual_bench_file = open( "b_g_TO_Wminus_t_$(nloop)Loop_visuals_benchmark/visual_diagram$(diagram_index).tex" ) 
    visual_list_bench = readlines( visual_bench_file )
    close( visual_bench_file )
    
    visual_file = open( "b_g_TO_Wminus_t_$(nloop)Loop_visuals/visual_diagram$(diagram_index).tex" )
    visual_list = readlines( visual_file )
    close( visual_file )

    @test visual_list == visual_list_bench 
  end # testset for diagram_index

end # testset



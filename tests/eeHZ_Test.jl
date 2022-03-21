using SymEngine, FeAmGen, Test, BenchmarkTools, YAML, JLD, Pipe, Dates, Logging

io = open("eeHZ_Test.log", "w+")
logger = SimpleLogger(io)
global_logger(logger)

@info "eeHZ_Test starts @ $(now())"




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

# min ep power in the amplitude
Amp_Min_Ep_Xpt: $(-2*nloop)
# max ep power in the amplitude
Amp_Max_Ep_Xpt: 0

# incoming and outgoing information
incoming: [ "eplus", "eminus" ]          # incoming particles
outgoing: [ "H", "Z" ]               # outgoing particles 
"""

for nloop in [0,1,2]

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



@info "eeHZ_Test ends @ $(now())"

close(io)


using SymEngine, FeAmGen, Test, YAML, JLD, Pipe, Dates

@info "DrellYan_Test starts @ $(now())"

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
# order of special coupling in the amplitude
Amp_SPC_order: 0 

# min ep power in the amplitude
Amp_Min_Ep_Xpt: $(-2*nloop)
# max ep power in the amplitude
Amp_Max_Ep_Xpt: 0

# incoming and outgoing information
incoming: [ "dbar", "u" ]          # incoming particles
outgoing: [ "Wplus" ]               # outgoing particles 

# whether to check the consistency between two versions of amplitudes
check_consistency: true

"""



for nloop in [0,1,2]

  open( "dy_seed_proc_$(nloop)Loop.yaml", "w" ) do infile
    write( infile, generic_dy_seed_proc_yaml_str(nloop=nloop) )
  end # close

  digest_seed_proc( "dy_seed_proc_$(nloop)Loop.yaml" )

  generate_amp( "dbar_u_TO_Wplus_$(nloop)Loop/dbar_u_TO_Wplus.yaml" )

end # for nloop

for nloop in [0,1,2]

  n_diagram = @pipe readdir( "dbar_u_TO_Wplus_$(nloop)Loop_amplitudes" ) |> filter( name->name[end-4:end]==".jld2", _ ) |> length

  @testset "Drell-Yan $(nloop)-loop diagrams" begin
  for diagram_index in 1:n_diagram

    content_dict = load( "dbar_u_TO_Wplus_$(nloop)Loop_amplitudes/amplitude_diagram$(diagram_index).jld2" )
    content_dict_bench = load( "dbar_u_TO_Wplus_$(nloop)Loop_amplitudes_benchmark/amplitude_diagram$(diagram_index).jld2" )

    @test content_dict == content_dict_bench 

    visual_bench_file = open( "dbar_u_TO_Wplus_$(nloop)Loop_visuals_benchmark/visual_diagram$(diagram_index).tex" ) 
    visual_list_bench = readlines( visual_bench_file )
    close( visual_bench_file )
    
    visual_file = open( "dbar_u_TO_Wplus_$(nloop)Loop_visuals/visual_diagram$(diagram_index).tex" )
    visual_list = readlines( visual_file )
    close( visual_file )

    @test visual_list == visual_list_bench 
  end # for diagram_index
  end # testset 

end # for nloop

@info "DrellYan_Test ends @ $(now())"



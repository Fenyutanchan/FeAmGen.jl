using SymEngine, FeAmGen, Test, YAML, JLD2, Dates, Pipe

@info "ggttbar_Test starts @ $(now())"


##########################################################################
function generate_amp_debug( proc_file::String )::Nothing
######################################################################

  #------------------------------------------------------------------
  @assert isfile(proc_file) "The first argument is not a file!"
  input = YAML.load_file( proc_file )
  #------------------------------------------------------------------


  #------------------------------------------------------------------
  @info "Choose model" model=input["model_name"]
  model = FeAmGen.readin_model( input )
  FeAmGen.generate_QGRAF_model( model )
  FeAmGen.logging_model( model )
  #------------------------------------------------------------------

  @info "Usage of unitary gauge" unitary_gauge=input["unitary_gauge"]

  @info "Drop Tadpole" DropTadpole=input["DropTadpole"]

  @info "Drop WFcorrection" DropWFcorrection=input["DropWFcorrection"]

  @info "Number of loops" n_loop=input["n_loop"]

  @info "QCD CT-order" QCDCT_order=input["QCDCT_order"]

  @info "Order of QCD coupling gs in the amplitude" Amp_QCD_order=input["Amp_QCD_order"]

  @info "Order of QED coupling gs in the amplitude" Amp_QED_order=input["Amp_QED_order"]

  @info "Min ep power in the amplitude" Amp_Min_Ep_Xpt=input["Amp_Min_Ep_Xpt"]
  @info "Max ep power in the amplitude" Amp_Max_Ep_Xpt=input["Amp_Max_Ep_Xpt"]

  @info "Incoming" input["incoming"]
  @info "Outgoing" input["outgoing"]

  @info "Coupling factor" couplingfactor=input["couplingfactor"]

  #----------------------------------------------------------------------
  # Run the QGRAF
  FeAmGen.generate_Feynman_diagram( model, input )

  generate_amplitude_debug( model, input )

  return nothing

end # function generate_amp_debug




##########################################################################
function generate_amplitude_debug( model::FeAmGen.Model, input::Dict{Any,Any} )::Nothing
##########################################################################

  proc_str = join( [ input["incoming"]; "TO"; input["outgoing"]; "$(input["n_loop"])Loop" ], "_" )

  n_loop = input["n_loop"]
  couplingfactor = Basic(input["couplingfactor"]) 

  qgraf_out = YAML.load_file( "qgraf_out.dat" )

  qgraf_list = qgraf_out["FeynmanDiagrams"]

  #------------------------------------------------  
  # Convert qgraf to Graph
  graph_list = @pipe qgraf_list |> 
               map( q -> FeAmGen.convert_qgraf_TO_Graph( q, model ), _ ) |>
               convert( Array{FeAmGen.Graph,1}, _ ) |>
               filter( !isnothing, _ ) |>
               sort( _, by= g->g.property[:diagram_index] )
  #------------------------------------------------  


  #------------------------------------------------  
  # Generate Gauge choice
  gauge_choice = generate_gauge_choice_debug( graph_list )

  return nothing

end 

####################################################################
function generate_gauge_choice_debug( graph_list::Vector{FeAmGen.Graph} )::Dict{Basic,Basic}
####################################################################

  # Only the external fields are needed
  graph0 = first( graph_list )

  ext_edge_list = filter( e_ -> ( e_.property[:style]=="External" ), graph0.edge_list )

  null_ext_edge_list = filter( e_ -> ( FeAmGen.is_massless(e_.property[:particle]) ), ext_edge_list )
 # @assert length(null_ext_edge_list) >= 2

  gauge_choice = Dict{Basic,Basic}()
  push!( gauge_choice, null_ext_edge_list[1].property[:ref2_MOM] => null_ext_edge_list[2].property[:null_MOM] )

  not1st_ext_edge_list = filter( e_ -> ( e_.property[:mark] != null_ext_edge_list[1].property[:mark] ), ext_edge_list )

  for edge in not1st_ext_edge_list
    if FeAmGen.is_massive_fermion(edge.property[:particle])
      push!( gauge_choice, edge.property[:ref2_MOM] => Basic("barK$(edge.property[:mark])") )
    else
      push!( gauge_choice, edge.property[:ref2_MOM] => null_ext_edge_list[1].property[:null_MOM] )
    end # if
  end # for edge

@show gauge_choice
error("DEBUG")

  return gauge_choice

end # function generate_gauge_choice_debug




















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

# min ep power in the amplitude
Amp_Min_Ep_Xpt: $(-2*nloop)
# max ep power in the amplitude
Amp_Max_Ep_Xpt: 0

# incoming and outgoing information
incoming: [ "g", "g" ]          # incoming particles
outgoing: [ "t", "tbar" ]               # outgoing particles 

# whether to check the consistency between two versions of amplitudes
check_consistency: true

"""

for nloop in [0,1]

  open( "ggttbar_seed_proc_$(nloop)Loop.yaml", "w" ) do infile
    write( infile, generic_ggttbar_seed_proc_yaml_str(nloop=nloop) )
  end # close

  digest_seed_proc( "ggttbar_seed_proc_$(nloop)Loop.yaml" )

  generate_amp_debug( "g_g_TO_t_tbar_$(nloop)Loop/g_g_TO_t_tbar.yaml" )

end # for nloop

#@testset "gg->ttbar" begin
for nloop in [0,1]

  n_diagram = @pipe readdir( "g_g_TO_t_tbar_$(nloop)Loop_amplitudes" ) |> filter( name->name[end-4:end]==".jld2", _ ) |> length

  @testset "gg->ttbar $(nloop)-loop diagrams" begin
  for diagram_index in 1:n_diagram

    content_dict = load( "g_g_TO_t_tbar_$(nloop)Loop_amplitudes/amplitude_diagram$(diagram_index).jld2" )
    content_dict_bench = load( "g_g_TO_t_tbar_$(nloop)Loop_amplitudes_benchmark/amplitude_diagram$(diagram_index).jld2" )

    @test content_dict == content_dict_bench 

    visual_bench_file = open( "g_g_TO_t_tbar_$(nloop)Loop_visuals_benchmark/visual_diagram$(diagram_index).tex" ) 
    visual_list_bench = readlines( visual_bench_file )
    close( visual_bench_file )
    
    visual_file = open( "g_g_TO_t_tbar_$(nloop)Loop_visuals/visual_diagram$(diagram_index).tex" )
    visual_list = readlines( visual_file )
    close( visual_file )

    @test visual_list == visual_list_bench 
  end # for diagram_index
  end # testset for diagram_index

end # for nloop
#end # testset


@info "ggttbar_Test ends @ $(now())"



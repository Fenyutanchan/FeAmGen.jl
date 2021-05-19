
######################################################################
"""
    generate_amp( proc_file::String )::Nothing

Read-in the card for the specific process and produce the relevant amplitude.
This function is one of the front-end functions in FeAmGen.jl.
The directory of model files are supposed in `(dirname∘dirname∘pathof∘Base.moduleroot)(FeAmGen)*"/Models"`.
"""
function generate_amp( proc_file::String )::Nothing
######################################################################

  #------------------------------------------------------------------
  @assert isfile(proc_file) "The first argument is not a file!"
  input = YAML.load_file( proc_file )
  #------------------------------------------------------------------


  #------------------------------------------------------------------
  @info "Choose model: $(input["model_name"])"
  model = readin_model( input )
  generate_QGRAF_model( model )
  logging_model( model )
  #------------------------------------------------------------------

  @info "Usage of unitary gauge: $(input["unitary_gauge"])"

  @info "drop Tadpole: $(input["DropTadpole"])"

  @info "drop WFcorrection: $(input["DropWFcorrection"])"

  @info "# of loops: $(input["n_loop"])"

  @info "QCD CT-order: $(input["QCDCT_order"])"

  @info "order of QCD coupling gs in the amplitude: $(input["Amp_QCD_order"])"

  @info "order of QED coupling gs in the amplitude: $(input["Amp_QED_order"])"

  @info "min eps power in the amplitude: $(input["Amp_Min_Eps_Xpt"])"
  @info "max eps power in the amplitude: $(input["Amp_Max_Eps_Xpt"])"

  @info "incoming: $(input["incoming"])"
  @info "outgoing: $(input["outgoing"])"

  @info "coupling factor: $(input["couplingfactor"])"

  #----------------------------------------------------------------------
  # Run the QGRAF
  generate_Feynman_diagram( model, input )

  generate_amplitude( model, input )

  return nothing

end # function generate_amp





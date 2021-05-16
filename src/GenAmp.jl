
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
  green_message( "Choose model: ", string(input["model_name"]) )
  model = readin_model( input )
  generate_QGRAF_model( model )
  logging_model( model )
  #------------------------------------------------------------------

  green_message( "Usage of unitary gauge: ", string(input["unitary_gauge"]) )

  green_message( "drop Tadpole: ", string(input["DropTadpole"]) )

  green_message( "drop WFcorrection: ", string(input["DropWFcorrection"]) )

  green_message( "# of loops: ", string(input["n_loop"]) )

  green_message( "QCD CT-order: ", string(input["QCDCT_order"]) )

  green_message( "order of QCD coupling gs in the amplitude: ", string(input["Amp_QCD_order"]) )

  green_message( "order of QED coupling gs in the amplitude: ", string(input["Amp_QED_order"]) )

  green_message( "min eps power in the amplitude: ", string(input["Amp_Min_Eps_Xpt"]) )
  green_message( "max eps power in the amplitude: ", string(input["Amp_Max_Eps_Xpt"]) )

  green_message( "incoming: ", string(input["incoming"]) )
  green_message( "outgoing: ", string(input["outgoing"]) )

  green_message( "coupling factor: ", string(input["couplingfactor"]) )

  #----------------------------------------------------------------------
  # Run the QGRAF
  generate_Feynman_diagram( model, input )

  generate_amplitude( model, input )

  return nothing

end # function generate_amp





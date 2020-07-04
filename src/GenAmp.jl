
######################################################################
function generate_amp( proc_file::String, model_dir::String )::Nothing
######################################################################

  #------------------------------------------------------------------
  @assert isfile(proc_file) "The first argument is not a file!"
  file_stream = open(proc_file)
  input = YAML.load( file_stream )
  close( file_stream )
  #------------------------------------------------------------------


  #------------------------------------------------------------------
  green_message( "Choose model: ", string(input["model_name"]) )
  model = readin_model( input, model_dir )
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

  generate_Feynman_diagram( model, input )

  generate_amplitude( model, input )

  return nothing

end # function generate_amp





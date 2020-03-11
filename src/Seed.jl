
##########################################################################
function digest_seed_proc( seed_file::String, model_dir::String )::Nothing
##########################################################################



  #----------------------------------------------------------------------------------------------
  input = YAML.load( open(seed_file) )
  #----------------------------------------------------------------------------------------------


  #-------------------------------------------------------
  green_message( "Choose model: ", string(input["model_name"]) )
  particle_dict = simple_readin_model( input["model_name"], model_dir )

  printstyled( "All the particles: \n", color=:green )
  for name_part in particle_dict
    println( "  ", name_part )
  end # for part
  #-------------------------------------------------------


  #----------------------------------------------------------------------------------------------
  green_message( "unitary gauge: ", string(input["unitary_gauge"]) )
  green_message( "parton content: ", string(input["partons"]) )
  green_message( "Incoming: ", string(input["incoming"]) )
  green_message( "Outgoing: ", string(input["outgoing"]) )
  proc_list = expand_parton( input["incoming"], input["outgoing"], input["partons"] )
  #----------------------------------------------------------------------------------------------


  #----------------------------------------------------------------------------------------------
  n_inc = length(input["incoming"])
  proc_list = map( s_ -> filter_charge( s_, n_inc, particle_dict ), proc_list )

  if input["AllowLeptonNumberViolation"] == false
    proc_list = map( s_ -> filter_lepton_generations( s_, n_inc, particle_dict ), proc_list )
  end # if

  if input["AllowQuarkGenerationViolation"] == false
    proc_list = map( s_ -> filter_quark_generations( s_, n_inc, particle_dict ), proc_list )
  end # if

  proc_list = map( s_ -> sort_proc_str( s_, n_inc ), proc_list )

  proc_set = delete!( Set(proc_list), nothing )
  println()
  println( "Filtered subprocesses: " )
  for proc_str in proc_set 
    part_str_list = split( proc_str, "," )
    pretty_str = join( part_str_list[1:n_inc], "," )*" => "*join( part_str_list[n_inc+1:end], "," )
    printstyled( "  ", pretty_str, "\n", color=:green )
  end # for proc_str
  #----------------------------------------------------------------------------------------------


  #----------------------------------------------------------------------------------------------
  parton_proc_str = join( [ input["incoming"]; "TO"; input["outgoing"] ], "_" )
  if isdir( parton_proc_str ) == true
    rm( parton_proc_str, recursive=true )
  end # if
  mkdir( parton_proc_str )
  cd( parton_proc_str ) 

  print( "Writing subprocesses cards in " )
  printstyled( parton_proc_str, "\n", color=:green )
  for proc_str in proc_set 
    write_card( proc_str, n_inc, input )
  end # for proc_str
  println( "Done" )

  cd( ".." )
  #----------------------------------------------------------------------------------------------

  return nothing

end # function digest_seed_proc








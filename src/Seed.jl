
##########################################################################
"""
    digest_seed_proc( seed_file::String, model_dir::String )::Nothing

Read-in the card for seed process, and write-out the cards for the relevant specific processes.
This function is one of the front-end functions in FeAmGen.jl.
The directory of model files are supposed in ".".
"""
function digest_seed_proc( seed_file::String )::Nothing
##########################################################################



  #----------------------------------------------------------------------------------------------
  input = YAML.load_file( seed_file )
  #----------------------------------------------------------------------------------------------
  parton_proc_str = join( [ input["incoming"]; "TO"; input["outgoing"]; "$(input["n_loop"])Loop" ], "_" )

  box_message( "Read in model \"$(input["model_name"])\" for $(parton_proc_str)" )

  #-------------------------------------------------------
  @info "Choose model" model=input["model_name"]
  particle_dict = simple_readin_model( input["model_name"] )

  @info "All the particles" join( particle_dict, " " ) 
  #-------------------------------------------------------


  #----------------------------------------------------------------------------------------------
  @info "Unitary gauge" input["unitary_gauge"]
  @info "Parton content" input["partons"]
  @info "Incoming" input["incoming"]
  @info "Outgoing" input["outgoing"]
  proc_list = expand_parton( input["incoming"], input["outgoing"], input["partons"] )
  #----------------------------------------------------------------------------------------------


  #----------------------------------------------------------------------------------------------
  n_inc = length(input["incoming"])
  proc_list = map( x -> filter_charge( x, n_inc, particle_dict ), proc_list )

  if input["AllowLeptonNumberViolation"] == false
    proc_list = map( x -> filter_lepton_generations( x, n_inc, particle_dict ), proc_list )
  end # if

  if input["AllowQuarkGenerationViolation"] == false
    proc_list = map( x -> filter_quark_generations( x, n_inc, particle_dict ), proc_list )
  end # if

  proc_list = map( x -> sort_proc_str( x, n_inc ), proc_list )

  proc_set = delete!( Set(proc_list), nothing )
  @info "Filtered subprocesses" 
  for proc_str in proc_set 
    part_str_list = split( proc_str, "," )
    inc_str = join( part_str_list[1:n_inc], "," )
    out_str = join( part_str_list[n_inc+1:end], "," )
    pretty_str = "$(inc_str) => $(out_str)"
    @info "  "*pretty_str
  end # for proc_str
  #----------------------------------------------------------------------------------------------


  #----------------------------------------------------------------------------------------------
  if isdir( parton_proc_str ) == true
    rm( parton_proc_str, recursive=true )
  end # if
  mkdir( parton_proc_str )
  cd( parton_proc_str ) 

  @info "[ Generate subprocesses cards in $(parton_proc_str) ]"
  for proc_str in proc_set 
    write_card( proc_str, n_inc, input )
  end # for proc_str
  @info "Done" 

  cd( ".." )
  #----------------------------------------------------------------------------------------------

  return nothing

end # function digest_seed_proc








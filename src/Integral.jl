


#---------------------------------------------------
function generate_integral( yaml_file::String )::Nothing
#---------------------------------------------------

  file_stream = open( yaml_file )
  file_dict = YAML.load( file_stream ) 
  close( file_stream )

  name_str = file_dict["name"]
  n_loop = file_dict["n_loop"]::Int64
  min_eps_xpt = file_dict["min_eps_xpt"]::Int64
  max_eps_xpt = file_dict["max_eps_xpt"]::Int64
  ext_mom_list = file_dict["external_momenta"]::Vector{String}
  kin_relation = Dict{Basic,Basic}( map( p_->(Basic(p_[1]),Basic(p_[2])), file_dict["kin_relation"] ) )
  loop_den_list = map( Basic, file_dict["den_list"] )
  loop_den_xpt_list = file_dict["den_xpt_list"]

  #-----------------------
  ver_mass_list = free_symbols( vcat( collect( values(kin_relation) ), map(d_->get_args(d_)[2],loop_den_list) ) )

  scale2_list = Vector{Basic}( undef, length(ver_mass_list) )
  for index in 1:length(ver_mass_list)
    ver_mass_str = string(ver_mass_list[index])
    if ver_mass_str[1] == 'm'
      scale2_list[index] = ver_mass_list[index]^2
    else
      scale2_list[index] = ver_mass_list[index]
    end # if
  end # for index


  positive_pos_list = findall( x_ -> x_ > 0, loop_den_xpt_list )
  positive_loop_den_list = loop_den_list[positive_pos_list]
  positive_loop_den_xpt_list = loop_den_xpt_list[positive_pos_list]

  negative_pos_list = findall( x_ -> x_ < 0, loop_den_xpt_list )
  irreducible_numerator = one(Basic)
  for negative_pos in negative_pos_list
    negative_loop_den = loop_den_list[negative_pos]
    negative_xpt = loop_den_xpt_list[negative_pos]

    den_arg_list = get_args(negative_loop_den)
    irreducible_numerator *= ( make_SP(den_arg_list[1],den_arg_list[1]) - den_arg_list[2]^2 )^(-negative_xpt) 
  end # for negative_pos
  # need to further expand SP^n into (FV*FV)^n with different dummy indices, maybe use FORM script.

  numerator_expr = irreducible_numerator * Basic( file_dict["numerator"] )

  file_name = "numerator_contraction"
  file = open( "$(file_name).frm", "w" )
  write( file, make_amp_contraction_script( numerator_expr, file_name ) )
  close( file )

  file = open( "kin_relation.frm", "w" )
  write( file, (join∘map)( ele_->"id $(ele_[1]) = $(ele_[2]);\n", collect(kin_relation) ) )
  close(file)

  file = open( "model_parameters.frm", "w" )
  write( file, "symbol "*join( map( k_->string(k_), ver_mass_list ), "," )*";\n" )
  close(file)

  file = open( "contractor.frm", "w" )
  write( file, make_contractor_script() )
  close(file)

  file = open( "color.frm", "w" )
  write( file, make_color_script() )
  close(file)

  # baseINC only needs information from the external fields.
  touch( "baseINC.frm" )

  printstyled( "[ form $(file_name).frm ]\n", color=:yellow )
  run( pipeline( `form $(file_name).frm`, file_name*".log" ) )

  file = open( file_name*".out", "r" )
  result_str = read( file, String ) 
  close( file )

  # remove intermediate files
  rm( "baseINC.frm" )
  rm( "contractor.frm" )
  rm( "color.frm" )
  rm( "kin_relation.frm" )
  rm( "model_parameters.frm" )

  rm( file_name*".frm" )
  rm( file_name*".out" )
  rm( file_name*".log" )

  # write out
  jldopen( "integral_$(name_str).jld", "w" ) do file 
    write( file, "Generator", "FeAmGen.jl" )
    write( file, "n_loop", n_loop )
    write( file, "min_eps_xpt", min_eps_xpt )
    write( file, "max_eps_xpt", max_eps_xpt )
    write( file, "couplingfactor", "1" )
    write( file, "ext_mom_list", ext_mom_list )
    write( file, "scale2_list", map( string, scale2_list ) )
    write( file, "loop_den_list",  map( string, positive_loop_den_list ) )
    write( file, "loop_den_xpt_list", positive_loop_den_xpt_list )
    write( file, "kin_relation", map( p_->(string(p_[1]),string(p_[2])), collect(kin_relation) ) )
    write( file, "baseINC_script_str", string() )
    write( file, "model_parameter_dict", map( v_->(string(v_),"0"), ver_mass_list ) )
    write( file, "amp_color_list",  String["1"] )
    write( file, "amp_lorentz_list",  String[result_str] )
  end # file

  return nothing

end # function write_out_amplitude

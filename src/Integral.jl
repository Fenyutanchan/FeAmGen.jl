####################################
function has_ieta(
    den_index::Int64,
    den_list::Vector{Basic},
    ext_mom_list::Vector{Basic},
    kin_relation::Dict{Basic,Basic},
    ieta_scheme::Int64
)::Bool
####################################

  vanish_ext_mom_map = Dict{Basic,Basic}( map( m_ -> m_ => zero(Basic), ext_mom_list ) )
  den_ex = den_list[den_index]
  mom = get_args(den_ex)[1]
  loop_mom = subs( mom, vanish_ext_mom_map... )
  ext_mom = expand( mom - loop_mom )
  mass = get_args(den_ex)[2]

  if ieta_scheme == 0
    return false 
  elseif ieta_scheme == 1
    return true
  elseif ieta_scheme == 2
    if iszero(mass)  
      return false
    else
      return true
    end # if 
  elseif ieta_scheme == 3
    mass_list = map( den -> get_args(den)[2], den_list )
    kin_SP = subs( make_SP(ext_mom,ext_mom), kin_relation... )
    has_mass = !isempty( free_symbols( kin_SP) ∩ mass_list )
    if !has_mass && iszero(mass)
      return false
    else
      return true
    end # if
  elseif ieta_scheme > 10
    ieta_scheme_corr = ieta_scheme - 10
    bin_str = Base.bin( Unsigned(ieta_scheme_corr), 2, false )
    bin_str = "0"^(length(den_list)-length(bin_str)) * bin_str
    if bin_str[den_index] == '0'
      return false
    else
      return true
    end # if
  else
    error("Exception on ieta_scheme")
  end # if

end # function has_ieta




###################################
function generate_integral( 
    yaml_file::String )::Nothing
###################################

  dir_path_list = splitpath( yaml_file )[1:end-1]
  dir_path = isempty( dir_path_list ) ? string() : joinpath( dir_path_list... ) 

  file_dict = YAML.load_file( yaml_file; dicttype=OrderedDict{String,Any} ) 
  @assert collect(keys(file_dict)) == [
      "name", "n_loop", "min_ep_xpt", "max_ep_xpt", 
      "external_momenta", "kin_relation", "den_list", 
      "den_xpt_list", "numerator", "ieta_scheme", "comment"]

  name_str = file_dict["name"]
  n_loop = file_dict["n_loop"]::Int64
  min_ep_xpt = file_dict["min_ep_xpt"]::Int64
  max_ep_xpt = file_dict["max_ep_xpt"]::Int64
  ext_mom_list = map( Basic, file_dict["external_momenta"] )
  kin_relation = Dict{Basic,Basic}( map( p_->(Basic(p_[1]),Basic(p_[2])), file_dict["kin_relation"] ) )
  loop_den_list = map( Basic, file_dict["den_list"] )
  loop_den_xpt_list = file_dict["den_xpt_list"]
  ieta_scheme = file_dict["ieta_scheme"]

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


  @vars ieta, sqrteta
  @funs Den
  n_den = length(loop_den_list)
  positive_loop_den_list = Vector{Basic}()
  positive_loop_den_xpt_list = Vector{Int64}()
  irreducible_numerator = one(Basic)
  for den_index in 1:n_den
    one_den = loop_den_list[den_index]
    den_xpt = loop_den_xpt_list[den_index]

    den_arg_list = get_args(one_den)
    neg_mass2 = -den_arg_list[2]^2
    if has_ieta( den_index, loop_den_list, ext_mom_list, kin_relation, ieta_scheme )
      den_arg_list[3] = ieta
      neg_mass2 += im*sqrteta^2
    end # if
    new_den = Den( den_arg_list... )
  
    if den_xpt > 0 
      push!( positive_loop_den_list, new_den )
      push!( positive_loop_den_xpt_list, den_xpt )
    end # if  

    if den_xpt < 0
      mom2 = make_SP(den_arg_list[1],den_arg_list[1])
      irreducible_numerator *= ( mom2 + neg_mass2 )^(-den_xpt) 
    end # if
  end # for den_index



##positive_pos_list = findall( x_ -> x_ > 0, loop_den_xpt_list )
##positive_loop_den_list = loop_den_list[positive_pos_list]
##positive_loop_den_xpt_list = loop_den_xpt_list[positive_pos_list]

##negative_pos_list = findall( x_ -> x_ < 0, loop_den_xpt_list )
##irreducible_numerator = one(Basic)
##for negative_pos in negative_pos_list
##  negative_loop_den = loop_den_list[negative_pos]
##  negative_xpt = loop_den_xpt_list[negative_pos]

##  den_arg_list = get_args(negative_loop_den)
##  neg_mass2 = -den_arg_list[2]^2
##  # Now we still keep the i*eta in the irreducible numerator to be consistent with IBP reduction.
##  if ieta_scheme == 1 ||
##     ( ieta_scheme == 2 && !iszero(neg_mass2) )
##    @vars sqrteta
##    neg_mass2 += im*sqrteta^2
##  end # if

##  mom2 = make_SP(den_arg_list[1],den_arg_list[1])
##  irreducible_numerator *= ( mom2 + neg_mass2 )^(-negative_xpt) 
##end # for negative_pos
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

  run( pipeline( `form $(file_name).frm`, file_name*".log" ) )
  @info "[ Done FROM script execution ]" script="$(file_name).frm"

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
  jldopen( joinpath( dir_path, "integral$(name_str).jld" ), "w" ) do file 
    write( file, "Generator", "FeAmGen.jl" )
    write( file, "ieta_scheme", ieta_scheme )
    write( file, "n_loop", n_loop )
    write( file, "min_ep_xpt", min_ep_xpt )
    write( file, "max_ep_xpt", max_ep_xpt )
    write( file, "couplingfactor", "1" )
    write( file, "ext_mom_list", map( string, ext_mom_list ) )
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

end # function generate_integral








#-------------------------------------------------------------
function generate_multi_yaml( 
    original_yaml::String, 
    indices_list::Vector{Vector{Int64}}, 
    target_dir::String 
)::Vector{String}
#-------------------------------------------------------------

  if isdir( target_dir ) || isfile( target_dir )
    mv( target_dir, "$(target_dir)_$(now())" )
  end # if
  mkdir( target_dir )

  file_dict = YAML.load_file( original_yaml; dicttype=OrderedDict{String,Any} )
  @assert collect(keys(file_dict)) == [
      "name", "n_loop", "min_ep_xpt", "max_ep_xpt", 
      "external_momenta", "kin_relation", "den_list", 
      "den_xpt_list", "numerator", "ieta_scheme", "comment"]

  original_name = split( file_dict["name"], "_" )[1]
  new_multi_yaml_list = Vector{String}()
  for one_indices in indices_list
    one_name = "$(original_name)_$( join( map( string, one_indices ), "," ) )"
    file_dict["name"] = one_name
    file_dict["den_xpt_list"] = one_indices
    file_dict["comment"] = "Generated by FeAmGen function generate_multi_yaml @ $(now())"
    new_yaml = joinpath( target_dir, "$(one_name).yaml" )
    YAML.write_file( new_yaml, file_dict )
    push!( new_multi_yaml_list, new_yaml )
  end # for indices

  return new_multi_yaml_list

end # function generate_multi_yaml





#---------------------------------------------------------------------
"""
    generate_shiftUP_yaml( scalar_yaml_list::Vector{String}, target_dir::String )::Vector{String}

This is used to generate the derivative to the set of master integrals.
[0,-1,1,1] => [ [0,-1,2,1], [0,-1,1,2] ]
"""
function generate_shiftUP_yaml( scalar_yaml_list::Vector{String}, target_dir::String )::Vector{String}
#---------------------------------------------------------------------

  n_scalar = length(scalar_yaml_list)
  @assert n_scalar > 0

  first_scalar_yaml = scalar_yaml_list[1]
  #------------------------------------------------------------------
  @assert isfile(first_scalar_yaml) "The first argument is not a file!"
  first_dict = YAML.load_file( first_scalar_yaml )
  #------------------------------------------------------------------

  delete!( first_dict, "name" )
  delete!( first_dict, "den_xpt_list" )
  delete!( first_dict, "comment" )

  shiftUP_indices_set = Set{Vector{Int64}}()
  origin_indices_set = Set{Vector{Int64}}()
  for one_yaml in scalar_yaml_list
    one_dict = YAML.load_file( one_yaml )
    den_xpt_list = one_dict["den_xpt_list"]

    delete!( one_dict, "name" )
    delete!( one_dict, "den_xpt_list" )
    delete!( one_dict, "comment" )
    @assert one_dict == first_dict

    push!( origin_indices_set, den_xpt_list )

    n_xpt = length(den_xpt_list)
    for shift_pos in 1:n_xpt
      if den_xpt_list[shift_pos] <= 0 
        continue
      end # if
      new_xpt_list = copy( den_xpt_list )
      new_xpt_list[shift_pos] += 1
      push!( shiftUP_indices_set, new_xpt_list )
    end # for shift_pos
  end # for one_yaml

  new_indices_set = setdiff( shiftUP_indices_set, origin_indices_set )

  new_multi_yaml_list = generate_multi_yaml( first_scalar_yaml, collect(new_indices_set), target_dir )

  return new_multi_yaml_list

end # function generate_shiftUP_yaml

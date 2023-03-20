

####################################
# ieta_scheme 
# 0: none has iη 
# 1: all have iη
# 2: massive has iη
# >10: e.g. Int64(0b010100100)+10, indexing position of iη via binary number
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
  loop_mom = subs( mom, vanish_ext_mom_map )
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
    kin_SP = subs( make_SP(ext_mom,ext_mom), kin_relation )
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





######################################
function gen_vac_reduction_ieta(
    topology::Vector{Basic}
)::Tuple{ Vector{Vector{Basic}}, Vector{Basic} }
######################################

  n_loop = get_n_loop(topology)
  qi_list = [ Basic("q$index") for index in 1:n_loop ]
  mom_list = map( first∘get_args, topology )
  mass_list = (free_symbols∘map)( x->get_args(x)[2], topology )
  var_list = free_symbols( mom_list )
  ext_mom_list = setdiff( var_list, qi_list )

  null_dict = Dict( union( ext_mom_list .=> zeros(Basic), mass_list .=> zero(Basic) ) )
  vac_den_list = (unique∘map)( x->subs(x,null_dict), topology )
  vac_top_list = get_vac_top_list(vac_den_list)
  vac_top_list = filter( top -> any( !iszero, (last∘get_args).(top) ), vac_top_list )

  bk_mkdir( "vac_reduction" )
  cd( "vac_reduction" )
  mkdir( "config" )

  sector_str = join( map(string,ones( Int64, n_loop+div(n_loop*(n_loop-1),2) )) )
  file = open( "config/integralfamilies.yaml", "w" )
  write( file, """ 
  integralfamilies:
  """ )
  for vac_index in 1:length(vac_top_list)
  vac_top = vac_top_list[vac_index]
    write( file, """
      - name: "vac$(n_loop)loopT$(vac_index)"
        loop_momenta: [$(join(map(string,qi_list),","))]
        top_level_sectors: [b$(sector_str)]
        propagators:   
    """ )
  for one_den in vac_top
    mom, mass, width = get_args(one_den)
    write( file, """
          - [ "$mom", "$( iszero(width) ? "0" : "nim" )" ]
    """ )
  end # for one_den
  end # for vac_index
  close( file )

  file = open( "config/kinematics.yaml", "w" )
  write( file, """
  kinematics :
    incoming_momenta: []
    outgoing_momenta: []
    momentum_conservation: []
    kinematic_invariants:
      - [nim, 2]
    scalarproduct_rules: []
  """ )
  close( file )

  file = open( "masters.yaml", "w" )
  write( file, """
  # rmax: the maximal sum of positive propagator powers in the seed.
  # smax: the maximal negative sum of negative propagator powers in the seed.
  
  jobs:
    - reduce_sectors:
        reduce:
          - {sectors: [b$(sector_str)], r: 8, s: 2}
        run_initiate: true
        run_triangular: false
        run_back_substitution: false
        integral_ordering: 2
  """ )
  close( file )

  #-------------------
  run( `kira --parallel=$(Threads.nthreads()) masters.yaml` )
  #-------------------

  dir_list = (collect∘walkdir)("results")
  masters_list = filter( x->last(x)==["masters"], dir_list )
  file_list = map( x->first(x)*"/masters", masters_list )
  union_master_list = Vector{Basic}()
  for file_name in file_list
    file = open( file_name, "r" )
    content_str = read( file, String )
    close( file )
    content_str = seq_replace( content_str, "["=>"(", "]"=>")" )
    split_list = filter( !isempty, split( content_str, "\n" ) )
    master_list = map( x->(Basic∘string∘first∘split)(x,"#"), split_list )
    union_master_list = union( union_master_list, master_list )
  end # for file_name

  println( "[ vaccum masters ]" )
  map( println, union_master_list )
  println()

  cd( ".." )

  println( "[ Vaccum Master Integrals ]")
  map( println, union_master_list )
  println()
  for vac_index in 1:length(vac_top_list)
    vac_top = vac_top_list[vac_index]
    println( "Vaccum Topology #$(vac_index)" )
    map( println, vac_top )
    println()
  end # for vac_top

  return vac_top_list, union_master_list

end # function gen_vac_reduction_ieta




##################################
function get_vac_top_list(
    vac_den_list::Vector{Basic}
)::Vector{ Vector{Basic} }
##################################

  n_den = length(vac_den_list)
  vac_mom_list = map( first∘get_args, vac_den_list )
  unique_vac_mom_list = unique(vac_mom_list )

  if length(vac_mom_list) == length(unique_vac_mom_list)
    return [ copy(vac_den_list) ]
  end # if

  pos_list = findall( x->count(==(x),vac_mom_list)==2, unique_vac_mom_list )
  if length(pos_list) == 1
    split_mom = unique_vac_mom_list[first(pos_list)]
    split_mom_pos_list = findall( ==(split_mom), vac_mom_list )
    @assert length(split_mom_pos_list) == 2
    split_pos1 = split_mom_pos_list[1]
    split_pos2 = split_mom_pos_list[2]
    return [ vac_den_list[setdiff(1:n_den,split_pos1)],
             vac_den_list[setdiff(1:n_den,split_pos2)] ]
  end # if

  # more than one
  split_mom_list = unique_vac_mom_list[pos_list]
  split_mom_pos_list_collect = map( x->findall(==(x),vac_mom_list), split_mom_list )
  drop_choice_list = (vec∘collect∘Base.product)( split_mom_pos_list_collect... )
  @assert length(drop_choice_list) == 2^length(split_mom_list)
  @assert all( x->length(x)==length(split_mom_list), drop_choice_list )

  vac_top_list = Vector{Vector{Basic}}()
  for drop_choice in drop_choice_list
    push!( vac_top_list, vac_den_list[setdiff(1:n_den,drop_choice)] )
  end # for drop_choice

  return vac_top_list

end # function get_vac_top_list















###################################
function generate_integral( 
    yaml_file::String, 
    vac_top_list::Vector{Vector{Basic}}, 
    vac_master_list::Vector{Basic} 
)::Nothing
###################################

  dir_path_list = splitpath( yaml_file )[1:end-1]
  dir_path = isempty( dir_path_list ) ? string() : joinpath( dir_path_list... ) 

  file_dict = YAML.load_file( yaml_file; dicttype=OrderedDict{String,Any} ) 
  @assert collect(keys(file_dict)) == [
      "name", "n_loop", "min_ep_xpt", "max_ep_xpt", 
      "external_momenta", "kin_relation", "den_list", 
      "den_xpt_list", "numerator", "comment"]

  name_str = file_dict["name"]
  n_loop = file_dict["n_loop"]::Int64
  min_ep_xpt = file_dict["min_ep_xpt"]::Int64
  max_ep_xpt = file_dict["max_ep_xpt"]::Int64
  ext_mom_list = to_Basic( file_dict["external_momenta"] )
  kin_relation = to_Basic_dict( file_dict["kin_relation"] ) 
  loop_den_list = to_Basic( file_dict["den_list"] )
  loop_den_xpt_list = file_dict["den_xpt_list"]

  #-----------------------
  ver_mass_list = (free_symbols∘vcat)( 
                      (collect∘values)(kin_relation), 
                      map(x->get_args(x)[2],loop_den_list) ) 


  qi_list = [ Basic("q$index") for index in 1:n_loop ]
  n_ext = (length∘setdiff)( (free_symbols∘map)( first∘get_args, loop_den_list ), qi_list )
  n_den = length(loop_den_list)
  @assert n_den == n_loop + div(n_loop*(n_loop-1),2) + n_loop*n_ext

  #------------------------------------------------------------------
  # Prepare variable list and kinematic combo list
  var_list = (union∘free_symbols∘collect∘filter)( !is_FunctionSymbol, keys(kin_relation) )
  ver_mass_list = setdiff( ver_mass_list, var_list )

  var_list = union( var_list, ver_mass_list )
  var_list = sort( collect(var_list), by=gen_sorted_str )
  var_str_list = map( string, var_list )

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

    mom, mass, width = get_args(one_den)
    neg_mass2 = -mass^2
    if width == ieta
      neg_mass2 += im*sqrteta^2
    end # if
    #new_den = one_den
  
    if den_xpt > 0 
      push!( positive_loop_den_list, one_den )
      push!( positive_loop_den_xpt_list, den_xpt )
    end # if  

    if den_xpt < 0
      mom2 = make_SP(mom,mom)
      mom2 = subs( mom2, kin_relation )
      irreducible_numerator *= ( mom2 + neg_mass2 )^(-den_xpt) 
    end # if
  end # for den_index

  numerator_expr = irreducible_numerator * Basic( file_dict["numerator"] )

  #--------------
  form_script_str = make_amp_contraction_script( numerator_expr, kin_relation, ver_mass_list )

  result_io = IOBuffer()

  art_dir = Pkg.Artifacts.artifact"FeAmGen"
  cp( "$(art_dir)/scripts/contractor.frm", "contractor.frm", force=true )
  cp( "$(art_dir)/scripts/color.frm", "color.frm", force=true )

  try
    run( pipeline( `$(form()) -q -`; stdin=IOBuffer(form_script_str), stdout=result_io ) )
  catch
    file_name = "numerator_contraction"
    write( "$(file_name).frm", form_script_str )
    rethrow()
  end

  result_str = replace( (String∘take!)(result_io), "Coeff" => "" )

  # remove intermediate files
  rm( "contractor.frm" )
  rm( "color.frm" )


  # write out
  jldopen( joinpath( dir_path, "integral$(name_str).jld2" ), "w" ) do file 
    write( file, "Generator", "FeAmGen.jl function generate_integral" )
    # n_inc is used to generate mom_conserv. 
    # And this is supposed to be scalar integral, 
    # where the mom_conserv has been implemented already.
    write( file, "n_inc", 0 ) 
    write( file, "n_loop", n_loop )
    write( file, "min_ep_xpt", min_ep_xpt )
    write( file, "max_ep_xpt", max_ep_xpt )
    write( file, "couplingfactor", "1" )
    write( file, "ext_mom_list", to_String(ext_mom_list) )
    write( file, "scale2_list", to_String(scale2_list) )
    write( file, "topology", to_String(loop_den_list) )
    write( file, "loop_den_list", to_String(positive_loop_den_list) )
    write( file, "loop_den_xpt_list", positive_loop_den_xpt_list )
    write( file, "kin_relation", to_String_dict(kin_relation) ) 
    write( file, "baseINC_script_str", string() )
    write( file, "model_parameter_dict", (Dict∘map)( x->string(x)=>"0", ver_mass_list ) )
    write( file, "amp_color_list",  String["1"] )
    write( file, "amp_lorentz_list",  String[result_str] )

    write( file, "vac_master_list", to_String(vac_master_list) )
    write( file, "n_vac_top", length(vac_top_list) )
    for vac_index in 1:length(vac_top_list)
      vac_top = vac_top_list[vac_index]
      write( file, "vac_top$(vac_index)", to_String(vac_top ) )
    end # for vac_top
  end # file

  return nothing

end # function generate_integral








#########################################
function generate_multi_yaml( 
    original_yaml::String, 
    indices_list::Vector{Vector{Int64}}, 
    target_dir::String 
)::Vector{String}
#########################################

  if isdir( target_dir ) || isfile( target_dir )
    mv( target_dir, "$(target_dir)_$(now())" )
  end # if
  mkdir( target_dir )

  file_dict = YAML.load_file( original_yaml; dicttype=OrderedDict{String,Any} )
  @assert collect(keys(file_dict)) == [
      "name", "n_loop", "min_ep_xpt", "max_ep_xpt", 
      "external_momenta", "kin_relation", "den_list", 
      "den_xpt_list", "numerator", "comment"]

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
    generate_shiftUP_yaml( 
        scalar_yaml_list::Vector{String}, 
        target_dir::String 
    )::Vector{String}

This is used to generate the derivative to the set of master integrals.
[0,-1,1,1] => [ [0,-1,2,1], [0,-1,1,2] ]
"""
function generate_shiftUP_yaml( 
    scalar_yaml_list::Vector{String}, 
    target_dir::String 
)::Vector{String}
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

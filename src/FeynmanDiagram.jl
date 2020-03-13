
################################################################
function prepare_qgraf_dat( model::Model, input::Dict{Any,Any} )
################################################################

  inc_part_list = map( s_ -> model.particle_name_dict[s_], input["incoming"] )
  inc_idx_str_list = map( i_ -> string(i_), eachindex(inc_part_list) )
  inc_mom_str_list = map( (i_,p_) -> is_massless(p_) ? p_.name*"[k"*i_*"]" : p_.name*"[K"*i_*"]", 
                          inc_idx_str_list, inc_part_list )

  n_inc = length(inc_part_list)
  out_part_list = map( s_ -> model.particle_name_dict[s_], input["outgoing"] )
  out_idx_str_list = map( i_ -> string(i_+n_inc), eachindex(out_part_list) )
  out_mom_str_list = map( (i_,p_) -> is_massless(p_) ? p_.name*"[k"*i_*"]" : p_.name*"[K"*i_*"]", 
                          out_idx_str_list, out_part_list )

  file = open( "qgraf.dat", "w" )
  write( file, 
    "output='qgraf_out.dat';\n"*
    "style='miracle.sty';\n"*
    "model='model.qgraf';\n"*
    "in = "*join( inc_mom_str_list, "," )*";\n"*
    "out = "*join( out_mom_str_list, "," )*";\n"*
    "loops = "*string(input["n_loop"]+input["QCDCT_order"])*";\n"*
    "loop_momentum = q;\n" )
  if input["DropTadpole"] == true && input["DropWFcorrection"] == true 
    write( file, 
    "options = notadpole, onshell;\n" )
  elseif input["DropTadpole"] == true && input["DropWFcorrection"] == false 
    write( file, 
    "options = notadpole;\n" )
  elseif input["DropTadpole"] == false && input["DropWFcorrection"] == true 
    write( file, 
    "options = onshell;\n" )
  else 
    write( file, 
    "options =;\n" )
  end # if
  write( file, 
    "true = vsum[ gspow, "*string(input["Amp_QCD_order"])*", "*string(input["Amp_QCD_order"])*"];\n"*
    "true = vsum[ epow, "*string(input["Amp_QED_order"])*", "*string(input["Amp_QED_order"])*"];\n"*
    "true = vsum[ qcdctpow, "*string(input["QCDCT_order"])*", "*string(input["QCDCT_order"])*"];\n" )
  close(file)

  file = open( "miracle.sty", "w" )
  write( file, """
<prologue>
#
# file generated by <program>
#
<command_loop><command_line_loop>#<command_data><end><end>#

FeynmanDiagrams:
<diagram>
- # Diagram<diagram_index>
  diagram_index: <diagram_index>
  symmetry_factor: <symmetry_factor>
  sign: <sign>1
  incoming_propagators:<in_loop>
    - # incoming particle <in_index>
      in_index: <in_index>
      vertex_index: <vertex_index>
      field: "<field>"
      field_index: <field_index>
      momentum: "<momentum>"<end>
  outgoing_propagators:<out_loop>
    - # outgoing particle <out_index>
      out_index: <out_index>
      vertex_index: <vertex_index>
      field: "<field>"
      field_index: <field_index>
      momentum: "<momentum>"<end>
  remnant_propagators:<propagator_loop> 
    - # internal or loop propagator <propagator_index>
      propagator_index: <propagator_index>
      birth_index: <dual-vertex_index>
      death_index: <vertex_index>
      field: "<field>"
      momentum: "<momentum>"<end>
  vertices:<vertex_loop>
    - # vertex <vertex_index>
      vertex_index: <vertex_index>
      momentum_product: [[<ray_loop>"(<momentum>)",<end><back>]]
      fields: [[<ray_loop>"<field>",<end><back>]]
      propagator_index_list: [[<ray_loop><propagator_index>,<end><back>]]<end>
<epilogue>


<exit>

  """)
  close(file)


end # function prepare_qgraf_dat


#######################################################################
function generate_Feynman_diagram( model::Model, input::Dict{Any,Any} )
#######################################################################

  println()
  printstyled( "[Generate Feynman diagrams]\n", color=:green )

  prepare_qgraf_dat( model, input )

  if isfile( "qgraf_out.dat" ) == true
    rm( "qgraf_out.dat" )
  end # if
  run(`qgraf`)
  @assert isfile( "qgraf_out.dat" )

  rm( "qgraf.dat" )
  rm( "model.qgraf" )
  rm( "miracle.sty" )

end # function generate_Feynman_diagram





###################################################################################################
"""
Get the interaction according to the link field name list.
"""
###################################################################################################
function get_interaction( field_name_list::Vector{String}, model::Model )::Tuple{Interaction,Int64}
###################################################################################################

  QCDct_link_name_list = filter( s_ -> s_ in ["QCDct1","QCDct2"], field_name_list )
  QCDct_order = length( QCDct_link_name_list )

  normal_link_name_list = filter( s_ -> (s_ in ["QCDct1","QCDct2"]) == false, field_name_list ) 
  normal_link_kf_list = map( s_ -> model.particle_name_dict[s_].kf, normal_link_name_list )
  sorted_normal_link_kf_list = sort( normal_link_kf_list )
  inter = model.sorted_kf_list_dict[sorted_normal_link_kf_list]

  return inter, QCDct_order
end # function get_interaction
#################################################




###########################################################################################################
function get_incoming_couplings_lorentz_list( part::Particle, mark::Int64, momentum::Basic )::Vector{Basic}
###########################################################################################################

  if part.spin == :fermion && part.kf > 0 
    return [ Basic(" SpU( $mark, spb$mark, $momentum, r$mark, $(part.mass) ) ") ]
  elseif part.spin == :fermion && part.kf < 0
    return [ Basic(" SpVB( $mark, spb$mark, $momentum, r$mark, $(part.mass) ) ") ]
  elseif part.spin == :vector
    return [ Basic(" VecEps( $mark, spb$mark, $momentum, r$mark, $(part.mass) ) ") ]
  elseif part.spin == :scalar
    return [ Basic("1") ]
  else
    @assert false "We should have not considered ghost in external field."
  end # if

end # function get_incoming_couplings_lorentz_list


###########################################################################################################
function get_outgoing_couplings_lorentz_list( part::Particle, mark::Int64, momentum::Basic )::Vector{Basic}
###########################################################################################################

  if part.spin == :fermion && part.kf > 0
    return [ Basic(" SpUB( $mark, spa$mark, $momentum, r$mark, $(part.mass) ) ") ]
  elseif part.spin == :fermion && part.kf < 0
    return [ Basic(" SpV( $mark, spa$mark, $momentum, r$mark, $(part.mass) ) ") ]
  elseif part.spin == :vector
    return [ Basic(" VecEpsC( $mark, spa$mark, $momentum, r$mark, $(part.mass) ) ") ]
  elseif part.spin == :scalar
    return [ Basic("1") ]
  else
    printstyled( "We should have not considered ghost in external field.\n", color=:red )
    exit()
  end # if

end # function get_outgoing_couplings_lorentz_list


##########################################################################################################
function get_remnant_couplings_lorentz_list( part::Particle, mark::Int64, momentum::Basic, use_unitary_gauge::Bool )::Vector{Basic}
##########################################################################################################

  if part.spin == :fermion
    return [ Basic(" I*( GAij(spb$mark,spa$mark,$momentum)+ONEij(spb$mark,spa$mark)*$(part.mass) )*Den($momentum,$(part.mass),$(part.width)) ") ]
  elseif part.spin == :vector && is_massless(part) == true 
    return [ Basic(" (-1)*I*LMT(mua$mark,mub$mark)*Den($momentum,0,0) ") ]
  elseif part.spin == :vector && is_massless(part) == false 
    if use_unitary_gauge == true 
      return [ Basic(" (-1)*I*Den($momentum,$(part.mass),$(part.width))*( LMT(mua$mark,mub$mark)-FV($momentum,mua$mark)*FV($momentum,mub$mark)*$(part.mass)^(-2) ) ") ]
    else
      return [ Basic(" (-1)*I*LMT(mua$mark,mub$mark)*Den($momentum,$(part.mass),$(part.width)) ") ]
    end # if
  elseif part.spin in [:scalar, :ghost]
    return [ Basic(" I*Den( $momentum, $(part.mass), $(part.width) ) ") ]
  else 
    printstyled( "We should have not considered ghost in external field.\n", color=:red )
    exit()
  end # if

end # function get_remnant_couplings_lorentz_list



###########################################################################################
function edge_from_link_index( vert::ExVertex, link_index::Int64, g::GenericGraph )::ExEdge
###########################################################################################

  edge_list = edges(g)
  propagator_index = vert.attributes["propagator_index_list"][link_index]
  edge_idx = findfirst( e_ -> e_.attributes["propagator_index"] == propagator_index, edge_list )
  edge = edge_list[edge_idx]

  return edge
end # function edge_from_link_index


###########################################################################################
function vertex_from_mark( mark::Int64, g::GenericGraph )::ExVertex
###########################################################################################

  vertex_list = vertices(g)
  vertex_idx = findfirst( v_ -> v_.attributes["mark"] == mark, vertex_list )
  vert = vertex_list[vertex_idx]

  return vert
end # function vertex_from_mark

###########################################################################################
function vertex_from_label( label::String, g::GenericGraph )::ExVertex
###########################################################################################

  vertex_list = vertices(g)
  vertex_idx = findfirst( v_ -> v_.label == label, vertex_list )
  vert = vertex_list[vertex_idx]

  return vert
end # function vertex_from_label



#################################################################################
function get_link_color( g::GenericGraph, vert::ExVertex, link_index::Int64 )::Basic
################################################################################

  in_edge_list = in_edges( vert, g )

  if link_index < 0
    link_color = Basic("clv$(vert.attributes["mark"]*100+abs(link_index))")
  else 
    link_edge = edge_from_link_index( vert, link_index, g )
    link_color = link_edge in in_edge_list ? link_edge.attributes["death_COLOR"] : link_edge.attributes["birth_COLOR"]
  end # if 

  return link_color
end # function get_link_color


#########################################################################################
function translate_color_factor( one_color::Basic, vert::ExVertex, g::GenericGraph, )::Basic
#########################################################################################

  @funs Identity DeltaFun DeltaAdj T SUNT f SUNF

  color_str = string(one_color)
  new_color = one_color

  range_list = findall( r"Identity\([+-]*\d+, [+-]*\d+\)", color_str )
  Identity_str_list = map( r_ -> color_str[r_], range_list ) 

  for one_Identity_str in Identity_str_list 
    args = get_args( Basic(one_Identity_str) )
    link1_index, link2_index = convert(Int64,args[1]), convert(Int64,args[2])

    # It is not possible that both of link1_index and link2_index are negative (dummy).
    @assert link1_index > 0 || link2_index > 0 
    if link1_index > 0 
      edge1or2 = edge_from_link_index( vert, link1_index, g )
    elseif link2_index > 0
      edge1or2 = edge_from_link_index( vert, link2_index, g )
    else 
      @assert false
    end # if

    # This first case is considering the tadpole sunset diagram.
    if link1_index == link2_index
      color1, color2 = edge1or2.attributes["birth_COLOR"], edge1or2.attributes["death_COLOR"]
    # Then in the other case, we need to consider the case with dummy index.
    else
      color1 = get_link_color( g, vert, link1_index )
      color2 = get_link_color( g, vert, link2_index )
    end # if
    
    edge1or2_color = edge1or2.attributes["particle"].color
    @assert edge1or2_color != :singlet
    new_color = subs( new_color, Basic(one_Identity_str), 
                      edge1or2_color == :triplet ? DeltaFun(color1,color2) : DeltaAdj(color1,color2) )
  end # for one_Identity_str

  range_list = findall( r"[Tf]+\([+-]*\d+, [+-]*\d+, [+-]*\d+\)", color_str )
  Tf_str_list = map( r_ -> color_str[r_], range_list ) 

  for one_Tf_str in Tf_str_list 
    args = get_args( Basic(one_Tf_str) )
    link1_index, link2_index, link3_index = convert(Int64,args[1]), convert(Int64,args[2]), convert(Int64,args[3])

    # It is not possible that all of link1_index and link2_index and link3_index are negative (dummy).
    @assert link1_index > 0 || link2_index > 0 || link3_index > 0

    # two of [link1_index,link2_index,link3_index] are same 
    if link1_index == link2_index
      edge12 = edge_from_link_index( vert, link1_index, g )
      color1, color2 = edge12.attributes["birth_COLOR"], edge12.attributes["death_COLOR"]
      color3 = get_link_color( g, vert, link3_index )
    elseif link1_index == link3_index
      edge13 = edge_from_link_index( vert, link1_index, g )
      color1, color3 = edge13.attributes["birth_COLOR"], edge13.attributes["death_COLOR"]
      color2 = get_link_color( g, vert, link2_index )
    elseif link2_index == link3_index
      edge23 = edge_from_link_index( vert, link2_index, g )
      color2, color3 = edge23.attributes["birth_COLOR"], edge23.attributes["death_COLOR"]
      color1 = get_link_color( g, vert, link1_index )
    # Then none of them are same.
    else
      color1 = get_link_color( g, vert, link1_index )
      color2 = get_link_color( g, vert, link2_index )
      color3 = get_link_color( g, vert, link3_index )
    end # if
    
    new_color = subs( new_color, Basic(one_Tf_str), Basic("SUN"*one_Tf_str[1]*"($color1,$color2,$color3)") )
  end # for one_T_str

  return new_color
end # function translate_color_factor



######################################################################################
function get_link_lorentz( g::GenericGraph, vert::ExVertex, link_index::Int64 )::Basic
######################################################################################

  in_edge_list = in_edges( vert, g )

  if link_index < 0
    link_lor = Basic("muv$(vert.attributes["mark"]*100+abs(link_index))")
  else 
    link_edge = edge_from_link_index( vert, link_index, g )
    link_lor = link_edge in in_edge_list ? link_edge.attributes["death_LORENTZ"] : link_edge.attributes["birth_LORENTZ"]
  end # if 

  return link_lor
end # function get_link_lorentz


#######################################################################################
function get_link_momentum( g::GenericGraph, vert::ExVertex, link_index::Int64 )::Basic
#######################################################################################

  in_edge_list = in_edges( vert, g )

  @assert link_index > 0
  link_edge = edge_from_link_index( vert, link_index, g )
  # momentum is flow out in UFO Feynman rules
  link_momentum = ( link_edge in in_edge_list ? (-1) : 1 ) * link_edge.attributes["momentum"]

  return link_momentum
end # function get_link_momentum



#####################################################################################
function get_link_spinor( g::GenericGraph, vert::ExVertex, link_index::Int64 )::Basic
#####################################################################################

  in_edge_list = in_edges( vert, g )

  if link_index < 0
    link_spinor = Basic("spv$(vert.attributes["mark"]*100+abs(link_index))")
  else 
    link_edge = edge_from_link_index( vert, link_index, g )
    link_spinor = link_edge in in_edge_list ? link_edge.attributes["death_SPINOR"] : link_edge.attributes["birth_SPINOR"]
  end # if 

  return link_spinor
end # function get_link_spinor





################################################################################################
function translate_lorentz_factor( one_lorentz::Basic, vert::ExVertex, g::GenericGraph, )::Basic
################################################################################################

  @funs Metric LMT P FV Gamma GAij ProjP PRij ProjM PLij Identity ONEij 

  lorentz_str = string(one_lorentz)
  new_lorentz = one_lorentz

  range_list = findall( r"Metric\([+-]*\d+, [+-]*\d+\)", lorentz_str )
  Metric_str_list = map( r_ -> lorentz_str[r_], range_list ) 
  for one_Metric_str in Metric_str_list
    args = get_args( Basic(one_Metric_str) )
    link1_index, link2_index = convert(Int64,args[1]), convert(Int64,args[2])
    @assert link1_index > 0 && link2_index > 0
    if link1_index > 0 
      edge1or2 = edge_from_link_index( vert, link1_index, g )
    elseif link2_index > 0
      edge1or2 = edge_from_link_index( vert, link2_index, g )
    else 
      @assert false
    end # if

    if link1_index == link2_index
      lor1, lor2 = edge1or2.attributes["birth_LORENTZ"], edge1or2.attributes["death_LORENTZ"]
    else
      lor1 = get_link_lorentz( g, vert, link1_index )
      lor2 = get_link_lorentz( g, vert, link2_index )
    end # if

    new_lorentz = subs( new_lorentz, Basic(one_Metric_str), LMT(lor1,lor2) )
  end # for one_Metric_str

  range_list = findall( r"P\([+-]*\d+, [+-]*\d+\)", lorentz_str )
  P_str_list = map( r_ -> lorentz_str[r_], range_list ) 
  for one_P_str in P_str_list
    args = get_args( Basic(one_P_str) )
    link1_index, link2_index = convert(Int64,args[1]), convert(Int64,args[2])

    mom1 = get_link_momentum( g, vert, link1_index )
    lor2 = get_link_lorentz( g, vert, link2_index )

    new_lorentz = subs( new_lorentz, Basic(one_P_str), FV(mom1,lor2) )
  end # for one_P_str

  range_list = findall( r"Gamma\([+-]*\d+, [+-]*\d+, [+-]*\d+\)", lorentz_str )
  Gamma_str_list = map( r_ -> lorentz_str[r_], range_list ) 
  for one_Gamma_str in Gamma_str_list
    args = get_args( Basic(one_Gamma_str) )
    link1_index, link2_index, link3_index = convert(Int64,args[1]), convert(Int64,args[2]), convert(Int64,args[3])

    lor1 = get_link_lorentz( g, vert, link1_index )
    sp2 = get_link_spinor( g, vert, link2_index )
    sp3 = get_link_spinor( g, vert, link3_index )

    new_lorentz = subs( new_lorentz, Basic(one_Gamma_str), GAij(sp2,sp3,lor1) )
  end # for one_Gamma_str

  range_list = findall( r"ProjP\([+-]*\d+, [+-]*\d+\)", lorentz_str )
  ProjP_str_list = map( r_ -> lorentz_str[r_], range_list ) 
  for one_ProjP_str in ProjP_str_list
    args = get_args( Basic(one_ProjP_str) )
    link1_index, link2_index = convert(Int64,args[1]), convert(Int64,args[2])

    sp1 = get_link_spinor( g, vert, link1_index )
    sp2 = get_link_spinor( g, vert, link2_index )

    new_lorentz = subs( new_lorentz, Basic(one_ProjP_str), PRij(sp1,sp2) )
  end # for one_ProjP_str

  range_list = findall( r"ProjM\([+-]*\d+, [+-]*\d+\)", lorentz_str )
  ProjM_str_list = map( r_ -> lorentz_str[r_], range_list ) 
  for one_ProjM_str in ProjM_str_list
    args = get_args( Basic(one_ProjM_str) )
    link1_index, link2_index = convert(Int64,args[1]), convert(Int64,args[2])

    sp1 = get_link_spinor( g, vert, link1_index )
    sp2 = get_link_spinor( g, vert, link2_index )

    new_lorentz = subs( new_lorentz, Basic(one_ProjM_str), PLij(sp1,sp2) )
  end # for one_ProjM_str

  range_list = findall( r"Identity\([+-]*\d+, [+-]*\d+\)", lorentz_str )
  Identity_str_list = map( r_ -> lorentz_str[r_], range_list ) 
  for one_Identity_str in Identity_str_list
    args = get_args( Basic(one_Identity_str) )
    link1_index, link2_index = convert(Int64,args[1]), convert(Int64,args[2])

    sp1 = get_link_spinor( g, vert, link1_index )
    sp2 = get_link_spinor( g, vert, link2_index )

    new_lorentz = subs( new_lorentz, Basic(one_Identity_str), ONEij(sp1,sp2) )
  end # for one_Identity_str

  return new_lorentz
end # function translate_lorentz_factor




######################################################################################################
function convert_qgraf_TO_Graph( one_qgraf::Dict{Any,Any}, model::Model )::Union{GenericGraph,Nothing}
######################################################################################################

  # Check if there is CT vertex and if the QCDct field is properly used for CT graphs.
  # QCDct1 and QCDct2 should have same end-points as designed.
  QCDct_propagators = filter( p_ -> p_["field"] in ["QCDct1","QCDct2"], one_qgraf["remnant_propagators"] )
  invalid_propagator_pos = findfirst( p_ -> p_["birth_index"] != p_["death_index"], QCDct_propagators )
  if invalid_propagator_pos != nothing 
    printstyled( "Found one invalid diagram!\n", color=:red )
    return nothing
  end # if

  g = graph( ExVertex[], ExEdge[] )  
  v0 = ExVertex(0,"graph property")
  add_vertex!(g,v0)
  ## Diagram index
  v0.attributes["diagram_index"] = one_qgraf["diagram_index"] 
  ## Read-in symmetry_factor.
  v0.attributes["symmetry_factor"] = one_qgraf["symmetry_factor"] 
  ## Read-in sign of this diagram. 
  ## The sign from QGRAF is combination of anti-commutivative field (fermion,ghost,QCDct1,QCDct2...) loops,
  ##   and the way open fermionic/ghost lines connect those external fields.
  ## Also correct the sign if there is QCDct1 or QCDct2, 
  ##   since they are defined as anti-commutative fields to avoid changing symmetry factor.
  n_QCDct = length( QCDct_propagators )
  v0.attributes["sign"] = one_qgraf["sign"]*(-1)^n_QCDct
  ## For sake of function vertex_from_mark, we also need "mark" for v0
  v0.attributes["mark"] = 0
  ## For filtering the "style"
  v0.attributes["style"] = ""

  qgraf_incoming_propagators = one_qgraf["incoming_propagators"]
  n_inc = length(qgraf_incoming_propagators)
  qgraf_outgoing_propagators = one_qgraf["outgoing_propagators"]
  n_out = length(qgraf_outgoing_propagators)

  ## For distinguishing incoming and outgoing edges
  v0.attributes["n_inc"] = n_inc
  v0.attributes["n_out"] = n_out



  ## Add incoming vertices.
  for one_inc in qgraf_incoming_propagators
    vert_idx = num_vertices(g) # There is v0 ahead already.
    inc_v = ExVertex( vert_idx, "v$vert_idx" )
    inc_v.attributes = Dict(
      "mark" => one_inc["in_index"],
      "style" => "External",
      "interaction" => nothing,
      "QCDct_order" => 0,
      "color_list" => [1],
      "couplings_lorentz_list" => [1]
    ) # end Dict
    add_vertex!( g, inc_v ) 
  end # for one_inc


  # Add outgoing vertices.
  for one_out in qgraf_outgoing_propagators
    vert_idx = num_vertices(g) # There is v0 ahead already.
    out_v = ExVertex( vert_idx, "v$vert_idx" )
    out_v.attributes = Dict(
      "mark" => one_out["out_index"]+n_inc,
      "style" => "External",
      "interaction" => nothing,
      "QCDct_order" => 0,
      "color_list" => [1],
      "couplings_lorentz_list" => [1] 
    ) # end Dict
    add_vertex!( g, out_v ) 
  end # for one_out


  #  Add internal or loop vertices.
  qgraf_vertices = one_qgraf["vertices"]
  for one_vert in qgraf_vertices
    vert_idx = num_vertices(g) # There is v0 ahead already.
    int_v = ExVertex( vert_idx, "v$vert_idx" )

    # Get interaction and QCDct_order
    inter, QCDct_order = get_interaction( one_vert["fields"], model )
    
    #-----------------------------------------------------------------------------------------------
    # Now we need to re-order the list "propagator_index_list" according to the vertex Feynman rule.
    #-----------------------------------------------------------------------------------------------
    # Combine "fields" and "propagator_index_list" in to Dict.
    # Example:
    #     "fields" => ["tbar", "b", "Wplus"] and "propagator_index_list" => [-4, -1, 1] 
    #       convert into [(tbar,-4), (b,-1), (Wplus,1)], where "fields" are also converted into Particle type.
    # Here the "fields" can contain "QCDct1" or "QCDct2".
    # Then use Dict{Particle,Int64}( collect(xx) ) to convert.
    link_index_pair_list = map( (s_, i_) -> (model.particle_name_dict[s_], i_), one_vert["fields"], one_vert["propagator_index_list"] )
    link_index_dict = Dict{Particle,Int64}( collect(link_index_pair_list) )
    # Then use Dict to find the relevant index list according to the particle order of Feynman rules vertex.
    # # Example:
    #     Dict(b=>-1, Wplus=>1, tbar=>-4) for the vertex [tbar,b,Wplus] should give index list [-4,-1,1].
    feynrules_index_list = map( p_ -> link_index_dict[p_], inter.link_list )

    # color_list and couplings_lorentz_list will be generated later.
    int_v.attributes = Dict(
      "mark" => one_vert["vertex_index"]+n_inc+n_out,
      "style" => "Internal",
      "interaction" => inter,
      "QCDct_order" => QCDct_order,
      "propagator_index_list" => feynrules_index_list,
      "color_list" => nothing,
      "couplings_lorentz_list" => nothing 
    ) # end Dict
    add_vertex!( g, int_v ) 
  end # for one_vert


  # Add incoming propagators.
  # It could be u(kf>0) or vbar(kf<0) for fermion, epsilon_\mu for vector, 1 for scalar.
  for one_inc in qgraf_incoming_propagators
    in_index = one_inc["in_index"]
    vert_index = one_inc["vertex_index"]+n_inc+n_out
    field_name = one_inc["field"]
    field_part = model.particle_name_dict[field_name]

    in_v = vertex_from_mark( in_index, g )
    vert_v = vertex_from_mark( vert_index, g )

    new_edge = ExEdge( num_edges(g)+1, in_v, vert_v )

    mark = in_index
    momentum = Basic(one_inc["momentum"])
    couplings_lorentz_list = get_incoming_couplings_lorentz_list( field_part, mark, momentum )

    new_edge.attributes = Dict( 
      "mark" => mark,
      "particle" => field_part,
      "style" => "External",
      "propagator_index" => one_inc["field_index"],
      "momentum" => momentum,
      "ref2_MOM" => Basic("r$mark"),
      "null_MOM" => Basic("k$mark"),
      "birth_LORENTZ" => Basic("mua$mark"), "birth_COLOR" => Basic("cla$mark"), "birth_SPINOR" => Basic("spa$mark"),
      "death_LORENTZ" => Basic("mub$mark"), "death_COLOR" => Basic("clb$mark"), "death_SPINOR" => Basic("spb$mark"),
      "color_list" => [ Basic("1") ],
      "couplings_lorentz_list" => couplings_lorentz_list 
    ) # end Dict

    add_edge!( g, new_edge )
  end # for one_inc


  # Add outgoing propagaotrs.
  # It could be ubar(kf>0) or v(kf<0) for fermion, epsilon_\mu^* for vector, 1 for scalar.
  for one_out in qgraf_outgoing_propagators
    out_index = one_out["out_index"]+n_inc
    vert_index = one_out["vertex_index"]+n_inc+n_out
    field_name = one_out["field"]
    field_part = model.particle_name_dict[field_name]

    vert_v = vertex_from_mark( vert_index, g )
    out_v = vertex_from_mark( out_index, g )

    new_edge = ExEdge( num_edges(g)+1, vert_v, out_v )

    mark = out_index
    momentum = Basic(one_out["momentum"])
    couplings_lorentz_list = get_outgoing_couplings_lorentz_list( field_part, mark, momentum )

    new_edge.attributes = Dict(
      "mark" => mark,
      "particle" => field_part,
      "style" => "External",
      "propagator_index" => one_out["field_index"],
      "momentum" => momentum,
      "ref2_MOM" => Basic("r$mark"),
      "null_MOM" => Basic("k$mark"),
      "birth_LORENTZ" => Basic("mua$mark"), "birth_COLOR" => Basic("cla$mark"), "birth_SPINOR" => Basic("spa$mark"),
      "death_LORENTZ" => Basic("mub$mark"), "death_COLOR" => Basic("clb$mark"), "death_SPINOR" => Basic("spb$mark"),
      "color_list" => [ Basic("1") ],
      "couplings_lorentz_list" => couplings_lorentz_list 
    ) # end Dict

    add_edge!( g, new_edge )
  end # for one_out


  # Filter out the "QCDct1" and "QCDct2" propagators.
  qgraf_remnant_propagators = one_qgraf["remnant_propagators"]
  true_remnant_propagators = filter( rem_ -> (rem_["field"] in ["QCDct1","QCDct2"]) == false, qgraf_remnant_propagators )
  # Add internal and loop propagators.
  for one_rem in true_remnant_propagators
    mark = one_rem["propagator_index"]+n_inc+n_out
    birth_index = one_rem["birth_index"]+n_inc+n_out
    death_index = one_rem["death_index"]+n_inc+n_out
    field_name = one_rem["field"]
    field_part = model.particle_name_dict[field_name]

    birth_v = vertex_from_mark( birth_index, g )
    death_v = vertex_from_mark( death_index, g )

    new_edge = field_part.kf > 0 ? ExEdge( num_edges(g)+1, birth_v, death_v ) : ExEdge( num_edges(g)+1, death_v, birth_v )

    momentum = sign(field_part.kf)*Basic(one_rem["momentum"])

    id_color_dict = Dict( :triplet => [ Basic(" DeltaFun(clb$mark,cla$mark) ") ], 
                          :octet   => [ Basic(" DeltaAdj(clb$mark,cla$mark) ") ],
                          :singlet => [ Basic("1") ] )

    couplings_lorentz_list = get_remnant_couplings_lorentz_list( field_part, mark, momentum, model.unitary_gauge )

    new_edge.attributes = Dict( 
      "mark" => mark,
      "particle" => model.particle_kf_dict[abs(field_part.kf)],
      "style" => findfirst("q",one_rem["momentum"]) == nothing ? "Internal" : "Loop",
      "propagator_index" => one_rem["propagator_index"], 
      "momentum" => momentum,
      "ref2_MOM" => Basic("0"),
      "null_MOM" => Basic("0"),
      "birth_LORENTZ" => Basic("mua$mark"), "birth_COLOR" => Basic("cla$mark"), "birth_SPINOR" => Basic("spa$mark"),
      "death_LORENTZ" => Basic("mub$mark"), "death_COLOR" => Basic("clb$mark"), "death_SPINOR" => Basic("spb$mark"),
      "color_list" => id_color_dict[field_part.color],
      "couplings_lorentz_list" => couplings_lorentz_list 
    ) # end Dict
    
    add_edge!( g, new_edge )
  end # for one_rem



# #-------------------------------------------------------------------
# # For convenience, we create indexing Dict for the edges according to :propagator_index
# propagator_index_pair_list = map( e_ -> get_prop(mg,e_,:propagator_index) => e_ , edges( mg ) )
# propagator_index_dict = Dict{Int64,Edge}( collect(propagator_index_pair_list) )
# set_prop!( mg, :propagator_index_dict, propagator_index_dict )
  #-------------------------------------------------------------------
  # Now the structure of this graph has been digested into Graph.
  # Then we can evaluate color_row_list and couplings_lorentz_col_list for the internal or loop vertices.
  #-------------------------------------------------------------------
  internal_vertex_list = filter( v_ -> v_.attributes["style"] == "Internal", vertices(g) )
  for vert in internal_vertex_list
    inter = vert.attributes["interaction"]

    new_color_row_list = map( color_ -> translate_color_factor(color_,vert,g), inter.color_row_list )

    new_lorentz_col_list = map( lor_ -> translate_lorentz_factor(lor_,vert,g), inter.lorentz_col_list )

    n_row, n_col = size( inter.couplings_matrix )
    new_couplings_lorentz_list = Array{Basic,1}(undef,n_row)
    @vars CTorder
    for r_ in 1:n_row
      new_couplings_lorentz_list[r_] = 0
      for c_ in 1:n_col
        if vert.attributes["QCDct_order"] == 0 
          new_couplings_lorentz_list[r_] += subs(inter.couplings_matrix[r_,c_],CTorder,0)*new_lorentz_col_list[c_]
        else 
          new_couplings_lorentz_list[r_] += coeff( expand(inter.couplings_matrix[r_,c_]), CTorder^(vert.attributes["QCDct_order"]) )*new_lorentz_col_list[c_]
        end # if
      end # for c_
    end # for r_

    vert.attributes["color_list"] = new_color_row_list
    vert.attributes["couplings_lorentz_list"] = new_couplings_lorentz_list
  end # for vert 

  return g

end # function convert_qgraf_TO_Graph



###########################################
"""
tensor_product( x, y )
x and y are two string lists/arrays.
This function calculate the tensor production of two arrays.
"""
function tensor_product( 
    ex_list1::Union{Vector{Basic},Vector{Int64}}, 
    ex_list2::Union{Vector{Basic},Vector{Int64}}
)::Vector{Basic}
###########################################
  res = Vector{Basic}()
  for ex1 in ex_list1, ex2 in ex_list2
    push!(res,ex1*ex2)
  end
  return res

end # function tensor_product




###################################################################
"""
Now this graph can be evaluated according to the values of the propagators and vertices.
"""
function assemble_amplitude( g::GenericGraph, model::Model )::Tuple{Vector{Basic},Vector{Basic}}
###################################################################

  amp_color_list = Basic[1]
  amp_couplings_lorentz_list = Basic[1]
  for vert in vertices(g)
    if vert.label == "graph property"
      continue
    end # if
    amp_color_list = tensor_product( amp_color_list, vert.attributes["color_list"] )
    amp_couplings_lorentz_list = tensor_product( amp_couplings_lorentz_list, vert.attributes["couplings_lorentz_list"] )
  end # for vert

  for edge in edges(g)
    amp_color_list = tensor_product( amp_color_list, edge.attributes["color_list"] )
    amp_couplings_lorentz_list = tensor_product( amp_couplings_lorentz_list, edge.attributes["couplings_lorentz_list"] )
  end # for edge

  return amp_color_list, amp_couplings_lorentz_list;

end # function assemble_amplitude





##################################################################################################
function canonicalize_loop_PowDen( lorentz_expr_list::Vector{Basic}, model::Model )::Vector{Basic}
##################################################################################################

  file = open( "model_parameters.frm", "w" )
  write( file, "symbol "*join( map( k_->string(k_), collect(keys(model.parameter_dict)) ), "," )*";\n" )
  close(file)

  file = open( "contractor.frm", "w" )
  write( file, make_contractor_script() )
  close(file)

  new_lorentz_expr_list = Vector{Basic}( undef, length(lorentz_expr_list) )
  for index in 1:length(lorentz_expr_list)
    lorentz_expr = lorentz_expr_list[index]
    file_name = "lorentz_expr$(index)"
    form_script_str = make_canonicalize_amplitude_script( lorentz_expr, file_name )

    file = open( file_name*".frm", "w" )
    write( file, form_script_str )
    close(file)

    printstyled( "[ form $(file_name).frm ]\n", color=:yellow )
    run( pipeline( `form $(file_name).frm`, file_name*".log" ) )

    file = open( file_name*".out", "r" )
    result_str = read( file, String )
    result_expr = Basic(result_str)
    new_lorentz_expr_list[index] = result_expr

    run( `rm $(file_name).frm $(file_name).out $(file_name).log` )
  end # for lorentz_expr

  return new_lorentz_expr_list

end # function canonicalize_loop_PowDen


##################################################################################
function contract_Dirac_indices( lorentz_expr_list::Vector{Basic}, graph::GenericGraph )::Vector{Basic}
##################################################################################

  # model_parameters.frm and contractor.frm have been created in canonicalize_loop_PowDen

  file = open( "baseINC.frm", "w" )
  write( file, make_baseINC_script( graph ) )
  close(file)

  new_lorentz_expr_list = Vector{Basic}( undef, length(lorentz_expr_list) )
  for index in 1:length(lorentz_expr_list)
    lorentz_expr = lorentz_expr_list[index]
    file_name = "contract_lorentz_expr$(index)"
    form_script_str = make_amp_contraction_script( lorentz_expr, file_name )

    file = open( file_name*".frm", "w" )
    write( file, form_script_str )
    close(file)

    printstyled( "[ form $(file_name).frm ]\n", color=:yellow )
    run( pipeline( `form $(file_name).frm`, file_name*".log" ) )

    file = open( file_name*".out", "r" )
    result_str = read( file, String )
    result_expr = Basic(result_str)
    new_lorentz_expr_list[index] = result_expr

    run( `rm $(file_name).frm $(file_name).out $(file_name).log` )
  end # for index
  rm( "contractor.frm" )

  return new_lorentz_expr_list

end # function contract_Dirac_indices




##########################################################################
function generate_amplitude( model::Model, input::Dict{Any,Any} )::Nothing
##########################################################################

  couplingfactor = Basic(input["couplingfactor"]) 

  qgraf_out = YAML.load( open("qgraf_out.dat") )

  qgraf_list = qgraf_out["FeynmanDiagrams"]

  #------------------------------------------------  
  # Convert qgraf to GenericGraph
  graph_set = Set{GenericGraph}()
  for one_qgraf in qgraf_list 

    g = convert_qgraf_TO_Graph( one_qgraf, model )
    if g == nothing 
      continue
    end # if

    push!( graph_set, g )
  end # for one_qgraf
  graph_list = sort( collect( graph_set ), by= g_->vertex_from_label("graph property",g_).attributes["diagram_index"] )


  #------------------------------------------------  
  # Generate Gauge choice
  gauge_choice = generate_gauge_choice( graph_list )
  # Generate kinematics relation
  kin_relation = generate_kin_relation( graph_list, gauge_choice )
  file = open( "kin_relation.frm", "w" )
  write( file, join( map( ele_->"id $(ele_[1]) = $(ele_[2]);", collect(kin_relation) ), "\n" ) )
  close(file)

  #------------------------------------------------  
  # Calculate amplitude for each graph
  amp_list = Vector{Amplitude}( undef, length(graph_list) )
  for index in 1:length(graph_list)
    g = graph_list[index]

    amp_color_list, amp_lorentz_list = assemble_amplitude( g, model )
    amp_lorentz_list = canonicalize_loop_PowDen( amp_lorentz_list, model )
    amp_lorentz_list = contract_Dirac_indices( amp_lorentz_list, g )

    amp_list[index] = Amplitude( g, amp_color_list, amp_lorentz_list )
  end # for g


  #------------------------------------------------  
  # Write out amplitude
  printstyled( "[ Generate amplitude.out ]\u264e\n", color=:green, bold=true )
  amp_file = open( "amplitude.out", "w" )
  write( amp_file, 
    "couplingfactor: $(couplingfactor)\n\n" )
  for amp in amp_list
    g = amp.graph

    diagram_index = vertex_from_label("graph property",g).attributes["diagram_index"]

    write( amp_file, 
      "Diagram #$(diagram_index): \n" )

    for ii in 1:length(amp.color_list)
      one_color = amp.color_list[ii]
      write( amp_file, 
      "  amp_color #$(ii): \n"*
      "    $(one_color); \n" )
    end # for ii

    for ii in 1:length(amp.lorentz_list)
      one_val = amp.lorentz_list[ii]
      write( amp_file, 
      "  amp_couplings_lorentz #$(ii): \n"* 
      "    $( one_val/couplingfactor ); \n" )
    end # for ii
    write( amp_file, "\n\n" )

  end # for g
  close( amp_file )

  #------------------------------------------------  
  # Write out visual graph
  printstyled( "[ Generate visual_graphs.tex ]\u264e\n", color=:green, bold=true )
  visual_file = open( "visual_graphs.tex", "w" )
  write( visual_file,
    "\\documentclass{article}\n"*
    "\\usepackage{tikz-feynman}\n"*
    "\n"*
    "\\begin{document}\n"*
    "\n" )
  for g in graph_list
    result_str = generate_visual_graph( g, model )
    write( visual_file, result_str ) 
  end # for one_qgraf

  write( visual_file, 
    "\\end{document} \n"*
    "\n" )
  close( visual_file )


  print( "\nRun \"lualatex visual_graphs.tex\" to generate PDF file.\n" )
  # run( `lualatex visual_graphs.tex` )

  return nothing

end # function generate_amplitude












#####################################################
function get_line_style_str( part::Particle )::String
#####################################################

 if part.spin == :scalar 
   if is_not_majorana(part) && part.kf > 0 
     return "charged scalar"
   elseif is_not_majorana(part) && part.kf < 0 
     return "anti charged scalar"
   else
     return "scalar"
   end # if
 elseif part.spin == :fermion 
   if is_not_majorana(part) && part.kf > 0 
     return "fermion"
   elseif is_not_majorana(part) && part.kf < 0 
     return "anti fermion"
   else
     return "majorana"
   end # if
 elseif part.spin == :vector
   if is_gluon(part) 
     return "gluon"
   else
     return "boson"
   end # if
 elseif part.spin == :ghost 
   return "ghost"
 else
   printstyled( "Spin exception!\n" )
   exit()
 end # if

end # function get_line_style_str



#######################################################################
function generate_visual_graph( g::GenericGraph, model::Model )::String
#######################################################################

  v0 = vertex_from_label( "graph property", g )

  result_str = 
    "\\begin{figure}[!htb] \n"*
    "\\begin{center} \n"*
    "\\feynmandiagram [large] { \n"*
    "  %$(string(v0)) \n"

  for att in v0.attributes
    result_str *= 
    "    %$(att[1]) => $(att[2]) \n"
  end # for att

  for edge in edges(g)
    result_str *= 
    "  %$(string(edge)) \n"
    for att in edge.attributes
      result_str *= 
      "    %$(att[1]) => $(att[2]) \n"
    end # for att
  end # for edge

  for vert in vertices(g)
    result_str *= 
    "  %$(string(vert)) \n"
    for att in vert.attributes
      result_str *= 
      "    %$(att[1]) => $(att[2]) \n"
    end # for att
  end # for vert

  QCDct_str_dict = Dict( 0 => "", 1 => " [crossed dot]", 2 => " [red, label = DOUBLE, crossed dot]" )

  edge_list = edges(g)
  for edge in edge_list
    src_v = source(edge)
    tgt_v = target(edge)
    src_mark = src_v.attributes["mark"]
    tgt_mark = tgt_v.attributes["mark"]
    mark_pair_set = Set([src_mark,tgt_mark])

    parallel_edge_list = filter( e_ -> Set([source(e_).attributes["mark"],target(e_).attributes["mark"]]) == mark_pair_set, edge_list )

    if length(parallel_edge_list) <= 1
      half_circle_option = ""
    else 
      parallel_edge_mark_list = map( e_ -> e_.attributes["mark"], parallel_edge_list )   
      max_mark = max(parallel_edge_mark_list...)
      min_mark = min(parallel_edge_mark_list...)
      if edge.attributes["mark"] == max_mark
        half_circle_option = src_mark > tgt_mark ? ", half left" : ", half right"
      elseif edge.attributes["mark"] == min_mark
        half_circle_option = src_mark > tgt_mark ? ", half right" : ", half left"
      else
        @assert false
      end # if
    end # end if

    src_QCDct_str = QCDct_str_dict[src_v.attributes["QCDct_order"]]
    src_mark = src_v.attributes["mark"]

    tgt_QCDct_str = QCDct_str_dict[tgt_v.attributes["QCDct_order"]]
    tgt_mark = tgt_v.attributes["mark"]

    edge_style_str = get_line_style_str(edge.attributes["particle"])

    tadpole_option = ""
    if src_mark == tgt_mark 
      tadpole_option = ", loop, min distance=3cm"
    end # if
   
    src_external_marking_str = ""
    if src_v.attributes["style"] == "External"
      src_external_marking_str = " [particle = \\($(src_mark)\\)]"
    end # if
    tgt_external_marking_str = ""
    if tgt_v.attributes["style"] == "External"
      tgt_external_marking_str = " [particle = \\($(tgt_mark)\\)]"
    end # if


    particle_name = edge.attributes["particle"].name
    particle_name = replace( particle_name, "plus" => "^{+}" )
    particle_name = replace( particle_name, "minus" => "^{-}" )
    particle_name = replace( particle_name, "ve" => "\\nu_{e}" )
    particle_name = replace( particle_name, "mu" => "\\mu" )
    particle_name = replace( particle_name, "ta" => "\\tau" )
    particle_name = replace( particle_name, "vm" => "\\nu_{\\mu}" )
    particle_name = replace( particle_name, "vt" => "\\nu_{\\tau}" )
    particle_name = replace( particle_name, "^a" => "\\gamma" )
    if length(particle_name) > 3 && particle_name[end-2:end] == "bar"
      particle_name = "\\overline{"*particle_name[1:end-3]*"}"
    end # if

    mom_str = replace( string(edge.attributes["momentum"]), r"([Kkq]+)(\d+)" => s"\1_{\2}" )

    result_str *= 
      "$(src_v.label)$(src_QCDct_str)$(src_external_marking_str) -- [$(edge_style_str)$(half_circle_option)$(tadpole_option), edge label' = \\($(particle_name)\\), momentum = \\($(mom_str)\\) ] $(tgt_v.label)$(tgt_QCDct_str)$(tgt_external_marking_str), \n"
  end # for edge

  result_str *= 
    "};\n"*
    "\\end{center}\n"*
    "\\caption{Diagram$(v0.attributes["diagram_index"]), Sign: $(v0.attributes["sign"]), Symmetry factor: $(v0.attributes["symmetry_factor"])}\n"*
    "\\end{figure}\n"*
    "\\newpage\n"

  return result_str

end # function generate_visual_graph


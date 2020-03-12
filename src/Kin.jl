
####################################################################
function generate_gauge_choice( graph_list::Vector{GenericGraph} )::Dict{Basic,Basic}
####################################################################

  # Only the external fields are needed
  graph0 = first( graph_list )

  ext_edge_list = filter( e_ -> ( e_.attributes["style"]=="External" ), edges(graph0) )

  null_ext_edge_list = filter( e_ -> ( is_massless(e_.attributes["particle"]) ), ext_edge_list )
  @assert length(null_ext_edge_list) >= 2

  gauge_choice = Dict{Basic,Basic}()
  push!( gauge_choice, null_ext_edge_list[1].attributes["ref2_MOM"] => null_ext_edge_list[2].attributes["null_MOM"] )

  not1st_ext_edge_list = filter( e_ -> ( e_.attributes["mark"] != null_ext_edge_list[1].attributes["mark"] ), ext_edge_list )

  for edge in not1st_ext_edge_list
    if is_massive_fermion(edge.attributes["particle"])
      push!( gauge_choice, edge.attributes["ref2_MOM"] => Basic("barK$(edge.attributes["mark"])") )
    else
      push!( gauge_choice, edge.attributes["ref2_MOM"] => null_ext_edge_list[1].attributes["null_MOM"] )
    end # if
  end # for edge

  return gauge_choice

end # function generate_gauge_choice



######################################################################################################################
function generate_kin_relation( graph_list::Vector{GenericGraph}, gauge_choice::Dict{Basic,Basic} )::Dict{Basic,Basic}
######################################################################################################################

  graph0 = first( graph_list )

  @funs SP
  v0 = vertex_from_label( "graph property", graph0 )
  n_inc = v0.attributes["n_inc"]
  n_out = v0.attributes["n_out"]

  kin_relation = Dict{Basic,Basic}()

  ext_edge_list = filter( e_ -> ( e_.attributes["style"]=="External" ), edges(graph0) )
  @assert n_inc+n_out == length(ext_edge_list)

  null_ext_edge_list = filter( e_ -> ( is_massless(e_.attributes["particle"]) ), ext_edge_list )
  massive_ext_edge_list = filter( e_ -> ( is_massive(e_.attributes["particle"]) ), ext_edge_list )
  massive_fermion_ext_edge_list = filter( e_ -> ( is_massive_fermion(e_.attributes["particle"]) ), massive_ext_edge_list )

  # SP(Ki,Ki) = m_i^2
  for edge in ext_edge_list
    mom = edge.attributes["momentum"]
    push!( kin_relation, SP(mom,mom) => edge.attributes["particle"].mass^2 )
  end # for edge 

  # SP(ki,ki) = 0
  for edge in massive_ext_edge_list 
    null_mom = edge.attributes["null_MOM"]
    push!( kin_relation, SP(null_mom,null_mom) => Basic(0) )
  end # for edge

  # SP(barKi,barKi) = 0
  for edge in massive_fermion_ext_edge_list
    ref_mom = subs( edge.attributes["ref2_MOM"], gauge_choice... )
    push!( kin_relation, SP(ref_mom,ref_mom) => Basic(0) )
  end # for edge

  # massive momentum Ki can be decomposed Ki=ki+m^2/(2*ki.qi)*qi, so Ki.ki=m^2/2.
  for edge in massive_ext_edge_list
    mom = edge.attributes["momentum"]
    null_mom = edge.attributes["null_MOM"]
    push!( kin_relation, SP(null_mom,mom) => edge.attributes["particle"].mass^2/Basic(2) )
  end # for edge

  sorted_ext_edge_list = sort( ext_edge_list, by=edge_index )
  edge1 = sorted_ext_edge_list[1]
  edge2 = sorted_ext_edge_list[2]
  edgeN = sorted_ext_edge_list[n_inc+n_out]

  # Special Case:
  mom1 = edge1.attributes["momentum"]
  mom2 = edge2.attributes["momentum"]
  if n_inc == 1 && n_out == 2
    mass1 = edge1.attributes["particle"].mass
    mass2 = edge2.attributes["particle"].mass
    massN = edgeN.attributes["particle"].mass
    push!( kin_relation, SP(mom1,mom2) => (mass1^2+mass2^2-massN^2)/Basic(2) )
  elseif n_inc == 2 && n_out == 1
    massN = edgeN.attributes["mass"]
    push!( kin_relation, SP(mom1,mom2) => massN^2/Basic(2) )
  elseif n_inc == 2 && n_out > 1
    push!( kin_relation, SP(mom1,mom2) => Basic("shat/2") )
  else 
    printstyled( "n_inc & n_out exception!\n", color=:yellow )
    exit()
  end # if

  # Start to create verI's
  ver_index = 1;

  # two massive/massless momentums ki and pj (i<j), ki.pj=verI.
  for index1 in 1:(n_inc+n_out-1)
    edge1 = sorted_ext_edge_list[index1]
    for index2 in (index1+1):(n_inc+n_out)
      edge2 = sorted_ext_edge_list[index2]

      if index1 == 1 && index2 == 2
        continue
      end # if

      mom1 = edge1.attributes["momentum"]
      mom2 = edge2.attributes["momentum"]
      push!( kin_relation, SP(mom1,mom2) => Basic("ver$(ver_index)") )
      ver_index += 1

    end # for index2
  end # for index1

  # barKi.Kj=verI. 
  # NB: Here i==j is also included since barKi is related to Ki but not linearly.
  for edge1 in massive_fermion_ext_edge_list
    ref_mom1 = subs( edge1.attributes["ref2_MOM"], gauge_choice... )
    for edge2 in ext_edge_list
      mom2 = edge2.attributes["momentum"]

      push!( kin_relation, SP(ref_mom1,mom2) => Basic("ver$(ver_index)") )
      ver_index += 1
    end # for edge2
  end # for edge1

  sorted_massive_fermion_ext_edge_list = sort( massive_fermion_ext_edge_list, by=edge_index )
  # barKi.barKj = verI, (i < j) 
  for index1 in 1:(length(sorted_massive_fermion_ext_edge_list)-1)
    ref_mom1 = subs( sorted_massive_fermion_ext_edge_list[index1].attributes["ref2_MOM"], gauge_choice... )
    for index2 in (index1+1):length(sorted_massive_fermion_ext_edge_list)
      ref_mom2 = subs( sorted_massive_fermion_ext_edge_list[index2].attributes["ref2_MOM"], gauge_choice... )
      
      push!( kin_relation, SP(ref_mom1,ref_mom2) => Basic("ver$(ver_index)") )
      ver_index += 1
    end # for index2
  end # for index1

  # barKi.kj=verI, (i == j)
  for edge in massive_fermion_ext_edge_list
    mom = edge.attributes["momentum"]
    ref_mom = edge.attributes["ref2_MOM"]
    null_mom = edge.attributes["null_MOM"]
    push!( kin_relation, SP(ref_mom,null_mom) => subs( SP(ref_mom,mom), kin_relation... ) )
  end # for edge

  # barKi.kj=verI, (i != j)
  for edge1 in massive_fermion_ext_edge_list
    ref_mom1 = edge1.attributes["ref2_MOM"]
    for edge2 in massive_ext_edge_list
      null_mom2 = edge2.attributes["null_MOM"]
      
      push!( kin_relation, SP(ref_mom1,null_mom2) => Basic("ver$(ver_index)") )
      ver_index += 1
    end # for edge2
  end # for edge1

  # massive momentum Ki and its decomposed null momentum ki, massive/massless Kj (i != j), ki.Kj=verI
  for edge1 in massive_ext_edge_list
    null_mom1 = edge1.attributes["null_MOM"]
    for edge2 in ext_edge_list
      if edge_index(edge1) == edge_index(edge2)
        continue
      end # if
      mom2 = edge2.attributes["momentum"]

      push!( kin_relation, SP(null_mom1,mom2) => Basic("ver$(ver_index)") )
      ver_index += 1
    end # for edge2
  end # for edge1


  sorted_massive_ext_edge_list = sort( massive_ext_edge_list, by=edge_index )
  # ki.kj = verI, (i<j)
  for index1 in 1:(length(sorted_massive_ext_edge_list)-1)
    null_mom1 = edge1.attributes["null_MOM"]
    for index2 in (index1+1):length(sorted_massive_ext_edge_list)
      null_mom2 = edge2.attributes["null_MOM"]

      push!( kin_relation, SP(null_mom1,null_mom2) => Basic("ver$(ver_index)") )
      ver_index += 1
    end # for index2
  end # for index1

  @funs FF, JJ
  # FF(ki,kj) and JJ(ki,kj) (i<j)
  for index1 in 1:(length(null_ext_edge_list)-1)
    null_mom1 = null_ext_edge_list[index1].attributes["null_MOM"]
    for index2 in (index1+1):length(null_ext_edge_list)
      null_mom2 = null_ext_edge_list[index2].attributes["null_MOM"]

      push!( kin_relation, FF(null_mom1,null_mom2) => Basic("ver$(ver_index)") )
      push!( kin_relation, FF(null_mom2,null_mom1) => Basic("-ver$(ver_index)") )
      ver_index += 1
      push!( kin_relation, JJ(null_mom1,null_mom2) => Basic("ver$(ver_index)") )
      push!( kin_relation, JJ(null_mom2,null_mom1) => Basic("-ver$(ver_index)") )
      ver_index += 1
    end # for index2
  end # for index1

  @funs Den
  for graph in graph_list 

    int_edge_list = filter( e_ -> ( e_.attributes["style"]=="Internal" ), edges(graph) )
    for edge in int_edge_list
      mom = edge.attributes["momentum"]
      mass = edge.attributes["particle"].mass
      width = edge.attributes["particle"].width
  
      den_expr = Den(mom,mass,width)
      if den_expr != subs( den_expr, kin_relation... )
        continue
      end # if
  
      push!( kin_relation, Den(mom,mass,width) => Basic("ver$(ver_index)") )
      push!( kin_relation, Den(-mom,mass,width) => Basic("ver$(ver_index)") )
      ver_index += 1
    end # for edge

  end # for graph
  
  return kin_relation

end # function generate_kin_relation




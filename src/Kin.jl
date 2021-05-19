
####################################################################
"""
    generate_gauge_choice( graph_list::Vector{GenericGraph} )::Dict{Basic,Basic}

Automatically prepare a gauge choice for this process.
"""
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



########################################################################################
# \begin{quote}
# \hypertarget{header-n3}{%
# \section{\texorpdfstring{Kinematic variables for (\(n-2\))-body final
# state}{Kinematic variables for (n-2)-body final state}}\label{header-n3}}
# \end{quote}
# 
# Considering \(2\to (n-2)\) process, we set up momenta
# \(k_1,k_2,k_3,\dots,k_n\).
# 
# \textbf{According to RAJENDRA KUMAR, Phys.Rev. D2 (1970) 1902-1914}
# 
# For (\(n-2\))-body final state, there Eq. (2.6) shows the \((3(n-2)-4)\)
# independent variables are
# 
# \[s_{i} = \left( k_1 -\sum_{j=3}^{i+2} k_j \right)^2, \quad \text{for}~~i=1,\dots,n-3\label{eq1}\]
# 
# \[s_{i+n-3} = \left( k_1+k_2-\sum_{j=3}^{i+2}k_j\right)^2, \quad \text{for}~~i=1,\dots,n-4\label{eq2}\]
# 
# \[s_{i+2n-7} = \left(\sum_{j=3}^{i+3} k_j\right)^2,\quad \text{for}~~i=1,\dots,n-4\label{eq3}\]
# 
# \[s =(k_1+k_2)^2\]
# 
# And those variables may depends on above variables but not linearly
# 
# \[(k_i\cdot k_j), \quad \text{for}~~i<j,~\text{and}~i,j\in\lbrace 4,\dots,n-1\rbrace\]
# 
# which can be used as \((n-4)(n-5)/2\) independent variables.
# 
# However, in total we have \(n(n-1)/2\) scalar products from
# \(k_1,k_2,\dots,k_n\).
# 
# \begin{enumerate}
# \def\labelenumi{\arabic{enumi}.}
# \item
#   We can take equations in \((\ref{eq1})\) by expansion
# 
#   \[s_i = k_1^2-2k_1\cdot\sum_{j=3}^{i+2}k_j +\left(\sum_{j=3}^{i+2}k_j\right)^2,\quad \text{for}~~i=1,\dots,n-3\]
# 
#   Then we can make subtraction \(s_i-s_{i-1}\) to obtain
# 
#   \[(k_1\cdot k_3) = \frac{1}{2}\left(k_1^2+k_3^2-s_1\right),\\(k_1\cdot k_4) = \frac{1}{2}\left(s_{1+2n-7}-k_3^2-s_2+s_1\right),\\(k_1\cdot k_{i+2}) = \frac{1}{2}\left(s_{i-1+2n-7}-s_{i-2+2n-7}-s_i+s_{i-1}\right),\quad \text{for}~~i=3,\dots,n-3\]
# 
#   And
# 
#   \[(k_1\cdot k_n) = k_1\cdot(k_1+k_2-k_3-\cdots-k_{n-1})\]
# 
#   Therefore we got \((n-2)\) scalar products
#   \((k_1\cdot k_3),~ (k_1\cdot p_4),~ \dots,~ (k_1\cdot k_n)\).
# \end{enumerate}
# 
# \begin{enumerate}
# \def\labelenumi{\arabic{enumi}.}
# \item
#   We can take equations in \((\ref{eq1})\) again by
# 
#   \[s_i = \left(k_2-\sum_{j=i+3}^n k_j\right)^2,\quad \text{for}~~i=1,\dots,n-3\]
# 
#   Then we can make subtraction \(s_i-s_{i+1}\) to obtain
# 
#   \[(k_2\cdot k_{i+3})=\frac{1}{2}(s_{i+n-3}-s_{i+1+n-3}-s_i+s_{i+1}),\quad \text{for}~~i=1,\dots,n-5\\(k_2\cdot k_{n-1})=\frac{1}{2}(s_{n-4+(n-3)}-k_n^2-s_i+s_{i+1}),\\(k_2\cdot k_n)=\frac{1}{2}(k_2^2+k_n^2-s_{n-3}),\]
# 
#   And according to \(s_{1+n-3}=(k_1+k_2-k_3)^2\) from \((\ref{eq2})\),
#   we can have
# 
#   \[(k_2\cdot k_3)=\frac{1}{2}(k_1^2+k_2^2+k_3^2+2k_1\cdot k_3-2k_1\cdot k_3-s_{1+n-3})\]
# 
#   Therefore, we got \(n-2\) scalar products
#   \((k_2\cdot k_3),~(k_2\cdot k_4),~\dots,~(k_2\cdot k_{n})\).
# \end{enumerate}
# 
# \begin{enumerate}
# \def\labelenumi{\arabic{enumi}.}
# \item
#   We can also take equations in \((\ref{eq3})\) by
# 
#   \[s_{i+2n-7} = \left(k_3+\sum_{j=4}^{i+3} k_j\right)^2,\quad \text{for}~~ i=1,\dots,n-4\]
# 
#   Then subtraction \(s_{i+2n-7}-s_{i-1+2n-7}\) can give us
# 
#   \[(k_3\cdot k_{i+3})=\frac{1}{2}\left[s_{i+2n-7}-s_{i-1+2n-7}-\left(\sum_{j=4}^{i+3}k_j\right)^2-\left(\sum_{j=4}^{i+2}k_j\right)^2\right],\quad \text{for}~~i=2,\dots,n-4\\(k_3\cdot k_4)=\frac{1}{2}(s_{1+(2n-7)}-k_3^2-k_4^2).\]
# 
#   And
# 
#   \[(k_3\cdot k_n) = k_3\cdot(k_1+k_2-k_3-\cdots-k_{n-1})\]
# 
#   Therefore, we got \(n-3\) scalar products
#   \((k_3\cdot k_4),~(k_3\cdot k_5),~\dots,~(k_3\cdot k_n)\).
# \item
#   We can take equations in \((\ref{eq2})\) by
# 
#   \[s_{i+(n-3)} = \left( k_n+\sum_{j=i+3}^{n-1}k_j\right)^2,\quad \text{for}~~i=1,\dots,n-4\]
# 
#   Then subtraction \(s_{i+(n-3)}-s_{i+1+(n-3)}\) can give us
# 
#   \[(k_n\cdot k_{i+3}) = \frac{1}{2}\left[s_{i+(n-3)}-s_{i+1+(n-3)}-\left(\sum_{j=i+3}^{n-1}k_j\right)^2+\left(\sum_{j=i+4}^{n-1}k_j\right)^2\right],\quad \text{for}~~i=1,\dots,n-5\\(k_n\cdot k_{n-1}) =\frac{1}{2}\left(s_{n-4+(n-3)}-k_n^2-k_{n-1}^2\right)\]
# 
#   Therefore we got \(n-4\) scalar products
#   \((k_n\cdot k_4),~(k_n\cdot k_5),~\dots,~(k_n\cdot k_{n-1})\).
# \end{enumerate}
# 
# \begin{enumerate}
# \def\labelenumi{\arabic{enumi}.}
# \item
#   Also we have \((k_1\cdot k_2) = \frac{1}{2} (s-k_1^2-k_2^2)\).
# \end{enumerate}
# 
# Finally in total we have
# \((n-2)+(n-2)+(n-3)+(n-4)+1+(n-4)(n-5)/2=n(n-1)/2\) solved scalar
# products, which is consistent with the total number.
# 
# \textbf{For the 1-(n-1) decay mode, we can simply change the sign of
# \(k_2\).}
#
#########################################################################################################################
"""
    generate_kin_relation( n_inc::Int64, n_out::Int64, mom::Vector{Basic}, mass2::Vector{Basic} )::Dict{Basic,Basic}

Generate the kinematic relations, e.g. Mandelstam variables, according to the external fields.
"""
function generate_kin_relation( n_inc::Int64, n_out::Int64, mom::Vector{Basic}, mass2::Vector{Basic} )::Dict{Basic,Basic}
#########################################################################################################################

  @funs SP
  @vars shat

  nn = n_inc+n_out
  k2_sign = n_inc == 2 ? (+1) : (-1)
  half = Basic(1)/Basic(2)

  kin_relation = Dict{Basic,Basic}()

  ver_index_pre = 3*(nn-2)-4-1 # pre-occupied index for 3(n-2)-4 variable and one of them is shat
  ver_index = ver_index_pre + 1 
  for ii in 4:(nn-1)
    for jj in (ii+1):(nn-1)
      push!( kin_relation, make_SP(mom[ii],mom[jj]) => Basic("ver$(ver_index)") )
      ver_index += 1
    end # for jj
  end # for ii

  # on-shell conditions 
  for ii in 1:nn
    push!( kin_relation, make_SP(mom[ii],mom[ii]) => mass2[ii] )
  end # for ii


  # (k_1\cdot k_2) = \frac{1}{2} (s-k_1^2-k_2^2)  
  if n_out == 1
    push!( kin_relation, make_SP(mom[1],mom[2]) => half*( mass2[3] - mass2[1] - mass2[2] ) )
  else 
    push!( kin_relation, make_SP(mom[1],mom[2]) => k2_sign*half*( shat - mass2[1] - mass2[2] ) )
  end # if

  # (k_1\cdot k_3) = \frac{1}{2}\left(k_1^2+k_3^2-s_1\right),
  @vars ver1 # s_1
  if n_out == 1
    push!( kin_relation, make_SP(mom[1],mom[3]) => half*( mass2[1] + mass2[3] - mass2[2] ) )
  else 
    push!( kin_relation, make_SP(mom[1],mom[3]) => half*( mass2[1] + mass2[3] - ver1 ) )
  end # if

  @vars ver2 # s_2
  if nn >= 4
    if nn == 4
      # (k_1\cdot k_4) = \frac{1}{2}\left(shat-k_3^2-k_2^2+s_1\right),
      push!( kin_relation, make_SP(mom[1],mom[4]) => half*( shat - mass2[3] - mass2[2] + ver1 ) )
    else 
      # (k_1\cdot k_4) = \frac{1}{2}\left(s_{1+2n-7}-k_3^2-s_2+s_1\right),
      push!( kin_relation, make_SP(mom[1],mom[4]) => half*( Basic("ver$(1+2*nn-7)") - mass2[3] - ver2 + ver1 ) )
    end # if
  end # if

  # (k_1\cdot k_{i+2}) = \frac{1}{2}\left(s_{i-1+2n-7}-s_{i-2+2n-7}-s_i+s_{i-1}\right),\quad \text{for}~~i=3,\dots,n-3
  for ii in 3:(nn-3)
    push!( kin_relation, make_SP(mom[1],mom[ii+2]) => half*Basic("ver$(ii-1+2*nn-7) - ver$(ii-2+2*nn-7) - ver$(ii) + ver$(ii-1)") )
  end # for ii

  if nn > 4
    # (k_1\cdot k_n) = k_1\cdot(k_1+k_2-k_3-\cdots-k_{n-1})
    rhs = mass2[1] + k2_sign*make_SP(mom[1],mom[2]) # k_1\cdot k_1 + k_1\cdot k_2
    for ii = 3:(nn-1)
      rhs += (-1)*make_SP(mom[1],mom[ii])
    end # for ii
    rhs = subs( rhs, kin_relation... )
    push!( kin_relation, make_SP(mom[1],mom[nn]) => rhs )
  end # if



  # (k_2\cdot k_{i+3})=\frac{1}{2}(s_{i+n-3}-s_{i+1+n-3}-s_i+s_{i+1}),\quad \text{for}~~i=1,\dots,n-5
  for ii in 1:(nn-5)
    push!( kin_relation, make_SP(mom[2],mom[ii+3]) => k2_sign*half*Basic("ver$(ii+nn-3) - ver$(ii+1+nn-3) - ver$(ii) + ver$(ii+1)") )
  end # for ii

  # (k_2\cdot k_n)=\frac{1}{2}(k_2^2+k_n^2-s_{n-3}),
  if n_out == 1
    push!( kin_relation, make_SP(mom[2],mom[3]) => half*(mass2[3]+mass2[2]-mass2[1]) )
  else 
    push!( kin_relation, make_SP(mom[2],mom[nn]) => k2_sign*half*(mass2[2]+mass2[nn]-Basic("ver$(nn-3)")) )
  end # if

  if nn >= 4
    if nn == 4
      # (k_2\cdot k_3) = k_1\cdot k_2 + k_2^2 - k_2\cdot k_4
      push!( kin_relation, make_SP(mom[2],mom[3]) => k2_sign*(expand∘subs)( k2_sign*make_SP(mom[1],mom[2]) + mass2[2] - k2_sign*make_SP(mom[2],mom[4]), kin_relation... ) )
    else 
      # (k_2\cdot k_{n-1})=\frac{1}{2}(s_{n-4+(n-3)}-k_n^2-s_{n-4}+s_{n-3}),
      push!( kin_relation, make_SP(mom[2],mom[nn-1]) => k2_sign*half*Basic("ver$(nn-4+nn-3) - $(mass2[nn]) - ver$(nn-4) + ver$(nn-3)") )
    end # if
  end # if

  if nn == 4
    # (k_2\cdot k_3)=\frac{1}{2}(k_1^2+k_2^2+k_3^2+2k_1\cdot k_2-2k_1\cdot k_3-k_4^2)
    push!( kin_relation, make_SP(mom[2],mom[3]) => k2_sign*half*(expand∘subs)( mass2[1] + mass2[2] + mass2[3] + k2_sign*2*make_SP(mom[1],mom[2]) - 2*make_SP(mom[1],mom[3]) - mass2[4], kin_relation... ) )
  elseif n_out == 1
    push!( kin_relation, make_SP(mom[2],mom[3]) => half*( mass2[3] + mass2[2] - mass2[1] ) )
  else 
    # (k_2\cdot k_3)=\frac{1}{2}(k_1^2+k_2^2+k_3^2+2k_1\cdot k_2-2k_1\cdot k_3-s_{1+n-3})
    push!( kin_relation, make_SP(mom[2],mom[3]) => k2_sign*half*(expand∘subs)( mass2[1] + mass2[2] + mass2[3] + k2_sign*2*make_SP(mom[1],mom[2]) - 2*make_SP(mom[1],mom[3]) - Basic("ver$(1+nn-3)"), kin_relation... ) )
  end # if


  if nn >= 4
    if nn == 4
      # (k_3\cdot k_4)=\frac{1}{2}(shat-k_3^2-k_4^2).
      push!( kin_relation, make_SP(mom[3],mom[4]) => half*( shat - mass2[3] - mass2[4] ) )
    else 
      # (k_3\cdot k_4)=\frac{1}{2}(s_{1+(2n-7)}-k_3^2-k_4^2).
      push!( kin_relation, make_SP(mom[3],mom[4]) => half*( Basic("ver$(1+2*nn-7)") - mass2[3] - mass2[4] ) )
    end # if
  end # if


  # (k_3\cdot k_{i+3})=\frac{1}{2}\left[s_{i+2n-7}-s_{i-1+2n-7}-\left(\sum_{j=4}^{i+3}k_j\right)^2-\left(\sum_{j=4}^{i+2}k_j\right)^2\right],\quad \text{for}~~i=2,\dots,n-4
  for ii in 2:(nn-4)
    rhs = Basic("ver$(ii+2*nn-7) - ver$(ii-1+2*nn-7)") - mass2[ii+3]
    for jj in 4:(ii+2)
      rhs += (-2)*make_SP(mom[jj],mom[ii+3])
    end # for jj
    push!( kin_relation, make_SP(mom[3],mom[ii+3]) => half*(expand∘subs)( rhs, kin_relation... ) )
  end # for ii

  # (k_3\cdot k_n) = k_3\cdot(k_1+k_2-k_3-\cdots-k_{n-1})
  rhs = make_SP(mom[1],mom[3])+k2_sign*make_SP(mom[2],mom[3])-mass2[3]
  for ii in 4:(nn-1)
    rhs += (-1)*make_SP(mom[3],mom[ii])
  end # for ii
  if n_out == 1
    push!( kin_relation, make_SP(mom[3],mom[3]) => mass2[3] )
  else
    push!( kin_relation, make_SP(mom[3],mom[nn]) => (expand∘subs)( rhs, kin_relation... ) )
  end # if


  if nn >= 4
    if nn == 4 
      # (k_n\cdot k_{n-1}) =\frac{1}{2}\left(shat-k_n^2-k_{n-1}^2\right)
      push!( kin_relation, make_SP(mom[nn-1],mom[nn]) => half*( shat - mass2[nn] - mass2[nn-1] ) )
    else 
      # (k_n\cdot k_{n-1}) =\frac{1}{2}\left(s_{n-4+(n-3)}-k_n^2-k_{n-1}^2\right)
      push!( kin_relation, make_SP(mom[nn-1],mom[nn]) => half*( Basic("ver$(nn-4+nn-3)") - mass2[n] - mass2[nn-1] ) )
    end # if
  end # if

  # (k_n\cdot k_{i+3}) = \frac{1}{2}\left[s_{i+(n-3)}-s_{i+1+(n-3)}-\left(\sum_{j=i+3}^{n-1}k_j\right)^2+\left(\sum_{j=i+4}^{n-1}k_j\right)^2\right],\quad \text{for}~~i=1,\dots,n-5
  for ii in 1:(nn-5)
    rhs = Basic("ver$(ii+nn-3) - ver$(ii+1+nn-3)") - mass2[ii+3]
    for jj in (ii+4):(nn-1)
      rhs += (-2)*make_SP(mom[ii+3],mom[jj])
    end # for jj
    push!( kin_relation, make_SP(mom[ii+3],mom[nn]) => half*(expand∘subs)( rhs, kin_relation... ) )
  end # for ii


  return kin_relation

end # function generate_kin_relation

########################################################################################
"""
    generate_kin_relation( graph_list::Vector{GenericGraph} )::Dict{Basic,Basic}

Generate the kinematic relations for this processes including the internal non-loop propagators.
"""
function generate_kin_relation( graph_list::Vector{GenericGraph} )::Dict{Basic,Basic}
########################################################################################

  graph0 = first( graph_list )

  v0 = vertex_from_label( "graph property", graph0 )
  n_inc = v0.attributes["n_inc"]
  n_out = v0.attributes["n_out"]
  nn = n_inc+n_out

  ext_edge_list = filter( e_ -> ( e_.attributes["style"]=="External" ), edges(graph0) )
  @assert n_inc+n_out == length(ext_edge_list)

  sorted_ext_edge_list = sort( ext_edge_list, by=e_->e_.attributes["mark"] )
  mom = map( e_->e_.attributes["momentum"], sorted_ext_edge_list ) 
  mass2 = map( e_->e_.attributes["particle"].mass^2, sorted_ext_edge_list )
  # here mom and mass2 in fact are Vector{Basic}, but for later convenience we do not explicitly show it.

  @info "external momenta: $mom"

  kin_relation = generate_kin_relation( n_inc, n_out, mom, mass2 )


  mom_n = Basic(0)
  for ii in 1:n_inc
    mom_n += mom[ii]
  end # for ii
  for ii in (n_inc+1):(nn-1)
    mom_n += (-1)*mom[ii]
  end # for ii

  @info "momentum conservation: $(mom[nn]) = $(mom_n)"


  @funs Den
  den_set = Set{Basic}()
  for g in graph_list
    int_edge_list = filter( e_->e_.attributes["style"] == "Internal", edges(g) )
    for edge in int_edge_list
      den_mom = subs( edge.attributes["momentum"], mom[nn] => mom_n )
      den_mass = edge.attributes["particle"].mass
      den_width = edge.attributes["particle"].width
      push!( den_set, Den(den_mom,den_mass,den_width) )
    end # for edge
  end # for g

  # we know at most n(n-1)/2-1 verI's have been occupied
  ver_index_pre = nn*(nn-1)/Basic(2) - 1
  ver_index = ver_index_pre + 1
  for one_den in den_set
    push!( kin_relation, one_den => Basic("ver$(ver_index)") )

    arg_list = get_args(one_den)
    den_mom = arg_list[1]
    den_mass = arg_list[2]
    den_width = arg_list[3]
    push!( kin_relation, Den(expand(-den_mom),den_mass,den_width) => Basic("ver$(ver_index)") )

    ver_index += 1 
  end # for one_den

  for g in graph_list
    int_edge_list = filter( e_->e_.attributes["style"] == "Internal", edges(g) )
    for edge in int_edge_list
      den_mom = edge.attributes["momentum"]
      den_mass = edge.attributes["particle"].mass
      den_width = edge.attributes["particle"].width

      subs_den_mom = subs( den_mom, mom[nn] => mom_n )
      if den_mom != subs_den_mom
        push!( kin_relation, Den(den_mom,den_mass,den_width) => subs(Den(subs_den_mom,den_mass,den_width),kin_relation...) )
        push!( kin_relation, Den(expand(-den_mom),den_mass,den_width) => subs(Den(subs_den_mom,den_mass,den_width),kin_relation...) )
      end # if

    end # for edge
  end # for g

  return kin_relation

end # function generate_kin_relation

#################################################################################
"""
    generate_ext_mom_list( graph_list::Vector{GenericGraph} )::Vector{Basic}

Generate the list of external momenta according to this process.
"""
function generate_ext_mom_list( graph_list::Vector{GenericGraph} )::Vector{Basic}
#################################################################################

  graph0 = first( graph_list )

  v0 = vertex_from_label( "graph property", graph0 )
  n_inc = v0.attributes["n_inc"]
  n_out = v0.attributes["n_out"]

  ext_edge_list = filter( e_ -> ( e_.attributes["style"]=="External" ), edges(graph0) )
  @assert n_inc+n_out == length(ext_edge_list)

  sorted_ext_edge_list = sort( ext_edge_list, by=e_->e_.attributes["mark"] )
  ext_mom_list = map( e_->e_.attributes["momentum"], sorted_ext_edge_list ) 

  return ext_mom_list

end # function generate_ext_mom_list







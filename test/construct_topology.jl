using AbstractAlgebra
using Combinatorics
using JLD2
using SymEngine
using AmpTools

###########################################
mutable struct DenTop 
  den_list::Vector{Basic} 
  coeff_mat::Matrix{Rational}  
end # struct DenTop 
###########################################
mutable struct TopSet 
  n_sp::Int64 
  incomplete_dentop_list::Vector{DenTop} 
  complete_dentop_list::Vector{DenTop} 
  den_coeff_mat_dict::Dict{Basic,Matrix{Rational}} 
end # struct TopSet 
###########################################




#######################
function is_complete(
    topset::TopSet 
)::Bool
#######################

  return isempty(topset.incomplete_dentop_list)

end # function is_complete



##################################
function get_coeff_mat(
    one_den::Basic, # Den(...)
    sp_dict::Dict{Basic,Basic} 
)::Matrix{Rational}
##################################

  sp_var_list = (collect∘values)(sp_dict)

  one_mom = (expand∘first∘get_args)(one_den)
  mom2 = subs( make_SP(one_mom,one_mom), sp_dict )

  n_sp = length(sp_var_list)
  coeff_mat = zeros( Rational, 1, n_sp )
  for sp_index in 1:n_sp
    sp_coeff = SymEngine.coeff( mom2, sp_var_list[sp_index] )
    coeff_mat[1,sp_index] = Rational( parse( Int64, string(sp_coeff) ) )
  end # for sp_index

  return coeff_mat

end # function get_coeff_mat

#################################
function make_dentop()::DenTop
#################################

  return DenTop( Vector{Basic}(), zeros(Rational,0,0) )

end # function make_dentop


#################################
function is_empty( 
    dentop::DenTop 
)::Bool 
#################################

  return (iszero∘length)( dentop.den_list )

end # function is_empty


##############################
function make_dentop( 
    den_list::Vector{Basic}, 
    den_coeff_mat_dict::Dict{Basic,Matrix{Rational}}, 
    sp_dict::Dict{Basic,Basic} 
)::DenTop
##############################

  row_list = map( x->get_coeff_mat(x,sp_dict), den_list )
  coeff_mat = reduce( vcat, row_list )

  return DenTop( den_list, coeff_mat )

end # function make_dentop


####################################
function make_dentop(
    dentop::DenTop,
    extra_den_list::Vector{Basic}, 
    den_coeff_mat_dict::Dict{Basic,Matrix{Rational}} 
)::DenTop
####################################

  if isempty(extra_den_list)
    return dentop
  end # if

  new_den_list = vcat( dentop.den_list, extra_den_list )
  extra_coeff_mat_list = map( x->den_coeff_mat_dict[x], extra_den_list )
  extra_coeff_mat = reduce( vcat, extra_coeff_mat_list )
  new_coeff_mat = vcat( dentop.coeff_mat, extra_coeff_mat )

  return DenTop( new_den_list, new_coeff_mat )

end # function make_dentop



####################################
# used to compare two dentop_collect to unique the multiverse
#  Den*Den*...+Den*Den*...+....
function get_dentop_collect_signature(
    dentop_collect::Vector{DenTop} 
)::Basic
####################################

  return (sum∘map)( x->reduce(*,x.den_list), dentop_collect )

end # function get_dentop_collect_signature







#########################################################
function intersection( 
    dentop1::DenTop, 
    dentop2::DenTop 
)::Tuple{Vector{Basic},Vector{Basic},Vector{Basic}}
#########################################################

  com_den_list = intersect( dentop1.den_list, dentop2.den_list )
  rem1_den_list = setdiff( dentop1.den_list, com_den_list )
  rem2_den_list = setdiff( dentop2.den_list, com_den_list )

  return com_den_list, rem1_den_list, rem2_den_list

end # function intersection


############################################
function gen_sp_dict(
    qi_list::Vector{Basic},
    indep_mom_list::Vector{Basic} 
)::Dict{Basic,Basic} 
############################################

  n_loop = length(qi_list)

  sp_index = 1::Int64
  sp_dict = Dict{Basic,Basic}()

  for i1 in 1:n_loop, i2 in i1:n_loop
    qi = qi_list[i1] 
    qj = qi_list[i2]
    sp_var = Basic("sp$(sp_index)")
    push!( sp_dict, make_SP(qi,qj) => sp_var )
    sp_index += 1
  end # for qi_index, qj_index

  for one_mom in indep_mom_list, qi_index in 1:n_loop
    qi = qi_list[qi_index]
    sp_var = Basic( "sp$(sp_index)" )
    push!( sp_dict, make_SP(one_mom,qi) => sp_var )
    sp_index += 1
  end # for one_mom

  return sp_dict

end # function gen_sp_dict

######################
# dentop1 is child of dentop2 ?
function is_child( 
    dentop1::DenTop, 
    dentop2::DenTop
)::Bool
######################

  return dentop1.den_list ⊆ dentop2.den_list

end # function is_child


######################
# dentop1 is parent of dentop2 ?
function is_parent( 
    dentop1::DenTop, 
    dentop2::DenTop
)::Bool
######################

  return dentop2.den_list ⊆ dentop1.den_list

end # function is_parent

###############################################
function find_first_parent(
    dentop::DenTop, 
    dentop_collect::Vector{DenTop} 
)::Union{Int64,Nothing}
###############################################

  # Search for the parent of dentop
  pos = findfirst( x -> is_parent(x,dentop), dentop_collect )
  return pos

end # function find_first_parent


###############################################
function find_all_children(
    dentop::DenTop, 
    dentop_collect::Vector{DenTop} 
)::Vector{Int64}
###############################################

  # Search for the children of dentop
  pos_list = findall( x -> is_child(x,dentop), dentop_collect )
  return pos_list

end # function find_all_children


#############################################
function get_superior_dentop_collect(
    dentop_collect::Vector{DenTop}  
)::Vector{DenTop} 
#############################################

  new_dentop_collect = Vector{DenTop}()
  for dentop in dentop_collect

    pos = find_first_parent( dentop, new_dentop_collect )
    if pos != nothing 
      # Found parent and check next.
      continue
    end # if

    # Remove all children and then insert as parent.
    pos_list = find_all_children( dentop, new_dentop_collect )
    for pos in pos_list
      new_dentop_collect[pos] = make_dentop()
    end # for pos

    push!( new_dentop_collect, dentop )
  end # for dentop 
  new_dentop_collect = filter( !is_empty, new_dentop_collect )

  return new_dentop_collect 

end # function get_superior_dentop_collect





###############################
function make_parent_dentop( 
    n_sp::Int64, 
    dentop1::DenTop,
    dentop2::DenTop,
    den_coeff_mat_dict::Dict{Basic,Matrix{Rational}} 
)::Union{DenTop,Nothing} 
###############################

  com_den_list, rem1_den_list, rem2_den_list = intersection( dentop1, dentop2 )
  n_com = length(com_den_list)
  n_rem1 = length(rem1_den_list)
  n_rem2 = length(rem2_den_list)

  if n_com+n_rem1+n_rem2 > n_sp
    # cannot make common parent dentop
    return nothing
  end # if

  com_row_list = map( x -> den_coeff_mat_dict[x], com_den_list )
  com_coeff_mat = isempty(com_row_list) ? zeros(Rational,0,n_sp) : reduce( vcat, com_row_list )

  rem1_row_list = map( x -> den_coeff_mat_dict[x], rem1_den_list )
  rem1_coeff_mat = isempty(rem1_row_list) ? zeros(Rational,0,n_sp) : reduce( vcat, rem1_row_list )

  rem2_row_list = map( x -> den_coeff_mat_dict[x], rem2_den_list )
  rem2_coeff_mat = isempty(rem2_row_list) ? zeros(Rational,0,n_sp) : reduce( vcat, rem2_row_list )

  the_coeff_mat = vcat( com_coeff_mat, rem1_coeff_mat, rem2_coeff_mat )

  # the common parent should contain all of their den's, so they should be independent.
  if rank(the_coeff_mat) < n_com+n_rem1+n_rem2
    return nothing
  end # if

  parent_dentop = make_dentop( dentop1, rem2_den_list, den_coeff_mat_dict )

  return parent_dentop

end # function make_parent_dentop



###############################
function make_parent_dentop( 
    n_sp::Int64, 
    dentop_list::Vector{DenTop},
    den_coeff_mat_dict::Dict{Basic,Matrix{Rational}} 
)::Union{DenTop,Nothing} 
###############################

  n_dentop = length(dentop_list)

  parent_dentop = first(dentop_list)
  for index in 2:n_dentop
    parent_dentop = make_parent_dentop( n_sp, parent_dentop, dentop_list[index], den_coeff_mat_dict )
  end # for index

  return parent_dentop

end # function make_parent_dentop





##################################
function get_union_den_list( 
    dentop_list::Vector{DenTop}
)::Vector{Basic}
##################################

  union_den_list = union( map( x->x.den_list, dentop_list )... ) 
  return union_den_list

end # function get_union_den_list

########################################
function get_coeff_mat_rank(
    den_list::Vector{Basic},
    den_coeff_mat_dict::Dict{Basic,Matrix{Rational}} 
)::Int64
########################################

  row_list = map( x->den_coeff_mat_dict[x], den_list ) 
  the_rank = (rank∘reduce)( vcat, row_list )
  return the_rank

end # function get_coeff_mat_rank

###################################################
function is_fullrank(
    den_list::Vector{Basic}, 
    den_coeff_mat_dict::Dict{Basic,Matrix{Rational}} 
)::Bool
###################################################

  return get_coeff_mat_rank(den_list,den_coeff_mat_dict) == length(den_list)

end # function is_fullrank

##############################
function find_max_cover(
    n_sp::Int64, 
    dentop_collect::Vector{DenTop}, 
    den_coeff_mat_dict::Dict{Basic,Matrix{Rational}} 
)::Vector{Vector{Int64}}
##############################

  n_dentop = length(dentop_collect)
  #prev_indices_list = [ Int64[] ]
  prev_indices_list = [ Int64[x] for x in 1:n_dentop ]
  for n_choice in 2:n_dentop
    #println( "============" )
    #@show n_choice

    #indices_list = (collect∘powerset)( 1:n_dentop, n_choice, n_choice )
    #@show length(indices_list)
    #indices_list = filter( x->any(y->y⊆x,prev_indices_list), indices_list )
    #@show length(indices_list)

    indices_list = Vector{Vector{Int64}}()
    for one_indices in prev_indices_list
      rem_indices = setdiff(collect(1:n_dentop),one_indices)
      for one_index in rem_indices
        push!( indices_list, (sort∘vcat)( one_indices, [one_index] ) )
      end # for one_index
    end # for one_indices
    unique!(indices_list)
    #@show length(indices_list)

    dentop_subset_list = map( x->dentop_collect[x], indices_list )
    #dentop_subset_list = (collect∘powerset)( dentop_collect, n_choice, n_choice )
    dentop_subset_union_list = map( get_union_den_list, dentop_subset_list )

    # number of allowed by simply counting den_list
    allowed_pos_list = findall( x -> length(x) <= n_sp, dentop_subset_union_list ) 
    n_allowed = length(allowed_pos_list) 
    #@show n_allowed 

    allowed_indices_list = indices_list[allowed_pos_list]
    dentop_subset_list = dentop_subset_list[allowed_pos_list]
    dentop_subset_union_list = dentop_subset_union_list[allowed_pos_list]
    fullrank_pos_list = findall( x->is_fullrank(x,den_coeff_mat_dict), dentop_subset_union_list )
    n_fullrank = length(fullrank_pos_list)
    #@show n_fullrank

    if n_fullrank == 0
      return prev_indices_list
    end # if

    fullrank_indices_list = allowed_indices_list[fullrank_pos_list]
    prev_indices_list = fullrank_indices_list
    #@show length(prev_indices_list) #prev_indices_list

  end # for n_choice

  error("Is this possible? All dentop_collect can be covered by one topology.")

end # function find_max_cover


############################
# greedy algorithm
function greedy( 
    indices_list::Vector{Vector{Int64}}
)::Vector{Vector{Int64}}
############################

  universe = union( indices_list... )

  new_indices_list = Vector{Vector{Int64}}()
  while !isempty(universe)
    n_cover, index = findmax( x->(length∘intersect)(x,universe), indices_list )
    the_indices = indices_list[index]

    push!( new_indices_list, the_indices )
    setdiff!( universe, the_indices )
    indices_list = indices_list[ setdiff(1:length(indices_list),index) ]
  end # while

  return new_indices_list

end # function greedy




##############################
function update_dentop_list(
    n_sp::Int64, 
    dentop_list::Vector{DenTop}, 
    den_coeff_mat_dict::Dict{Basic,Matrix{Rational}} 
)::Tuple{ Vector{DenTop}, Vector{DenTop} }
##############################

  max_cover_indices_list = find_max_cover( n_sp, dentop_list, den_coeff_mat_dict ) 
  @show max_cover_indices_list

  greedy_indices_list = greedy(max_cover_indices_list)

  complete_dentop_list = Vector{DenTop}()
  for indices in greedy_indices_list
    one_complete_dentop = make_parent_dentop( n_sp, dentop_list[indices], den_coeff_mat_dict )
    push!( complete_dentop_list, one_complete_dentop )
  end # for indices

  # remove completed dentop's
  rem_indices = setdiff( 1:length(dentop_list), (sort∘union)( greedy_indices_list... ) )
  incomplete_dentop_list = dentop_list[rem_indices]

  return incomplete_dentop_list, complete_dentop_list

end # function update_dentop_list



#############################
function finalize_complete_dentop(
    dentop::DenTop, 
    indep_mom_list::Vector{Basic}, 
    sp_dict::Dict{Basic,Basic} 
)::DenTop
#############################

  n_sp = length(sp_dict)
  rank_coeff_mat = rank(dentop.coeff_mat) 

  n_loop = (length∘unique∘filter)( x -> (first∘string)(x) == 'q', free_symbols(dentop.den_list) )

  @vars q1, q2, q3
  if n_loop == 1
    loop_mom_list = [q1]
  elseif n_loop == 2
    loop_mom_list = [q1,q2,q1+q2]
  else 
    error("Exception")
    loop_mom_list = Baisc[]
  end # if

  if rank_coeff_mat < n_sp
    @funs Den

    for ext_mom in vcat(zero(Basic),indep_mom_list), qi in loop_mom_list, sign in [-1,+1]
      loop_mom = expand( qi+sign*ext_mom )
      trial_den = Den( loop_mom, 0, 0 )

      row_mat = get_coeff_mat( trial_den, sp_dict )
      new_coeff_mat = vcat( dentop.coeff_mat, row_mat )
      if rank(new_coeff_mat) > rank(dentop.coeff_mat)
        dentop.coeff_mat = new_coeff_mat
        push!( dentop.den_list, trial_den )

        if rank(dentop.coeff_mat) == n_sp
          break;
        end # if
      end # if
  
    end # for ext_mom

    @assert n_sp == rank(dentop.coeff_mat) == length(dentop.den_list)
  end # if

  return dentop

end # function finalize_complete_dentop









###########################
function main()::Nothing
###########################


  dir = "b_g_TO_Wminus_t_2Loop_amplitudes"
  root, dirs, files = (first∘collect∘walkdir)(dir)
  file_list = filter( s->endswith(s,".jld2"), files )
  @show length(file_list)

  file = jldopen( "$(dir)/$(first(file_list))", "r" )
  n_loop = read( file, "n_loop" )
  ext_mom_list = (to_Basic∘read)( file, "ext_mom_list" )
  close( file )

  qi_list = map( index->Basic("q$(index)"), 1:n_loop )
  indep_mom_list = ext_mom_list[1:end-1] 
  sp_dict = gen_sp_dict( qi_list, indep_mom_list )
  n_sp = length(sp_dict) # equal to number den's for complete topology.

  #------
  dentop_collect = Vector{DenTop}() 
  den_coeff_mat_dict = Dict{Basic,Matrix{Rational}}()
  for one_file in file_list
    file = jldopen( "$(dir)/$(one_file)", "r" )
    den_list = (to_Basic∘read)( file, "loop_den_list" )
    close( file )

    #-------
    for one_den in den_list
      if haskey( den_coeff_mat_dict, one_den )
        continue
      end # if
      push!( den_coeff_mat_dict, one_den => get_coeff_mat(one_den,sp_dict) )
    end # for one_den
    #-------
    dentop = make_dentop( den_list, den_coeff_mat_dict, sp_dict )
    push!( dentop_collect, dentop )
    #-------
  end # for one_file

  # make backup for pointer to file
  file_dentop_collect = dentop_collect

  # uniqueness
  unique!( x->reduce(*,x.den_list), dentop_collect )
  @show length(dentop_collect)

  # parent or child
  dentop_collect = get_superior_dentop_collect( dentop_collect ) 
  @show length(dentop_collect)

  #------
  incomplete_dentop_list = dentop_collect
  complete_dentop_list = Vector{DenTop}()
  while !isempty(incomplete_dentop_list)

    new_incomplete_dentop_list, new_complete_dentop_list = 
        update_dentop_list( n_sp, incomplete_dentop_list, den_coeff_mat_dict )

    incomplete_dentop_list = new_incomplete_dentop_list
    union!( complete_dentop_list, new_complete_dentop_list )

    @show length(incomplete_dentop_list) length(complete_dentop_list)
  end # while

###------------------
### CHECK
##for one_dentop in complete_dentop_list
##  println()
##  for one_den in one_dentop.den_list
##    println( one_den )
##  end # for one_den
##end # for one_dentop
###------------------

  complete_dentop_list = map( x->finalize_complete_dentop(x,indep_mom_list,sp_dict), 
                              complete_dentop_list )
  n_dentop = length(complete_dentop_list)
  for index in 1:n_dentop
    complete_dentop = complete_dentop_list[index]
    pos_list = find_all_children( complete_dentop, file_dentop_collect )        

    println()
    println( "--------------" )
    println( "Complete Topology #$(index) covering files: " )
    map( println, file_list[pos_list] )
    println( "--------------" )
    map( println, complete_dentop.den_list )

  end # for index

  return nothing

end # function main


#########
main()
#########

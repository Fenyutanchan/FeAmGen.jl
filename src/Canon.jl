###################################
# const definition
const preferred_vac_den_mom_Dict() = Dict{Int,Vector{Vector{Basic}}}(
  1 => [ [ Basic("q1") ] ],
  2 => [ to_Basic( ["q1", "q2", "q1 + q2"] ) ],
  3 => [ to_Basic( ["q1", "q2", "q3", "q1 + q3", "q2 + q3", "q1 + q2 + q3"] ),
         to_Basic( ["q1", "q2", "q3", "q1 + q2", "q1 + q3", "q2 + q3"] ) ]
) # end preferred_vac_den_mom_Dict
###################################






###################################
# This function is used for the case before canonicalization, 
#   so we only check the coefficient sign of the leading qi.
function normalize_loop_mom_single( 
    loop_mom::Basic
)::Basic
###################################

  leading_qi = (first∘get_loop_momenta)( [ loop_mom ] )
  # qi_list = free_symbols(loop_mom)
  # filter!(sym -> (first∘string)(sym) == 'q', qi_list)
  # sort!(qi_list, by=q->parse(Int, string(q)[2:end]))
  # leading_qi = first(qi_list)

  loop_mom = expand(loop_mom)
  if SymEngine.coeff(loop_mom,leading_qi) < 0
    return expand(-loop_mom)
  end # if

  return loop_mom

end # function normalize_loop_mom_single


###################################
function normalize_loop_mom( 
    loop_den_list::Vector{Basic} 
)::Vector{Basic}
###################################

  @funs Den

  qi_list = get_loop_momenta(loop_den_list)
  # qi_list = free_symbols(loop_den_list)
  # filter!(sym -> (first∘string)(sym) == 'q', qi_list)
  # sort!(qi_list, by=q->parse(Int, string(q)[2:end]))

  new_loop_den_list = zeros( Basic, length(loop_den_list) )
  for index in eachindex(loop_den_list)
    mom, mass, width = get_args(loop_den_list[index])
    mom = expand(mom)
    coeff_list = SymEngine.coeff.( mom, qi_list )
    first_nonzero_index = findfirst( !iszero, coeff_list )
    @assert !isnothing(first_nonzero_index)
    if coeff_list[first_nonzero_index] < 0
      mom = expand(-mom)
    end # if
    new_loop_den_list[index] = Den( mom, mass, width )
  end # for index

  return new_loop_den_list

end # function normalize_loop_mom






#########################################################
# Created by Quan-feng Wu 
# Feb. 16 2023
function gen_loop_mom_canon_map(
    mom_list::Vector{Basic}
)::Dict{Basic,Basic}
#########################################################

  # preferred_vac_den_mom_Dict = Dict{Int, Vector{Vector{Basic}}}(
  #   3 =>  [
  #     map( Basic, ["q1", "q2", "q3", "q1 + q2", "q1 + q3", "q2 + q3"] ),
  #     map( Basic, ["q1", "q2", "q3", "q1 + q3", "q2 + q3", "q1 + q2 + q3"] )
  #   ]
  # )

  q_list = get_loop_momenta(mom_list)
  n_loop = isempty(q_list) ? 0 : (get_loop_index∘last)( q_list )
  @assert q_list == [Basic("q$ii") for ii ∈ 1:n_loop]

  vac_mom_list = mom_list - subs.( mom_list, Ref(Dict(q_list .=> 0)) )
  map!( expand, vac_mom_list, vac_mom_list )
  filter!( !iszero, vac_mom_list )
  map!( normalize_loop_mom_single, vac_mom_list, vac_mom_list )
  unique!(vac_mom_list)

  preferred_flag = n_loop ∈ keys(preferred_vac_den_mom_Dict())

  # if preferred_flag && any(
  #     preferred_mom_list -> vac_mom_list ⊆ preferred_mom_list,
  #     preferred_vac_den_mom_Dict[n_loop]
  #   ) # end any
  #     return Dict{Basic, Basic}()
  # elseif !preferred_flag &&
  #   (isempty∘setdiff)( q_list, vac_mom_list ) && 
  #   all( ≥(0), coefficient_matrix( vac_mom_list, q_list ) )
  #     return Dict{Basic, Basic}()
  # end # if
  
  prefered_all_possible_repl_rules = Dict{Basic, Basic}[]
  all_possible_repl_rules = Dict{Basic,Basic}[]
  
  for selected_mom_indices ∈ permutations(eachindex(vac_mom_list), length(q_list))
    for sign_list ∈ Iterators.product([(1, -1) for _ in q_list]...)
      if last(sign_list) == -1
        break
      end # if

      tmp_coeff_mat = coefficient_matrix( sign_list .* vac_mom_list[selected_mom_indices], q_list )
      if (iszero∘expand∘get_det)(tmp_coeff_mat)
        break
      end # if

      repl_rule = Dict( q_list .=> inv(tmp_coeff_mat) * q_list )
      new_vac_mom_list = (expand∘subs).( vac_mom_list, Ref(repl_rule) )
      @assert new_vac_mom_list[selected_mom_indices] == sign_list .* q_list

      map!( normalize_loop_mom_single, new_vac_mom_list, new_vac_mom_list )
      new_coeff_mat = coefficient_matrix( new_vac_mom_list, q_list )

      if preferred_flag
        if any(
          preferred_mom_list -> new_vac_mom_list ⊆ preferred_mom_list,
          preferred_vac_den_mom_Dict()[n_loop]
        ) # end any
          # sort_index = gen_repl_rule_sort_index( subs.(mom_list, Ref(repl_rule)) )
          push!( prefered_all_possible_repl_rules, repl_rule )
        elseif all( ≥(0), new_coeff_mat )
          # sort_index = gen_repl_rule_sort_index( subs.(mom_list, Ref(repl_rule)) )
          push!( all_possible_repl_rules, repl_rule )
        end # if
      elseif all( ≥(0), new_coeff_mat )
        # sort_index = gen_repl_rule_sort_index( subs.(mom_list, Ref(repl_rule)) )
        push!( all_possible_repl_rules, repl_rule )
      end # if
    end # for sign_list
  end # for selected_mom_indices

  if preferred_flag
    if isempty( prefered_all_possible_repl_rules )
      printstyled("Warning: There is no fetching preferred_flag loop momenta list.\n"; color=:yellow)

      sort!( all_possible_repl_rules; 
             by=repl->get_repl_rule_sort_index( subs.( mom_list, Ref(repl) ) ) )
      return first( all_possible_repl_rules )
    end

    sort!( prefered_all_possible_repl_rules; 
           by=repl->get_repl_rule_sort_index( subs.( mom_list, Ref(repl) ) ) )
    return first( prefered_all_possible_repl_rules )
  end # if preferred_flag

  sort!( all_possible_repl_rules;
         by=repl->get_repl_rule_sort_index( subs.( mom_list, Ref(repl) ) ) )
  return first( all_possible_repl_rules )

end # function gen_loop_mom_canon_map






#########################################################
# Created by Quanfeng Wu 
# Mar. 26 2023
# 
# A simple trial for sorting the replace rules generated in `gen_loop_mom_canon_map`.
function get_repl_rule_sort_index( 
    mom_list::Vector{Basic} 
)::Basic
#########################################################

  tmp_mom_list = unique( abs, mom_list )
  q_list = get_loop_momenta( tmp_mom_list )
  q_index_list = map( get_loop_index, q_list )
  k_list = get_ext_momenta( tmp_mom_list )
  k_index_list = map( get_ext_index, k_list )

  q_coeff_mat = coefficient_matrix( tmp_mom_list, q_list )
  k_coeff_mat = coefficient_matrix( tmp_mom_list, k_list )

  map!( abs, q_coeff_mat, q_coeff_mat )
  map!( abs, k_coeff_mat, k_coeff_mat )

  result = transpose( q_coeff_mat * q_index_list ) * ( k_coeff_mat * k_index_list )

  return result
end # function get_repl_rule_sort_index






#####################################################
function canonicalize_amp( 
    loop_den_list::Vector{Basic},  
    amp_lorentz_list::Vector{Basic} 
)::Tuple{Vector{Basic},Vector{Basic}}
######################################################

  # n_loop = get_n_loop( loop_den_list )
  q_list = get_loop_momenta( loop_den_list )
  n_loop = isempty(q_list) ? 0 : (get_loop_index∘last)( q_list )

  if n_loop == 0
    return loop_den_list, amp_lorentz_list 
  end # if

  mom_list = map( first∘get_args, loop_den_list )
  canon_map = gen_loop_mom_canon_map(mom_list) 

  if isempty(canon_map)
    new_loop_den_list = loop_den_list 
    new_amp_lorentz_list = amp_lorentz_list 
  else
    new_loop_den_list = map( den->subs(den,canon_map), loop_den_list )
    new_amp_lorentz_list = map( amp->subs(amp,canon_map), amp_lorentz_list )
  end # if

  # Normalize the leading coefficient sign of the loop momenta.
  new_loop_den_list = normalize_loop_mom( new_loop_den_list )

  # CHECK begin
  # qi_list = Basic[ Basic("q$ii") for ii in 1:n_loop ]
  new_mom_list = map( first∘get_args, new_loop_den_list )
  coeff_list = union( coefficient_matrix( new_mom_list, q_list ) )
  @assert all( ≥(0), coeff_list ) "$new_mom_list"
  # CHECK end

  return new_loop_den_list, new_amp_lorentz_list 

end # function canonicalize_amp



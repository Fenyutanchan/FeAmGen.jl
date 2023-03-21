
###################################
# This function is used for the case before canonicalization, 
#   so we only check the coefficient sign of the leading qi.
function normalize_loop_mom_single( 
    loop_mom::Basic
)::Basic
###################################

  qi_list = free_symbols(loop_mom)
  filter!(sym -> (first∘string)(sym) == 'q', qi_list)
  sort!(qi_list, by=q->parse(Int, string(q)[2:end]))
  leading_qi = first(qi_list)

  loop_mom = expand(loop_mom)
  if SymEngine.coeff(loop_mom,leading_qi) < 0
    return expand(-loop_mom)
  end # if

  return loop_mom

end # function normalize_loop_mom


###################################
function normalize_loop_mom( 
    loop_den_list::Vector{Basic} 
)::Vector{Basic}
###################################

  @funs Den

  qi_list = free_symbols(loop_den_list)
  filter!(sym -> (first∘string)(sym) == 'q', qi_list)
  sort!(qi_list, by=q->parse(Int, string(q)[2:end]))

  new_loop_den_list = zeros( Basic, length(loop_den_list) )
  for index in 1:length(loop_den_list)
    mom, mass, width = get_args(loop_den_list[index])
    mom = expand(mom)
    if any( <(0), SymEngine.coeff.(mom,qi_list) )
      mom = -mom
    end # if
    new_loop_den_list[index] = Den(expand(mom),mass,width)
  end # for index

  return new_loop_den_list

end # function normalize_loop_mom






#########################################################
# Created by Quanfeng Wu 
# Feb. 16 2023
function gen_loop_mom_canon_map(
    mom_list::Vector{Basic}
)::Dict{Basic,Basic}
#########################################################

  function is_sym_index_format(input::Basic, sym::Union{Basic, String})::Bool
    if SymEngine.get_symengine_class(input) != :Symbol
      return false
    end
    if SymEngine.get_symengine_class(Basic(sym)) != :Symbol
      return false
    end

    input_str = string(input)
    sym_str = string(sym)
    sym_len = length(sym_str)
    
    if !startswith(input_str, sym_str)
      return false
    end

    index_part_range  = findnext(r"[1-9]\d*", input_str, sym_len+1)
    if isnothing(index_part_range)
      return false
    end
    return sym_len + length(index_part_range) == length(input_str)
  end
  is_loop_mom(mom::Basic)::Bool = is_sym_index_format(mom, "q")
  is_ext_mom(mom::Basic)::Bool = is_sym_index_format(mom, "K")

  function get_sym_index(mom::Basic, sym::Union{Basic, String})::Int
    @assert is_sym_index_format(mom, sym)
    mom_str = string(mom)
    sym_len = (length∘string)(sym)
    return parse(Int, mom_str[sym_len+1:end])
  end
  get_loop_index(mom::Basic)::Int = get_sym_index(mom, "q")
  get_ext_index(mom::Basic)::Int = get_sym_index(mom, "K")

  function coefficient_matrix(mom_poly_list::Vector{Basic}, mom_list::Vector{Basic})::Matrix{Basic}
    @assert all(sym -> SymEngine.get_symengine_class(sym) == :Symbol, mom_list)

    return reduce(
      vcat,
      (transpose ∘ map)( q -> SymEngine.coeff(mom, q), mom_list )
        for mom ∈ mom_poly_list
    )
  end

  function get_loop_momenta(mom_list::Vector{Basic})::Vector{Basic}
    single_mom_list = free_symbols(mom_list)
    q_list = filter(is_loop_mom, single_mom_list)
    sort!(q_list; by=get_loop_index)
    return q_list
  end

  function get_ext_momenta(mom_list::Vector{Basic})::Vector{Basic}
    single_mom_list = free_symbols(mom_list)
    k_list = filter(is_ext_mom, single_mom_list)
    sort!(k_list; by=get_ext_index)
    return k_list
  end

  function gen_repl_rule_sort_index(mom_list::Vector{Basic})::Vector{Basic}
    local tmp_mom_list = unique( abs, mom_list )
    q_list = get_loop_momenta( tmp_mom_list )
    @assert q_list == [Basic("q$ii") for ii ∈ eachindex(q_list)]

    count_mat = zeros(Basic, length(tmp_mom_list), length(q_list))
    for (mom_index, tmp_mom) ∈ enumerate(tmp_mom_list)
      q_coeff_list = SymEngine.coeff.(tmp_mom, q_list)
      k_list = get_ext_momenta([tmp_mom])
      k_index_sum = if isempty(k_list)
        zero(Basic)
      else
        sum( get_ext_index, k_list )
      end
      count_mat[mom_index, findall(!iszero, q_coeff_list)] .+= k_index_sum
    end

    result_list = Basic[]
    for col ∈ eachcol(count_mat)
      tmp_col = filter( !iszero, col )
      if isempty(tmp_col)
        push!( result_list, (Basic∘typemax)(Int) )
      else
        push!( result_list, minimum(tmp_col) )
      end
    end
    return result_list
  end

  target_loop_mom_list = Dict{Int, Vector{Vector{Basic}}}(
    3 =>  [
      map( Basic, ["q1", "q2", "q3", "q1 + q2", "q1 + q3", "q2 + q3"] ),
      map( Basic, ["q1", "q2", "q3", "q1 + q3", "q2 + q3", "q1 + q2 + q3"] )
    ]
  )

  q_list = get_loop_momenta(mom_list)
  n_loop = length(q_list)
  @assert q_list == [Basic("q$ii") for ii ∈ eachindex(q_list)]

  # orig_mom_list = copy(mom_list)
  # sort!(
  #   orig_mom_list;
  #   by=mom->
  #     (
  #       try  
  #         sum(get_ext_index, filter(is_ext_mom, free_symbols(mom)))
  #       catch
  #         @assert (isempty∘filter)(is_ext_mom, free_symbols(mom))
  #         typemax(Int)
  #       end,
  #       findfirst(!iszero, SymEngine.coeff.(mom, q_list)),
  #       (length∘free_symbols)(mom)
  #     ) # end by
  # ) # end sort!
  tmp_mom_list = mom_list - subs.( mom_list, Ref(Dict(q => 0 for q ∈ q_list)) )
  map!( expand, tmp_mom_list, tmp_mom_list )
  filter!( !iszero, tmp_mom_list )
  map!(normalize_loop_mom_single, tmp_mom_list, tmp_mom_list)
  unique!(tmp_mom_list)

  target_loop_flag = n_loop ∈ keys(target_loop_mom_list)

  if target_loop_flag && any(
      target_mom_list -> tmp_mom_list ⊆ target_mom_list,
      target_loop_mom_list[n_loop]
    ) # end any
      return Dict{Basic, Basic}()
  elseif !target_loop_flag &&
    (isempty∘setdiff)( q_list, tmp_mom_list ) && 
    all( ≥(0), coefficient_matrix( tmp_mom_list, q_list ) )
      return Dict{Basic, Basic}()
  end # if
  
  target_all_repl_rules = Tuple{Vector{Basic}, Dict{Basic, Basic}}[]
  all_repl_rules = Tuple{Vector{Basic}, Dict{Basic, Basic}}[]
  
  for selected_mom_indices ∈ permutations(eachindex(tmp_mom_list), length(q_list))
    for sign_list ∈ Iterators.product([(1, -1) for _ in q_list]...)
      if last(sign_list) == -1
        break
      end # if

      tmp_coeff_mat = coefficient_matrix( sign_list .* tmp_mom_list[selected_mom_indices], q_list )
      if (iszero ∘ expand ∘ get_det)(tmp_coeff_mat)
        break
      end # if

      repl_rule = Dict(q_list .=> inv(tmp_coeff_mat) * q_list)
      new_mom_list = (expand ∘ subs).(tmp_mom_list, Ref(repl_rule))
      @assert new_mom_list[selected_mom_indices] == sign_list .* q_list

      map!( normalize_loop_mom_single, new_mom_list, new_mom_list )
      new_coeff_mat = coefficient_matrix( new_mom_list, q_list )

      if target_loop_flag
        if any(
          target_mom_list -> new_mom_list ⊆ target_mom_list,
          target_loop_mom_list[n_loop]
        ) # end any
          sort_index = gen_repl_rule_sort_index( subs.(mom_list, Ref(repl_rule)) )
          push!( target_all_repl_rules, (sort_index, repl_rule) )
        elseif all(≥(0), new_coeff_mat)
          sort_index = gen_repl_rule_sort_index( subs.(mom_list, Ref(repl_rule)) )
          push!( all_repl_rules, (sort_index, repl_rule) )
        end
      elseif all(≥(0), new_coeff_mat)
        sort_index = gen_repl_rule_sort_index( subs.(mom_list, Ref(repl_rule)) )
        push!( all_repl_rules, (sort_index, repl_rule) )
      end # if
    end # for sign_list
  end # for selected_mom_indices

  if target_loop_flag
    if isempty(target_all_repl_rules)
      printstyled("Warning: There is no fetching target loop momenta list.\n"; color=:yellow)
      sort!( all_repl_rules; by=first )

tmp_new_mom_list = map( mom -> subs(mom, (last∘first)(all_repl_rules) ), mom_list )
@show gen_repl_rule_sort_index(mom_list)
@show gen_repl_rule_sort_index(tmp_new_mom_list)

      return (last∘first)(all_repl_rules)
    end

    sort!( target_all_repl_rules; by=first )

tmp_new_mom_list = map( mom -> subs(mom, (last∘first)(target_all_repl_rules) ), mom_list )
@show gen_repl_rule_sort_index(mom_list)
@show gen_repl_rule_sort_index(tmp_new_mom_list)

    return (last∘first)(target_all_repl_rules)
  else
    sort!( all_repl_rules; by=first )

tmp_new_mom_list = map( mom -> subs(mom, (last∘first)(all_repl_rules) ), mom_list )
@show gen_repl_rule_sort_index(mom_list)
@show gen_repl_rule_sort_index(tmp_new_mom_list)

    return (last∘first)(all_repl_rules)
  end

end # function gen_loop_mom_canon_map






#####################################################
function canonicalize_amp( 
    loop_den_list::Vector{Basic},  
    amp_lorentz_list::Vector{Basic} 
)::Tuple{Vector{Basic},Vector{Basic}}
######################################################

  n_loop = get_n_loop( loop_den_list )

  if n_loop == 0
    return loop_den_list, amp_lorentz_list 
  end # if

  mom_list = (first∘get_args).(loop_den_list)
  canon_map = gen_loop_mom_canon_map(mom_list) 

  if isempty(canon_map)
    new_loop_den_list = loop_den_list 
    new_amp_lorentz_list = amp_lorentz_list 
  else
    new_loop_den_list = map( den -> subs( den, canon_map ), loop_den_list )
    new_amp_lorentz_list = map( amp -> subs( amp, canon_map ), amp_lorentz_list )
  end # if

  # Normalize the leading coefficient sign of the loop momenta.
  new_loop_den_list = normalize_loop_mom( new_loop_den_list )

  # CHECK begin
  qi_list = Basic[ Basic("q$ii") for ii in 1:n_loop ]
  den_mom_list = map( x -> (first∘get_args)(x), new_loop_den_list )
  coeff_list = union( map( mom -> SymEngine.coeff.( mom, qi_list ), den_mom_list )... )
  @assert all( x->x>=0, coeff_list )
  # CHECK end

  return new_loop_den_list, new_amp_lorentz_list 

end # function canonicalize_amp



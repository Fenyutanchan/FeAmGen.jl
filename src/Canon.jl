
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

  q_list = free_symbols(mom_list)
  filter!(q -> (first∘string)(q) == 'q', q_list)
  sort!(q_list, by=q->parse(Int, string(q)[2:end]))

  tmp_mom_list = mom_list - subs.(mom_list, Ref(Dict(q_ => 0 for q_ ∈ q_list)))
  filter!(!iszero, tmp_mom_list)
  unique!(mom -> abs(mom), tmp_mom_list)
  sort!( tmp_mom_list, by=mom->(findfirst(!iszero, SymEngine.coeff.(mom,q_list)), (length∘free_symbols)(mom)) )

  if (isempty∘setdiff)(q_list,tmp_mom_list)
    check_flag = all(
        mom -> begin
            coeff_list = SymEngine.coeff.(mom, q_list)
            all( x -> x ≥ 0, coeff_list ) || all( x -> x ≤ 0, coeff_list )
        end,
        tmp_mom_list
    ) # end all
    if check_flag
      return Dict{Basic,Basic}()
    end # if
  end # if

  for selected_mom_indices in combinations(eachindex(tmp_mom_list), length(q_list))
    counter = 0
    counter_target = 2^(length(q_list)-1)
    for sign_list in Base.product([(1, -1) for _ in q_list]...)
      counter += 1
      if counter > counter_target
        break
      end

      coeff_matrix = reduce(
          vcat,
          transpose(SymEngine.coeff.(mom_, q_list))
              for mom_ in (sign_list .* tmp_mom_list[selected_mom_indices])
      ) # end reduce
@show coeff_matrix typeof(coeff_matrix)
      if (iszero∘expand∘get_det)(coeff_matrix)
        break
      end # if
      replacement_rules = Dict(q_list .=> inv(coeff_matrix) * q_list)
      new_mom_list = (expand∘subs).(tmp_mom_list, Ref(replacement_rules))
      @assert new_mom_list[selected_mom_indices] == sign_list .* q_list

      check_flag = all(
          mom_ -> begin
              coeff_list = SymEngine.coeff.(mom_, q_list)
              all( x -> x ≥ 0, coeff_list) || all( x -> x ≤ 0, coeff_list)
          end,
          new_mom_list
      ) # end all
      if check_flag
        return replacement_rules
      end # if

    end # for sign_list
  end #for selected_mom_indices

  return Dict{Basic,Basic}()

end # function gen_loop_mom_canon_map






#####################################################
function canonicalize_amp( 
    n_loop::Int64, 
    loop_den_list::Vector{Basic},  
    amp_lorentz_list::Vector{Basic} 
)::Tuple{Vector{Basic},Vector{Basic}}
######################################################

  if n_loop == 0
    return loop_den_list, amp_lorentz_list 
  end # if

  mom_list = (first∘get_args).(loop_den_list)
  canon_map = gen_loop_mom_canon_map(mom_list) 

  println( "Canonicalization map:" )
  for x in canon_map
    println( "  ", x )
  end # for x
  println()

  if isempty(canon_map)
    new_loop_den_list = loop_den_list 
    new_amp_lorentz_list = amp_lorentz_list 
  else
    new_loop_den_list = map( den -> subs( den, canon_map... ), loop_den_list )
    new_amp_lorentz_list = map( amp -> subs( amp, canon_map... ), amp_lorentz_list )
  end # if

  # Normalize the leading coefficient sign of the loop momenta.
  new_loop_den_list = normalize_loop_mom( new_loop_den_list )

  println( "Old loop_den_list:" )
  for den in loop_den_list
    println( "  ", den )
  end # for den
  println()

  println( "New loop_den_list:" )
  for den in new_loop_den_list
    println( "  ", den )
  end # for den
  println()

  println( "Old numerator list:" )
  for one_lorentz in amp_lorentz_list
    println( "  ", one_lorentz )
  end # for one_lorentz
  println()

  println( "New numerator list:" )
  for one_lorentz in new_amp_lorentz_list
    println( "  ", one_lorentz )
  end # for one_lorentz
  println()

  # CHECK begin
  qi_list = Basic[ Basic("q$ii") for ii in 1:n_loop ]
  den_mom_list = map( x -> (first∘get_args)(x), new_loop_den_list )
  coeff_list = union( map( mom -> SymEngine.coeff.( mom, qi_list ), den_mom_list )... )
  @assert all( x->x>=0, coeff_list )
  # CHECK end

  return new_loop_den_list, new_amp_lorentz_list 

end # function canonicalize_amp



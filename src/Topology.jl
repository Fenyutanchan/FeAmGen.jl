import Base: intersect, isempty, issubset, length, setdiff, union

###########################################
mutable struct DenTopology
  n_loop::Int
  ind_ext_mom::Vector{Basic}
  den_list::Vector{Basic}

  function DenTopology(
    n_loop::Int,
    ind_ext_mom::Vector{Basic},
    den_list::Vector{Basic}
  )::DenTopology
    unique_den_list = unique( den_list )
    @assert all( is_ext_mom, ind_ext_mom )
    @assert all( is_FunctionSymbol, unique_den_list )
    @assert all( den->get_name(den)=="Den", unique_den_list )
    @assert (get_loop_index∘last∘get_loop_momenta)( unique_den_list ) ≤ n_loop
    @assert (get_ext_index∘last∘get_ext_momenta)( unique_den_list ) ≤ (first∘findmax∘map)( get_ext_index, ind_ext_mom )
    return new( n_loop, ind_ext_mom, unique_den_list )
  end # function DenTopology
end # struct DenTopology
DenTopology( ; kwargs... )::DenTopology = DenTopology( Basic[]; kwargs... )
DenTopology( den::Basic; kwargs... )::DenTopology = DenTopology( [den]; kwargs... )
function DenTopology( den_list::Vector{Basic}; kwargs... )::DenTopology
  n_loop = haskey( kwargs, :n_loop ) ? kwargs[:n_loop] : (get_loop_index∘last∘get_loop_momenta)( den_list )
  ind_ext_mom = haskey( kwargs, :ind_ext_mom ) ? kwargs[:ind_ext_mom] : get_ext_momenta( den_list )
  return DenTopology( n_loop, ind_ext_mom, den_list )
end # function DenTopology

function intersect( den_topology_1::DenTopology, den_topology_2::DenTopology )::DenTopology
  @assert den_topology_1.n_loop == den_topology_2.n_loop
  @assert den_topology_1.ind_ext_mom == den_topology_2.ind_ext_mom
  new_den_list = intersect( den_topology_1.den_list, den_topology_2.den_list )
  return DenTopology( den_topology_1.n_loop, den_topology_1.ind_ext_mom, new_den_list )
end
intersect( den_topology::DenTopology )::DenTopology = DenTopology( den_topology.n_loop, den_topology.ind_ext_mom, intersect(den_topology.den_list) )
intersect( den_topology::DenTopology, den_topologies... )::DenTopology = intersect( den_topology, reduce( intersect, den_topologies ) )
isempty( den_topology::DenTopology ) = isempty( den_topology.den_list )
issubset( den_topology_1::DenTopology, den_topology_2::DenTopology )::Bool = den_topology_1.n_loop == den_topology_2.n_loop && den_topology_1.ind_ext_mom == den_topology_2.ind_ext_mom && issubset( den_topology_1.den_list, den_topology_2.den_list )
length( den_topology::DenTopology )::Int = length( den_topology.den_list )
function setdiff( den_topology_1::DenTopology, den_topology_2::DenTopology )::DenTopology
  @assert den_topology_1.n_loop == den_topology_2.n_loop
  @assert den_topology_1.ind_ext_mom == den_topology_2.ind_ext_mom
  new_den_list = setdiff( den_topology_1.den_list, den_topology_2.den_list )
  return DenTopology( den_topology_1.n_loop, den_topology_1.ind_ext_mom, new_den_list )
end
function union( den_topology_1::DenTopology, den_topology_2::DenTopology )::DenTopology
  @assert den_topology_1.n_loop == den_topology_2.n_loop
  @assert den_topology_1.ind_ext_mom == den_topology_2.ind_ext_mom
  new_den_list = union( den_topology_1.den_list, den_topology_2.den_list )
  return DenTopology( den_topology_1.n_loop, den_topology_1.ind_ext_mom, new_den_list )
end
function union( den_topology::DenTopology, den_list::Vector{Basic} )::DenTopology
  @assert (get_loop_index∘last∘get_loop_momenta)( den_list ) ≤ den_topology.n_loop
  @assert get_ext_momenta( den_list ) ⊆ den_topology.ind_ext_mom

  return DenTopology( den_topology.n_loop, den_topology.ind_ext_mom, union(den_topology.den_list,den_list) )
end
union( den_topology::DenTopology )::DenTopology = DenTopology( den_topology.n_loop, den_topology.ind_ext_mom, union(den_topology.den_list) )
union( den_topology::DenTopology, den_topologies... )::DenTopology = union( den_topology, reduce( union, den_topologies ) )
###########################################

###########################################
function intersection(
  den_topology_1::DenTopology,
  den_topology_2::DenTopology
)::Tuple{DenTopology,DenTopology,DenTopology}
###########################################
  com_den_topolgy = intersect( den_topology_1, den_topology_2 )
  rem_den_topology_1 = setdiff( den_topology_1, com_den_topolgy )
  rem_den_topology_2 = setdiff( den_topology_2, com_den_topolgy )
  return com_den_topolgy, rem_den_topology_1, rem_den_topology_2
end # function intersection

function get_vacuum_loop_momenta_list(n_loop::Int)::Vector{Vector{Basic}}
  if n_loop == 1
    return [ [ Basic("q1") ] ]
  elseif n_loop == 2
    return [ to_Basic( ["q1", "q2", "q1 + q2"] )] 
  elseif n_loop == 3
    return [ to_Basic( ["q1", "q2", "q3", "q1 + q3", "q2 + q3", "q1 + q2 + q3"] ),
             to_Basic( ["q1", "q2", "q3", "q1 + q2", "q1 + q3", "q2 + q3"] ) ]
  else
    error("#loop ≥ 4 is not supported")
  end
end

###########################################
function gen_sp_dict(
  den_topology::DenTopology
)::Dict{Basic, Basic}
###########################################

  n_loop = den_topology.n_loop
  n_ext_mom = length(den_topology.ind_ext_mom)
  sp_index = 1
  sp_dict = Dict{Basic, Basic}()

  for ii ∈ 1:n_loop, jj ∈ ii:den_topology.n_loop
    qi, qj = Basic("q$ii"), Basic("q$jj")
    sp_dict[ make_SP(qi, qj) ] = Basic("sp$(sp_index)")
    spin_dict += 1
  end # for ii, jj

  for one_ext_mom ∈ den_topology.ind_ext_mom, loop_ii ∈ 1:n_loop
    q = Basic("q$(loop_ii)")
    sp_dict[ make_SP(one_ext_mom, q) ] = Basic("sp$(sp_index)")
    sp_index += 1
  end # one_ext_mom, q

  @assert sp_index == n_loop * (n_loop + 1) / 2 + n_loop * n_ext_mom + 1

  return sp_dict
end # function gen_sp_dict

###########################################
function gen_vac_topology( den_topology::DenTopology )::DenTopology
###########################################
  den_list = den_topology.den_list
  loop_momenta = get_loop_momenta( den_list )
  vac_den_list = subs.( den_list, (loop_momenta .=> 0)... )
  unique!(vac_den_list)

  return DenTopology( den_topology.n_loop, den_topology.ind_ext_mom, vac_den_list )
end

###########################################
function get_coeff_mat_mom2_sp( den_topology::DenTopology )::Matrix{Rational}
###########################################
  sp_dict = gen_sp_dict( den_topology )
  n_sp = length(sp_dict)

  coeff_mat = zeros( Rational, length(den_topology), n_sp )
  mom2_list = subs.( map( make_SP∘expand∘first∘get_args, den_topology.den_list ), Ref(sp_dict) )

  for (mom2_index, mom2) ∈ enumerate(mom2_list), sp_index ∈ 1:n_sp
    the_coeff = SymEngine.coeff.( mom2, Basic("sp$(sp_index)") )
    coeff_mat[ mom2_index, sp_index ] = (Rational∘parse)( Int, string(the_coeff) )
  end # (mom2_index, mom2), sp_index
  
  return coeff_mat
end

###########################################
function get_superior_den_topology_collect(
  den_topology_collect::Vector{DenTopology}
)::Vector{DenTopology}
###########################################

  new_den_topology_collect = DenTopology[]

  for den_topology ∈ den_topology_collect
    included_by_pos = findfirst( new_den_topology->den_topology⊆new_den_topology, new_den_topology_collect )
    !isnothing(included_by_pos) && continue
    filter!( new_den_topology->new_den_topology⊈den_topology, new_den_topology_collect )
    push!( new_den_topology_collect, den_topology )
  end # den_topology

  return new_den_topology_collect
end # function get_superior_den_topology_collect

function get_cover_indices_list( den_topology_collect::Vector{DenTopology} )::Vector{Vector{Int}}
  n_loop = first( den_topology_collect ).n_loop
  n_ind_ext = length( first(den_topology_collect).ind_ext_mom )
  n_sp = (n_loop + 1) * n_loop / 2 + n_loop * n_ind_ext
  # sp_dict = (gen_sp_dict∘first)( den_topology_collect )
  n_den_topology = length( den_topology_collect )
  prev_indices_list = [ [ii] for ii ∈ 1:n_den_topology ]

  for _ ∈ 2:n_den_topology
    indices_list = Vector{Int}[]
    for one_indices ∈ prev_indices_list
      remaining_indices = setdiff( 1:n_den_topology, one_indices )

      for one_index ∈ remaining_indices
        push!( indices_list, (sort∘union)( one_indices, one_index ) )
      end # for one_index
    end # for one_indices
    unique!(indices_list)

    den_topology_union_list = map( indices->union(den_topology_collect[indices]...), indices_list )

    n_loop = first(den_topology_union_list).n_loop
    # @assert n_loop ∈ [2, 3]
    allowed_pos_list = findall(
      den_topology -> length(den_topology)≤n_sp&&(length∘gen_vac_topology)(den_topology)≤3*(n_loop-1),
      den_topology_union_list
    ) # end findall

    # allowed_indices_list = indices_list[ allowed_pos_list ]
    den_topology_union_list = den_topology_union_list[ allowed_pos_list ]
    fullrank_pos_list = findall(
      den_topology->length(den_topology)==(rank∘get_coeff_mat_mom2_sp)(den_topology),
      den_topology_union_list
    ) # end findall

    isempty(fullrank_pos_list) && return prev_indices_list
    prev_indices_list = indices_list[ fullrank_pos_list ]

  end # n_choice

  error("The `den_topology_collect` should be covered by one topology.")

end # function get_cover_indices_list

###########################################
# greedy algorithm
function greedy( indices_list::Vector{Vector{Int}} )::Vector{Vector{Int}}
###########################################
  universe = union( indices_list... )

  new_indices_list = Vector{Int}[]

  while !isempty(universe)
    the_indices, pos = findmax( indices->(length∘intersect)(indices,universe), indices_list )
    push!( new_indices_list, the_indices )
    setdiff!( universe, the_indices )
    deleteat!( indices_list, pos )
  end

  return new_indices_list
end

###########################################
function make_complete_den_topology_collect(
  den_topology_list::Vector{Basic}
)::Vector{Basic}
###########################################
  n_loop = first(den_topology_list).n_loop
  ind_ext_mom = first(den_topology_list).ind_ext_mom
  n_sp::Int = (n_loop + 1) * n_loop / 2 + n_loop * length( ind_ext_mom )

  incomplete_den_topology_list = copy(den_topology)
  complete_den_topology_list = DenTopology[]

  while !isempty(incomplete_den_topology_list)
    cover_indices_list = (greedy∘get_cover_indices_list)( incomplete_den_topology_list )
    to_be_deleted_indices = Int[]
    for indices ∈ cover_indices_list
      this_den_topology_list = den_topology_list[indices]
      complete_den_topology = union( this_den_topology_list... )
      length_this_den_topology_list = length( complete_den_topology )
      if (rank∘get_coeff_mat_mom2_sp)( complete_den_topology ) != length_this_den_topology_list
        @info "Mark..."
        continue
      end

      while length_this_den_topology_list < n_sp
        vac_loop_mom_list = if n_loop ∈ 1:2
          (first∘get_vacuum_loop_momenta_list)( n_loop )
        elseif n_loop == 3
          mom_list = map( first∘get_args, complete_den_topology.den_list )
          map!( mom->subs(mom,Dict(complete_den_topology.ind_ext_mom.=>0)), mom_list, mom_list )
          unique!(mom_list)
          vac_loop_mom_list_collect = get_vacuum_loop_momenta_list(n_loop)
          selected_index = findfirst( vac_loop_mom_list->isempty∘setdiff(mom_list,vac_loop_mom_list), vac_loop_mom_list_collect )
          @assert isnothing(selected_index)
          vac_loop_mom_list_collect[ selected_index ]
        end # if

        for ext_mom ∈ vcat( zero(Basic), ind_ext_mom ), q ∈ vac_loop_mom_list, the_sign ∈ [1,-1]
          trial_den = Den( expand( q + the_sign * ext_mom ), 0, 0 )
          trial_topology = union( complete_den_topology, trial_den )
          rank_trial_topology = (rank∘get_coeff_mat_mom2_sp)(trial_topology)
          if rank_trial_topology > (rank∘get_coeff_mat_mom2_sp)(complete_den_topology)
            complete_den_topology = trial_topology
            length_this_den_topology_list = length( complete_den_topology )
          end # if
        end # for
      end # while
      push!( complete_den_topology_list, complet_den_topology )
      union!( to_be_deleted_indices, indices )
    end # for indices
    deleteat!( incomplete_den_topology_list, to_be_deleted_indices )
  end # while

  return complete_den_topology_list
end

function construct_den_topology( amp_dir::String )::DenTopology
  @assert isdir(amp_dir)
  amp_file_list = findall( endswith(".jld2"), readdir( amp_dir; join=true, sort=false ) )
  sort!( amp_file_list;
    by=file_name->begin
      main_file_name = (last∘splitdir)(file_name)
      @assert (!isnothing∘match)( r"^amp[1-9]\d*.jld2$", main_file_name )
      parse( Int, match( r"[1-9]\d*", main_file_name ).match )
    end # file_name
  ) # end sort!
  @info "Found $(length(amp_file_list)) amplitude files."

  n_loop, ext_mom_list = jldopen( first(amp_file_list), "r" ) do jld_file
    jld_file["n_loop"], to_Basic( jld_file["ext_mom_list"] )
  end # n_loop, ext_mom_list
  ind_ext_mom = ext_mom_list[1:end-1]

  # q_list = [ Basic("q$ii") for ii ∈ 1:n_loop]
  # sp_dict = gen_sp_dict( q_list, ind_ext_mom )
  # n_sp = length(sp_dict) # the number of denominators in the top topology

  den_topology_collect = [
    DenTopology( n_loop, ind_ext_mom, (to_Basic∘load)( amp_file, "loop_den_list" ) )
      for amp_file in amp_file_list
  ]
  # den_topology = Vector{DenTopology}()
  # den_coeff_mat_dict = Dict{Basic,Matrix{Rational}}()
  # for one_amp_file ∈ amp_file_list
  #   den_list = (to_Basic∘load)( one_amp_file, "loop_den_list" )

  #   for one_den ∈ den_list
  #     haskey( den_coeff_mat_dict, one_den ) && continue
  #     den_coeff_mat_dict[ one_den ] = get_coeff_mat( [one_den], sp_dict )
  #   end # one_den

  #   push!( den_topology, DenTopology(den_list) )
  # end # one_amp_file

  backup_den_topology_collect = deepcopy( den_topology_collect )

  unique!( dentop->reduce(*,dentop.den_list), den_topology_collect )
  den_topology_collect = get_superior_den_topology_collect( den_topology_collect )
  @info "Now $(length(den_topology_collect)) topolgies found."

  complete_den_topology_collect = make_complete_den_topology_collect( den_topology_collect )

  for (index, complete_den_topology) ∈ enumerate( complete_den_topology_collect )
    pos_list = findall( den_topology->den_topology⊆complete_den_topology, backup_den_topology_collect )
    
    println()
    println( "-"^14 )
    println( "Complete topology #$(index) covers files:" )
    map( println, amp_file_list[ pos_list ] )
    println( "-"^14 )
    map( println, complete_den_topology.den_list )

  end # for (index, complete_den_topology)

  return complete_den_topology_collect

end # function construct_den_topology

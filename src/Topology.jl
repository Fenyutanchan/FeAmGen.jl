import Base: intersect, isempty, issubset, length, setdiff, union

###########################################
mutable struct DenTopology
  n_loop::Int
  ind_ext_mom::Vector{Basic}
  den_list::Vector{Basic}
end # mutable struct DenTopology
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
union( den_topology::DenTopology, den::Basic ) = union( den_topology, [den] )
union( den_topology::DenTopology )::DenTopology = DenTopology( den_topology.n_loop, den_topology.ind_ext_mom, union(den_topology.den_list) )
union( den_topology::DenTopology, den_topologies... )::DenTopology = union( den_topology, reduce( union, den_topologies ) )
###########################################

###########################################
function is_valid_dentopology( den_topology::DenTopology )::Bool
###########################################
  n_loop = den_topology.n_loop
  ind_ext_mom = den_topology.ind_ext_mom
  den_list = den_topology.den_list

  !all( is_ext_mom, ind_ext_mom ) && return false
  if !isempty(den_list)
    unique_den_list = unique( den_list )
    !all( is_FunctionSymbol, unique_den_list ) && return false
    !all( den->get_name(den)=="Den", unique_den_list ) && return false

    loop_momenta = get_loop_momenta( unique_den_list )
    ext_momenta = get_ext_momenta( unique_den_list )

    max_loop_index = isempty( loop_momenta ) ? 0 : (get_loop_index∘last∘get_loop_momenta)( unique_den_list )
    max_ext_index = isempty( ext_momenta ) ? 0 : (get_ext_index∘last∘get_ext_momenta)( unique_den_list )
    max_ind_ext_index = isempty( ind_ext_mom ) ? 0 : (first∘findmax∘map)( get_ext_index, ind_ext_mom )

    max_loop_index > n_loop && return false
    max_ext_index > max_ind_ext_index && return false
  end # if

  return true
end # function is_valid_dentopology

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
    sp_index += 1
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
  ext_momenta = get_ext_momenta( den_list )
  vac_den_list = subs.( den_list, (ext_momenta .=> 0)... )
  unique!(vac_den_list)

  return DenTopology( den_topology.n_loop, den_topology.ind_ext_mom, vac_den_list )
end # function gen_vac_topology

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
  n_sp::Int = (n_loop + 1) * n_loop / 2 + n_loop * n_ind_ext

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
    unique!( sort, indices_list )

    den_topology_union_list = map( indices->union(den_topology_collect[indices]...), indices_list )

    allowed_pos_list = findall(
      den_topology -> length(den_topology)≤n_sp&&(length∘gen_vac_topology)(den_topology)≤3*(n_loop-1),
      den_topology_union_list
    ) # end findall

    indices_list = indices_list[ allowed_pos_list ]
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
  @show indices_list

  universe = union( indices_list... )

  new_indices_list = Vector{Int}[]

  while !isempty(universe)
    _, pos = findmax( indices->(length∘intersect)(indices,universe), indices_list )
    push!( new_indices_list, indices_list[pos] )
    setdiff!( universe, indices_list[pos] )
    deleteat!( indices_list, pos )
  end

  @show new_indices_list

  return new_indices_list
end # function greedy

###########################################
function make_complete_den_topology_collect(
  den_topology_list::Vector{DenTopology}
)::Vector{DenTopology}
###########################################
  n_loop = first(den_topology_list).n_loop
  ind_ext_mom = first(den_topology_list).ind_ext_mom
  n_sp::Int = (n_loop + 1) * n_loop / 2 + n_loop * length( ind_ext_mom )

  incomplete_den_topology_list = copy(den_topology_list)
  complete_den_topology_list = DenTopology[]

  while !isempty(incomplete_den_topology_list)
    cover_indices_list = (greedy∘get_cover_indices_list)( incomplete_den_topology_list )
    to_be_deleted_indices = Int[]
    for indices ∈ cover_indices_list
      this_den_topology_list = incomplete_den_topology_list[indices]
      to_be_complete_den_topology = union( this_den_topology_list... )

      while length(to_be_complete_den_topology) < n_sp
        println( "$to_be_complete_den_topology need to be completed." )
        @assert length(to_be_complete_den_topology) == (rank∘get_coeff_mat_mom2_sp)( to_be_complete_den_topology )
        vac_loop_mom_list = if n_loop ∈ 1:2
          (first∘get_vacuum_loop_momenta_list)( n_loop )
        elseif n_loop == 3
          mom_list = map( first∘get_args, to_be_complete_den_topology.den_list )
          map!( mom->subs(mom,Dict(ind_ext_mom.=>0)), mom_list, mom_list )
          unique!(mom_list)
          vac_loop_mom_list_collect = get_vacuum_loop_momenta_list( n_loop )
          selected_index = findfirst( vac_loop_mom_list->(isempty∘setdiff)(mom_list,vac_loop_mom_list), vac_loop_mom_list_collect )
          @assert !isnothing(selected_index)
          vac_loop_mom_list_collect[ selected_index ]
        end # if

        for ext_mom ∈ vcat( zero(Basic), ind_ext_mom ), q ∈ vac_loop_mom_list, the_sign ∈ [1,-1]
          trial_den = Basic( "Den( $(expand( q + the_sign * ext_mom )), 0, 0 )" )
          trial_topology = union( to_be_complete_den_topology, trial_den )
          rank_trial_topology = (rank∘get_coeff_mat_mom2_sp)(trial_topology)
          if rank_trial_topology > length(to_be_complete_den_topology)
            @show n_sp, rank_trial_topology, length(to_be_complete_den_topology)
            to_be_complete_den_topology = trial_topology
            rank_trial_topology == n_sp && break
          end # if
        end # for
        println()
      end # while

      @assert length( to_be_complete_den_topology ) == (rank∘get_coeff_mat_mom2_sp)( to_be_complete_den_topology ) == n_sp
      push!( complete_den_topology_list, to_be_complete_den_topology )
      union!( to_be_deleted_indices, indices )
    end # for indices
    deleteat!( incomplete_den_topology_list, sort!(to_be_deleted_indices) )
  end # while

  return complete_den_topology_list
end # function make_complete_den_topology_collect

function construct_den_topology( amp_dir::String )::Vector{DenTopology}
  @assert isdir(amp_dir)
  amp_file_list = filter( endswith(".jld2"), readdir( amp_dir; join=true, sort=false ) )
  sort!( amp_file_list;
    by=file_name->begin
      main_file_name = (last∘splitdir)(file_name)
      @assert (!isnothing∘match)( r"^amp[1-9]\d*.jld2$", main_file_name )
      parse( Int, match( r"[1-9]\d*", main_file_name ).match )
    end # file_name
  ) # end sort!
  @info "Found $(length(amp_file_list)) amplitude files at $amp_dir."

  n_loop, ext_mom_list = jldopen( first(amp_file_list), "r" ) do jld_file
    jld_file["n_loop"], to_Basic( jld_file["ext_mom_list"] )
  end # n_loop, ext_mom_list
  ind_ext_mom = ext_mom_list[1:end-1]

  den_topology_collect = [
    DenTopology( n_loop, ind_ext_mom, (to_Basic∘load)( amp_file, "loop_den_list" ) )
      for amp_file in amp_file_list
  ]
  backup_den_topology_collect = deepcopy( den_topology_collect )

  unique!( den_topology->reduce(*,den_topology.den_list), den_topology_collect )
  den_topology_collect = get_superior_den_topology_collect( den_topology_collect )
  @info "$(length(den_topology_collect)) topolgies found."

  complete_den_topology_collect = make_complete_den_topology_collect( den_topology_collect )
  @info "$(length(complete_den_topology_collect)) complete topologies found."

  for (index, complete_den_topology) ∈ enumerate( complete_den_topology_collect )
    @assert is_valid_dentopology(complete_den_topology)
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

import Base: intersect, isempty, issubset, length, setdiff, union

###########################################
mutable struct DenTop
  n_loop::Int
  ind_ext_mom::Vector{Basic}
  den_list::Vector{Basic}
end # mutable struct DenTop
###########################################




###########################################
function DenTop( 
    den_list::Vector{Basic}; 
    kwargs... 
)::DenTop
###########################################

  n_loop = haskey( kwargs, :n_loop ) ? kwargs[:n_loop] : (get_loop_index∘last∘get_loop_momenta)( den_list )
  ind_ext_mom = haskey( kwargs, :ind_ext_mom ) ? kwargs[:ind_ext_mom] : get_ext_momenta( den_list )
  return DenTop( n_loop, ind_ext_mom, den_list )

end # function DenTop

###########################################
function intersect( 
    dentop1::DenTop, 
    dentop2::DenTop 
)::DenTop
###########################################

  @assert dentop1.n_loop == dentop2.n_loop
  @assert dentop1.ind_ext_mom == dentop2.ind_ext_mom
  new_den_list = intersect( dentop1.den_list, dentop2.den_list )
  return DenTop( dentop1.n_loop, dentop1.ind_ext_mom, new_den_list )

end # function intersect

intersect( dentop::DenTop )::DenTop = 
    DenTop( dentop.n_loop, dentop.ind_ext_mom, intersect(dentop.den_list) )

intersect( dentop::DenTop, dentop_list... )::DenTop = 
    intersect( dentop, reduce( intersect, dentop_list ) )

isempty( dentop::DenTop ) = isempty( dentop.den_list )

issubset( dentop1::DenTop, dentop2::DenTop )::Bool = 
    dentop1.n_loop == dentop2.n_loop && 
    dentop1.ind_ext_mom == dentop2.ind_ext_mom && 
    issubset( dentop1.den_list, dentop2.den_list )

length( dentop::DenTop )::Int64 = length( dentop.den_list )

###########################################
function setdiff( 
    dentop1::DenTop, 
    dentop2::DenTop 
)::DenTop
###########################################
  @assert dentop1.n_loop == dentop2.n_loop
  @assert dentop1.ind_ext_mom == dentop2.ind_ext_mom
  new_den_list = setdiff( dentop1.den_list, dentop2.den_list )
  return DenTop( dentop1.n_loop, dentop1.ind_ext_mom, new_den_list )
end # function setdiff

###########################################
function union( 
    dentop1::DenTop, 
    dentop2::DenTop 
)::DenTop
###########################################
  @assert dentop1.n_loop == dentop2.n_loop
  @assert dentop1.ind_ext_mom == dentop2.ind_ext_mom
  new_den_list = union( dentop1.den_list, dentop2.den_list )
  return DenTop( dentop1.n_loop, dentop1.ind_ext_mom, new_den_list )
end # function union

###########################################
function union( 
    dentop::DenTop, 
    den_list::Vector{Basic} 
)::DenTop
###########################################

  @assert (get_loop_index∘last∘get_loop_momenta)( den_list ) ≤ dentop.n_loop
  @assert get_ext_momenta( den_list ) ⊆ dentop.ind_ext_mom

  return DenTop( dentop.n_loop, dentop.ind_ext_mom, union(dentop.den_list,den_list) )

end # function union

union( dentop::DenTop, den::Basic ) = union( dentop, [den] )

union( dentop::DenTop )::DenTop = DenTop( dentop.n_loop, dentop.ind_ext_mom, union(dentop.den_list) )

union( dentop::DenTop, dentop_list... )::DenTop = union( dentop, reduce( union, dentop_list ) )


###########################################
function is_valid_dentop( 
    dentop::DenTop 
)::Bool
###########################################

  n_loop = dentop.n_loop
  ind_ext_mom = dentop.ind_ext_mom
  den_list = dentop.den_list

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
end # function is_valid_dentop




###########################################
function get_vac_loop_momenta_list(
  ::Val{1} # n_loop
)::Vector{Vector{Basic}}
###########################################

  return [ [ Basic("q1") ] ]

end # function get_vac_loop_momenta_list

###########################################
function get_vac_loop_momenta_list(
  ::Val{2} # n_loop
)::Vector{Vector{Basic}}
###########################################

  return [ to_Basic( ["q1", "q2", "q1 + q2"] )] 

end # function get_vac_loop_momenta_list

###########################################
function get_vac_loop_momenta_list(
  ::Val{3} # n_loop
)::Vector{Vector{Basic}}
###########################################

  return [ to_Basic( ["q1", "q2", "q3", "q1 + q3", "q2 + q3", "q1 + q2 + q3"] ),
           to_Basic( ["q1", "q2", "q3", "q1 + q2", "q1 + q3", "q2 + q3"] ) ]

end # function get_vac_loop_momenta_list

###########################################
function gen_sp_dict(
  dentop::DenTop
)::Dict{Basic, Basic}
###########################################

  n_loop = dentop.n_loop
  n_ext_mom = length(dentop.ind_ext_mom)
  sp_index = 1
  sp_dict = Dict{Basic, Basic}()

  for ii ∈ 1:n_loop, jj ∈ ii:dentop.n_loop
    qi, qj = Basic("q$ii"), Basic("q$jj")
    sp_dict[ make_SP(qi, qj) ] = Basic("sp$(sp_index)")
    sp_index += 1
  end # for ii, jj

  for one_ext_mom ∈ dentop.ind_ext_mom, loop_ii ∈ 1:n_loop
    q = Basic("q$(loop_ii)")
    sp_dict[ make_SP(one_ext_mom, q) ] = Basic("sp$(sp_index)")
    sp_index += 1
  end # one_ext_mom, q

  @assert sp_index == n_loop * (n_loop + 1) / 2 + n_loop * n_ext_mom + 1

  return sp_dict

end # function gen_sp_dict

###########################################
function gen_vac_top( 
    dentop::DenTop 
)::DenTop
###########################################

  den_list = dentop.den_list
  ext_momenta = get_ext_momenta( den_list )
  vac_den_list = subs.( den_list, (ext_momenta .=> 0)... )
  unique!(vac_den_list)

  return DenTop( dentop.n_loop, dentop.ind_ext_mom, vac_den_list )

end # function gen_vac_top

###########################################
function get_coeff_mat_mom2_sp( 
    dentop::DenTop 
)::Matrix{Rational}
###########################################

  sp_dict = gen_sp_dict( dentop )
  n_sp = length(sp_dict)

  coeff_mat = zeros( Rational, length(dentop), n_sp )
  mom2_list = subs.( map( make_SP∘expand∘first∘get_args, dentop.den_list ), Ref(sp_dict) )

  for (mom2_index, mom2) ∈ enumerate(mom2_list), sp_index ∈ 1:n_sp
    the_coeff = SymEngine.coeff.( mom2, Basic("sp$(sp_index)") )
    coeff_mat[ mom2_index, sp_index ] = (Rational∘parse)( Int, string(the_coeff) )
  end # (mom2_index, mom2), sp_index
  
  return coeff_mat

end # function get_coeff_mat_mom2_sp

###########################################
function get_superior_dentop_collect(
  dentop_collect::Vector{DenTop}
)::Vector{DenTop}
###########################################

  new_dentop_collect = DenTop[]

  for dentop ∈ dentop_collect
    included_by_pos = findfirst( new_dentop -> dentop ⊆ new_dentop, new_dentop_collect )
    !isnothing(included_by_pos) && continue
    filter!( new_dentop -> new_dentop ⊈ dentop, new_dentop_collect )
    push!( new_dentop_collect, dentop )
  end # for dentop

  return new_dentop_collect

end # function get_superior_dentop_collect

###########################################
function get_cover_indices_list( 
    dentop_collect::Vector{DenTop} 
)::Vector{Vector{Int}}
###########################################
  n_loop = first( dentop_collect ).n_loop
  n_ind_ext = length( first(dentop_collect).ind_ext_mom )
  n_sp::Int = (n_loop + 1) * n_loop / 2 + n_loop * n_ind_ext

  n_dentop = length( dentop_collect )
  prev_indices_list = [ [ii] for ii ∈ 1:n_dentop ]

  for _ ∈ 2:n_dentop
    indices_list = Vector{Int}[]
    for one_indices ∈ prev_indices_list
      remaining_indices = setdiff( 1:n_dentop, one_indices )
      for one_index ∈ remaining_indices
        push!( indices_list, (sort∘union)( one_indices, one_index ) )
      end # for one_index
    end # for one_indices
    unique!( sort, indices_list )

    dentop_union_list = map( indices->union(dentop_collect[indices]...), indices_list )

    allowed_pos_list = findall(
      dentop -> length(dentop)≤n_sp&&(length∘gen_vac_top)(dentop)≤3*(n_loop-1),
      dentop_union_list
    ) # end findall

    indices_list = indices_list[ allowed_pos_list ]
    dentop_union_list = dentop_union_list[ allowed_pos_list ]
    fullrank_pos_list = findall(
      dentop->length(dentop)==(rank∘get_coeff_mat_mom2_sp)(dentop),
      dentop_union_list
    ) # end findall

    isempty(fullrank_pos_list) && return prev_indices_list
    prev_indices_list = indices_list[ fullrank_pos_list ]

  end # n_choice

  error("The `dentop_collect` should be covered by one topology.")

end # function get_cover_indices_list

###########################################
# greedy algorithm
function greedy( 
    indices_list::Vector{Vector{Int}} 
)::Vector{Vector{Int}}
###########################################
  @show indices_list

  universe = union( indices_list... )

  new_indices_list = Vector{Int}[]

  while !isempty(universe)
    _, pos = findmax( indices->(length∘intersect)(indices,universe), indices_list )
    push!( new_indices_list, indices_list[pos] )
    setdiff!( universe, indices_list[pos] )
    deleteat!( indices_list, pos )
  end # while

  @show new_indices_list

  return new_indices_list

end # function greedy

###########################################
function make_complete_dentop_collect(
  dentop_list::Vector{DenTop}
)::Vector{DenTop}
###########################################

  n_loop = first(dentop_list).n_loop
  ind_ext_mom = first(dentop_list).ind_ext_mom
  n_sp::Int = (n_loop + 1) * n_loop / 2 + n_loop * length( ind_ext_mom )

  incomplete_dentop_list = copy(dentop_list)
  complete_dentop_list = DenTop[]

  while !isempty(incomplete_dentop_list)
    cover_indices_list = (greedy∘get_cover_indices_list)( incomplete_dentop_list )
    to_be_deleted_indices = Int[]
    for indices ∈ cover_indices_list
      this_dentop_list = incomplete_dentop_list[indices]
      to_be_complete_dentop = union( this_dentop_list... )

      while length(to_be_complete_dentop) < n_sp
        println( "$to_be_complete_dentop need to be completed." )
        @assert length(to_be_complete_dentop) == (rank∘get_coeff_mat_mom2_sp)( to_be_complete_dentop )
        @assert n_loop in 1:3 
        vac_loop_mom_list = if n_loop ∈ 1:2
          (first∘get_vac_loop_momenta_list)( Val(n_loop) )
        elseif n_loop == 3
          mom_list = map( first∘get_args, to_be_complete_dentop.den_list )
          map!( mom->subs(mom,Dict(ind_ext_mom.=>0)), mom_list, mom_list )
          unique!(mom_list)
          vac_loop_mom_list_collect = get_vac_loop_momenta_list( Val(n_loop) )
          selected_index = findfirst( vac_loop_mom_list->(isempty∘setdiff)(mom_list,vac_loop_mom_list), vac_loop_mom_list_collect )
          @assert !isnothing(selected_index)
          vac_loop_mom_list_collect[ selected_index ]
        end # if

        for ext_mom ∈ vcat( zero(Basic), ind_ext_mom ), q ∈ vac_loop_mom_list, the_sign ∈ [1,-1]
          trial_den = Basic( "Den( $(expand( q + the_sign * ext_mom )), 0, 0 )" )
          trial_top = union( to_be_complete_dentop, trial_den )
          rank_trial_top = (rank∘get_coeff_mat_mom2_sp)(trial_top)
          if rank_trial_top > length(to_be_complete_dentop)
            @show n_sp, rank_trial_top, length(to_be_complete_dentop)
            to_be_complete_dentop = trial_top
            rank_trial_top == n_sp && break
          end # if
        end # for ext_mom
        println()
      end # while

      @assert length( to_be_complete_dentop ) == (rank∘get_coeff_mat_mom2_sp)( to_be_complete_dentop ) == n_sp
      push!( complete_dentop_list, to_be_complete_dentop )
      union!( to_be_deleted_indices, indices )
    end # for indices
    deleteat!( incomplete_dentop_list, sort!(to_be_deleted_indices) )
  end # while

  complete_dentop_list = get_superior_dentop_collect( complete_dentop_list )

  return complete_dentop_list

end # function make_complete_dentop_collect

###########################################
function construct_den_topology( 
    amp_dir::String 
)::Vector{DenTop}
###########################################

  @assert isdir(amp_dir) && endswith( amp_dir, "_amplitudes" )
  topology_name = begin
    tmp_dir, tmp_name = splitdir(amp_dir)
    tmp_dir = isempty(tmp_dir) ? pwd() : tmp_dir
    tmp_name = replace( tmp_name, "_amplitudes" => "_topology.out" )
    joinpath( tmp_dir, tmp_name )
  end # topology_name
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

  dentop_collect = [
    DenTop( n_loop, ind_ext_mom, (to_Basic∘load)( amp_file, "loop_den_list" ) )
      for amp_file in amp_file_list
  ]
  backup_dentop_collect = deepcopy( dentop_collect )

  unique!( dentop->reduce(*,dentop.den_list), dentop_collect )
  dentop_collect = get_superior_dentop_collect( dentop_collect )
  @info "$(length(dentop_collect)) topolgies found."

  complete_dentop_collect = make_complete_dentop_collect( dentop_collect )
  @info "$(length(complete_dentop_collect)) complete topologies found."

  file = open( topology_name, "w" )
  for (index, complete_dentop) ∈ enumerate( complete_dentop_collect )
    @assert is_valid_dentop(complete_dentop)
    pos_list = findall( dentop->dentop⊆complete_dentop, backup_dentop_collect )
    
    line_str = "-"^14
    println()
    println( line_str )
    println( "Complete topology #$(index) covers files:" )
    map( println, amp_file_list[ pos_list ] )
    println( line_str )
    map( println, complete_dentop.den_list )

    write( file, """
    $(line_str)
    Complete topology #$(index) covers files:
    $(join( map(string,amp_file_list[ pos_list ]), "\n" ))
    $(line_str)
    $(join( map(string,complete_dentop.den_list), "\n" ))

    """ )

  end # for (index, complete_dentop)
  close( file )

  box_message( "Information is in topology.out" )

  return complete_dentop_collect

end # function construct_den_topology




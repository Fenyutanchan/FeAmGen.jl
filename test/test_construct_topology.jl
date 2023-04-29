using FeAmGen

rootdir = dirname(@__FILE__)
amp_dir = joinpath( rootdir, "Wplus_t_TO_Wplus_t_3Loop_amplitudes" )

@assert isdir(amp_dir)
FeAmGen.construct_den_topology( amp_dir )



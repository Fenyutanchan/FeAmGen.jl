__precompile__()

"""
    module FeAmGen

FeAmGen represents `Fe`ynman diagram and `Am`plitude `Gen`erator.
"""
module FeAmGen

using AbstractAlgebra
using AmpTools
using Combinatorics
using Dates
using FORM_jll
using JLD2
using OrderedCollections
using Pipe
using Pkg
using PyCall
using SHA
using SymEngine
using Test
using YAML

export digest_seed_proc, generate_amp, generate_integral
export generate_multi_yaml, generate_shiftUP_yaml
export canonicalize_amp, gen_vac_reduction_ieta

include("Universe.jl")
include("Graph.jl")
include("Canon.jl")
include("Digest.jl")
include("FeynmanDiagram.jl")
include("FORMS.jl")
include("GenAmp.jl")
include("Integral.jl")
include("Kin.jl")
include("Seed.jl")
include("SimpleDigest.jl")
include("Visual.jl")


###################
function __init__()
###################
  return nothing
end # function __init__


end # module FeAmGen

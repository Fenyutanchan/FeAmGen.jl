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
using YAML

export digest_seed_proc
export generate_amp 
export canonicalize_amp 

include("Universe.jl")
include("Graph.jl")
include("Canon.jl")
include("Digest.jl")
include("FeynmanDiagram.jl")
include("FORMS.jl")
include("GenAmp.jl")
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

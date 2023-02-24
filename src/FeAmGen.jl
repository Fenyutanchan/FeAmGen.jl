__precompile__()

"""
    module FeAmGen

FeAmGen represents `Fe`ynman diagram and `Am`plitude `Gen`erator.
"""
module FeAmGen

using Combinatorics

using Downloads

# for reading input YAML file "SeedProcess.yaml"
using YAML

using OrderedCollections

using Pkg

# for read model python file
using PyCall

# for using @test and @testset
using Test

# for symbolic calculation
using SymEngine
using AmpTools

# for Feynman diagram
using Dates
using JLD2
using Pipe

using FORM_jll

export digest_seed_proc, generate_amp, generate_integral
export generate_multi_yaml, generate_shiftUP_yaml

include("Graph.jl")
include("Universe.jl")
include("SimpleDigest.jl")
include("Digest.jl")
include("Kin.jl")
include("FORMS.jl")
include("Visual.jl")
include("FeynmanDiagram.jl")
include("Seed.jl")
include("GenAmp.jl")
include("Integral.jl")

###################
function __init__()
###################
  return nothing
end # function __init__


end # module FeAmGen

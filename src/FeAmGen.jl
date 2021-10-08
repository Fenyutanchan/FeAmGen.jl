__precompile__()

"""
    module FeAmGen

FeAmGen represents `Fe`ynman diagram and `Am`plitude `Gen`erator.
"""
module FeAmGen

# for reading input YAML file "SeedProcess.yaml"
using YAML
using OrderedCollections
# for read model python file
using PyCall
# for using @test and @testset
using Test
# for symbolic calculation
using SymEngine
# for Feynman diagram
using Graphs
using Dates
using JLD
using Pipe

export digest_seed_proc, generate_amp, generate_integral, generate_multi_yaml, generate_shiftUP_yaml
export generate_SP_combo, gen_mma_str, box_message

include("Message.jl")
include("Extra.jl")
include("Universe.jl")
include("SimpleDigest.jl")
include("Digest.jl")
include("Kin.jl")
include("FORMS.jl")
include("Visual.jl")
include("Converter.jl")
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

__precompile__()

module FeAmGen

# for reading input YAML file "SeedProcess.yaml"
using YAML
# for read model python file
using PyCall
# for using @test and @testset
using Test
# for symbolic calculation
using SymEngine
# for Feynman diagram
using Graphs
using Dates

export digest_seed_proc, generate_amp

include("Message.jl")
include("Extra.jl")
include("Universe.jl")
include("SimpleDigest.jl")
include("Digest.jl")
include("Kin.jl")
include("FORMS.jl")
include("Visual.jl")
include("FeynmanDiagram.jl")
include("Seed.jl")
include("GenAmp.jl")

###################
function __init__()
###################

  show_welcome_message()

  return nothing
end # function __init__


end # module FeAmGen

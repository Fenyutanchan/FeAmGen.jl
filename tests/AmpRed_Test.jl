using SymEngine, FeAmGen, Test, BenchmarkTools, YAML, JLD, Pipe, Dates, Logging

io = open("AmpRed_Test.log", "w+")
logger = SimpleLogger(io)
global_logger(logger)

@info "AmpRed_Test starts @ $(now())"


#--------------------------------------------------------------------
# JLD file generation for single-top amplitude reduction.
#--------------------------------------------------------------------

yaml_str = """
name: Dia199SI1

n_loop: 2

min_eps_xpt: -4
max_eps_xpt: 0

external_momenta: [ k1, k2, K3, K4 ]

kin_relation: 
  - [ "SP(k1,k1)", "0" ]
  - [ "SP(k2,k2)", "0" ]
  - [ "SP(K3,K3)", "mw^2" ]
  - [ "SP(K4,K4)", "mt^2" ]
  - [ "SP(k1,k2)", "(1/2)*shat" ]
  - [ "SP(K3,k1)", "(1/2)*(-ver1 + mw^2)" ]
  - [ "SP(K4,k1)", "(1/2)*(shat + ver1 - mw^2)" ]
  - [ "SP(K3,k2)", "(1/2)*(shat + ver1 - mt^2)" ]
  - [ "SP(K3,K4)", "(1/2)*(shat - mt^2 - mw^2)" ]
  - [ "SP(K4,k2)", "(1/2)*(-ver1 + mt^2)" ]


den_list: [
"Den(q1,0,0)",                #D1
"Den(q2,mt,0)",               #D2
"Den(q1+q2,mt,0)",            #D3
"Den(q1+k1,0,0)",             #D4
"Den(q2-k1,mt,0)",            #D5
"Den(q2+k2-K3,0,0)",          #D6
"Den(q2-K3,0,0)",             #D7
"Den(q1-k2,0,0)",             #D8
"Den(q1+k1+k2-K3,0,0)"        #D9
]

den_xpt_list: [ 0, 0, 1, 0, 1, 0, 1, -2, 2 ]

numerator: "SP(q1,q2)"

comment: "For the tensor reduction of single-top amplitude."
"""

open( "scalar_integral.yaml", "w" ) do infile
  write( infile, yaml_str )
end 

generate_integral( "scalar_integral.yaml" )

rm( "scalar_integral.yaml" )

@testset "scalar_integral" begin

  content_dict = load( "integral_Dia199SI1.jld" )
  content_dict_bench = load( "integral_Dia199SI1_benchmark.jld" )

  @test content_dict == content_dict_bench 

end # @testset


@info "AmpRed_Test ends @ $(now())"

close(io)


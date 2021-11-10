using SymEngine, FeAmGen, Test, BenchmarkTools, YAML, JLD, Pipe, Dates, Logging

io = open("IRD_Test.log", "w+")
logger = SimpleLogger(io)
global_logger(logger)

@info "IRD_Test starts @ $(now())"

#--------------------------------------------------------------------
# Irreducible Integral (IRD)
#--------------------------------------------------------------------

IRD_origin_str = """
name: "IRD"

n_loop: 2

min_eps_xpt: -4
max_eps_xpt: 0

external_momenta: [ K1 ]

kin_relation:
  - [ "SP(K1,K1)", "mm^2" ]

den_list: [
"Den(q1,0,0)",
"Den(q1+K1,m0,0)",
"Den(q2,0,0)",
"Den(q2+K1,0,0)",
"Den(q1+q2,0,0)"
]

den_xpt_list: [ 0, 0, 0, 0, 0 ]

numerator: "1"

# 1: i*eta for all 
# 2: i*eta for only massive
ieta_scheme: 1

comment: "Seed yaml file for IRD"
"""

open( "IRD_original.yaml", "w" ) do infile
  write( infile, IRD_origin_str )
end 

# The last five (except the first one) integrals are master integarls.
indices_list = [ [0,1,1,-1,1], [0,1,3,-2,1] ]

multi_yaml_list = generate_multi_yaml( "IRD_original.yaml", indices_list, "IRD_integrals" )

rm( "IRD_original.yaml" )


for one_yaml in multi_yaml_list
  generate_integral( one_yaml ) 
end # for one_yaml

@testset "IRD" for one_yaml in multi_yaml_list
  file_name = last( splitpath( one_yaml )  )
  yaml_file0_dict = YAML.load_file( joinpath( "IRD_integrals_benchmark", file_name ) )
  yaml_file1_dict = YAML.load_file( joinpath( "IRD_integrals", file_name ) )
  delete!( yaml_file0_dict, "comment" )
  delete!( yaml_file1_dict, "comment" )
  delete!( yaml_file1_dict, "ieta_scheme" )
  @test yaml_file0_dict == yaml_file1_dict

  jld_name = "integral_$(file_name[1:end-5]).jld"
  jld_file0_dict = load( joinpath( "IRD_integrals_benchmark", jld_name ) )
  jld_file1_dict = load( joinpath( "IRD_integrals", jld_name ) )
  @test jld_file0_dict == jld_file1_dict
end # testset


@info "IRD_Test ends @ $(now())"

close(io)


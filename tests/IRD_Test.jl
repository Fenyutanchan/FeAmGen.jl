using SymEngine, FeAmGen, Test, YAML, JLD2, Dates

@info "IRD_Test starts @ $(now())"

#--------------------------------------------------------------------
# Irreducible Integral (IRD)
#--------------------------------------------------------------------

IRD_origin_str = """
name: "IRD"

n_loop: 2

min_ep_xpt: -4
max_ep_xpt: 0

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

# ieta_scheme 
# 0: none has iη 
# 1: all have iη
# 2: massive has iη
# >10: e.g. Int64(0b01010)+10, indexing position of iη via binary number
ieta_scheme: 1

comment: "Seed yaml file for IRD"
"""

#--------------------------------------------------------
# The integrals
indices_list = [ [0,1,1,-1,1], [0,1,3,-2,1] ]
# Generate original YAML file and then convert it into specific YAML files.
open( "IRD_original.yaml", "w" ) do infile
  write( infile, IRD_origin_str )
end # close
multi_yaml_list = generate_multi_yaml( "IRD_original.yaml", indices_list, "IRD_integrals" )
rm( "IRD_original.yaml" )
#--------------------------------------------------------


#-------------------------------
# Generate integral for each specific YAML file.
for one_yaml in multi_yaml_list
  generate_integral( one_yaml ) 
end # for one_yaml
#-------------------------------

#-------------------------------
# Compare to the prepared results.
@testset "IRD" begin
for one_yaml in multi_yaml_list
  file_name = last( splitpath( one_yaml )  )
  yaml_file0_dict = YAML.load_file( joinpath( "IRD_integrals_benchmark", file_name ) )
  yaml_file1_dict = YAML.load_file( joinpath( "IRD_integrals", file_name ) )
  delete!( yaml_file0_dict, "comment" )
  delete!( yaml_file1_dict, "comment" )
  @test yaml_file0_dict == yaml_file1_dict

  jld_name = "integral$(file_name[1:end-5]).jld2"
  jld_file0_dict = load( joinpath( "IRD_integrals_benchmark", jld_name ) )
  jld_file1_dict = load( joinpath( "IRD_integrals", jld_name ) )
  @test jld_file0_dict == jld_file1_dict
end # for one_yaml
end # testset


@info "IRD_Test ends @ $(now())"



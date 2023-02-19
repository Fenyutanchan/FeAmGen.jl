using SymEngine, FeAmGen, Test, YAML, JLD2, Dates 

@info "TSI_Test starts @ $(now())"

#--------------------------------------------------------------------
# Two-loop Self-energy Integral (TSI)
#--------------------------------------------------------------------

TSI_origin_str = """
name: "TSI"

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

comment: "Seed yaml file for TSI"
"""

#--------------------------------------------------------------------------
# The last five (except the first one) integrals are master integarls.
indices_list = [ [0,3,3,0,3], [0,1,2,0,1], [0,2,1,0,1], [0,1,1,0,1], [0,0,0,1,1], [0,1,0,0,1] ]
# Generate original YAML file and then convert it into specific YAML files.
open( "TSI_original.yaml", "w" ) do infile
  write( infile, TSI_origin_str )
end # close
multi_yaml_list = generate_multi_yaml( "TSI_original.yaml", indices_list, "TSI_integrals" )
rm( "TSI_original.yaml" )
#--------------------------------------------------------------------------


#-------------------------------
# Generate integral for each specific YAML file.
for one_yaml in multi_yaml_list
  generate_integral( one_yaml ) 
end # for one_yaml
#-------------------------------



#-------------------------------
# Generate YAML files for shifted integrals.
scalar_yaml_list = Vector{String}()
for one_indices in indices_list
  indices_str = join( map( string, one_indices ), "," )
  push!( scalar_yaml_list, "TSI_integrals/TSI_$(indices_str).yaml" )
end # for one_indices
new_yaml_list = generate_shiftUP_yaml( scalar_yaml_list, "TSI_shiftUP_integrals" )
#-------------------------------

#-------------------------------
# Generate integral for each specific YAML file.
for one_yaml in new_yaml_list
  generate_integral( one_yaml ) 
end # for one_yaml
#-------------------------------




@info "TSI_Test ends @ $(now())"


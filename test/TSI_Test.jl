using Dates, FeAmGen, AmpTools
using OrderedCollections, YAML

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
"Den(q1,0,ieta)",
"Den(q1+K1,m0,ieta)",
"Den(q2,0,ieta)",
"Den(q2+K1,0,ieta)",
"Den(q1+q2,0,ieta)"
]

den_xpt_list: [ 0, 0, 0, 0, 0 ]

numerator: "1"

comment: "Seed yaml file for TSI"
"""

#########################
function main()::Nothing
#########################

#--------------------------------------------------------------------------
# The last five (except the first one) integrals are master integarls.
indices_list = [ [0,3,3,0,3], [0,1,2,0,1], [0,2,1,0,1], [0,1,1,0,1], [0,0,0,1,1], [0,1,0,0,1] ]
# Generate original YAML file and then convert it into specific YAML files.
open( "TSI_original.yaml", "w" ) do infile
  write( infile, TSI_origin_str )
end # close
multi_yaml_list = generate_multi_yaml( "TSI_original.yaml", indices_list, "TSI_integrals" )
#--------------------------------------------------------------------------

#--------------------------------------------------------
file_dict = YAML.load_file( "TSI_original.yaml"; dicttype=OrderedDict{String,Any} ) 
loop_den_list = to_Basic( file_dict["den_list"] )
#--------------------------------------------------------
vac_top_list, vac_master_list = gen_vac_reduction_ieta( loop_den_list )
#--------------------------------------------------------
rm( "TSI_original.yaml" )
#--------------------------------------------------------


#-------------------------------
# Generate integral for each specific YAML file.
for one_yaml in multi_yaml_list
  generate_integral( one_yaml, vac_top_list, vac_master_list ) 
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
  generate_integral( one_yaml, vac_top_list, vac_master_list ) 
end # for one_yaml
#-------------------------------

return nothing

end # function main

########
main()
########


@info "TSI_Test ends @ $(now())"



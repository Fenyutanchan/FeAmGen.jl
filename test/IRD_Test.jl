using Dates, FeAmGen, AmpTools
using OrderedCollections, YAML

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
"Den(q1,0,ieta)",
"Den(q1+K1,m0,ieta)",
"Den(q2,0,ieta)",
"Den(q2+K1,0,ieta)",
"Den(q1+q2,0,ieta)"
]

den_xpt_list: [ 0, 0, 0, 0, 0 ]

numerator: "1"

comment: "Seed yaml file for IRD"
"""

##########################
function main()::Nothing
##########################

#--------------------------------------------------------
# The integrals
all_den_xpt_list = [ [0,1,1,-1,1], [0,1,3,-2,1] ]
# Generate original YAML file and then convert it into specific YAML files.
open( "IRD_original.yaml", "w" ) do infile
  write( infile, IRD_origin_str )
end # close
multi_yaml_list = gen_integral_multi_yaml( "IRD_original.yaml", all_den_xpt_list, "IRD_integrals" )
#--------------------------------------------------------

#--------------------------------------------------------
file_dict = YAML.load_file( "IRD_original.yaml"; dicttype=OrderedDict{String,Any} ) 
loop_den_list = to_Basic( file_dict["den_list"] )
#--------------------------------------------------------
vac_top_list, vac_master_list = gen_vac_reduction_ieta( loop_den_list )
#--------------------------------------------------------
rm( "IRD_original.yaml" )
#--------------------------------------------------------

#-------------------------------
# Generate integral for each specific YAML file.
for one_yaml in multi_yaml_list
  generate_integral( one_yaml, vac_top_list, vac_master_list ) 
end # for one_yaml
#-------------------------------

return nothing

end # function main

#######
main()
#######

@info "IRD_Test ends @ $(now())"



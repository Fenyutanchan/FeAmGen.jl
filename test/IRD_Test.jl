using Dates, FeAmGen

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

@info "IRD_Test ends @ $(now())"



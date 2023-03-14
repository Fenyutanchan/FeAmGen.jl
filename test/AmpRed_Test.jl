using Dates, FeAmGen, SymEngine, AmpTools 

@info "AmpRed_Test starts @ $(now())"

#--------------------------------------------------------------------
# JLD file generation for single-top amplitude reduction.
#--------------------------------------------------------------------

yaml_str( rank_str::String, index::Int64, num_str::String ) = """
name: Dia199_$(rank_str)_SI$(index)

n_loop: 2

min_ep_xpt: -4
max_ep_xpt: 0

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
"Den(q1,0,ieta)",                #D1
"Den(q2,mt,ieta)",               #D2
"Den(q1+q2,mt,ieta)",            #D3
"Den(q1+k1,0,ieta)",             #D4
"Den(q2-k1,mt,ieta)",            #D5
"Den(q2+k2-K3,0,ieta)",          #D6
"Den(q2-K3,0,ieta)",             #D7
"Den(q1-k2,0,ieta)",             #D8
"Den(q1+k1+k2-K3,0,ieta)"        #D9
]

den_xpt_list: [ 0, 0, 1, 0, 1, 0, 1, -2, 2 ]

numerator: "$(num_str)"

comment: "For the tensor reduction of single-top amplitude."
"""

@vars k1, k2, K3
rank_str = "q1q1q1"
num_list = generate_SPcombo( rank_str, [k1,k2,K3] )
num_list = sort( num_list, by=gen_sorted_str )
n_num = length( num_list )

target_dir = "AmpRed_integrals"
if isdir( target_dir ) || isfile( target_dir )
  mv( target_dir, "$(target_dir)_$(now())" )
end # if
mkdir( target_dir )


for index in 1:n_num
  one_num_str = string( num_list[index] )
  file_name = "$(target_dir)/scalar_integral_$(rank_str)_SI$(index).yaml"

  open( file_name, "w" ) do infile
    write( infile, yaml_str(rank_str,index,one_num_str) )
  end 

  generate_integral( file_name )
  rm( file_name )
end # for index


@info "AmpRed_Test ends @ $(now())"



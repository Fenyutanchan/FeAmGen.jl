using SymEngine, FeAmGen, Test, YAML, JLD2, Pipe, Dates

file = open( "eeHZ_mixed_proc.yaml", "w" )
write( file, """

# input file for calculation details
# model related information

# model name
model_name: "sm_CKMdiag_Haa"

# use unitary_gauge for internal vector boson
unitary_gauge: true

# content of "quark-parton", the "parton" will also contain gluon. Only for seed program.
#const partons = Array{String,1}( [ "g", "u", "d", "ubar", "dbar", "s", "c", "b", "sbar", "cbar", "bbar" ] )   
# for single top @ NNLO
#const partons = Array{String,1}( [ "g", "u", "d", "ubar", "dbar", "b", "bbar" ] )   
# for Higgs+Jet @ NLO
partons: [ "g", "u", "ubar", "d", "dbar", "b", "bbar" ] 

# only for seed program
AllowLeptonNumberViolation: false
AllowQuarkGenerationViolation: false

# process information
DropTadpole: true              # drop tadpole?
DropWFcorrection: true         # drop WFcorrection?

# number of loops
n_loop: 2    
# order of QCD counter-term vertices
QCDCT_order: 0   

# order of QCD coupling gs in the amplitude
Amp_QCD_order: 2 
# order of QED coupling ee in the amplitude
Amp_QED_order: 4 

# min ep power in the amplitude
Amp_Min_Ep_Xpt: -4
# max ep power in the amplitude
Amp_Max_Ep_Xpt: 0

# incoming and outgoing information
incoming: [ "eplus", "eminus" ]          # incoming particles
outgoing: [ "H", "Z" ]               # outgoing particles 

# whether to check the consistency between two versions of amplitudes
check_consistency: true

""")
close( file )

digest_seed_proc( "eeHZ_mixed_proc.yaml" )

generate_amp( "eplus_eminus_TO_H_Z_2Loop/eminus_eplus_TO_H_Z.yaml" )


using SymEngine, FeAmGen, Test, YAML, JLD2, Pipe, Dates, Downloads, SHA 

@info "eeHZ_Test starts @ $(now())"


#----------------------------------------------------------------------------
# ee->HZ 0-loop, 1-loop, 2-loop tests
#----------------------------------------------------------------------------
generic_eeHZ_seed_proc_yaml_str( ; nloop::Int64 = 2::Int64 ) = """
# input file for calculation details
# model related information

# model name
model_name: "sm"

# use unitary_gauge for internal vector boson
unitary_gauge: false

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
n_loop: $(nloop)
# order of QCD counter-term vertices
QCDCT_order: 0   

# order of QCD coupling gs in the amplitude
Amp_QCD_order: 0 
# order of QED coupling ee in the amplitude
Amp_QED_order: $(2+2*nloop)
# order of special coupling in the amplitude
Amp_SPC_order: 0 

# min ep power in the amplitude
Amp_Min_Ep_Xpt: $(-2*nloop)
# max ep power in the amplitude
Amp_Max_Ep_Xpt: 0

# incoming and outgoing information
incoming: [ "eminus", "eplus" ]          # incoming particles
outgoing: [ "H", "Z" ]               # outgoing particles 

# whether to check the consistency between two versions of amplitudes
check_consistency: false

"""

#-------------------------------
# Fetch the Model files.
if isdir("sm") && 
  (bytes2hex∘sha256)("sm") == 
      "5af308bec132bd499660594a6c127c57f910a25c18968eb4f9b976f48df8cb49"
  println( "sm has been found." )
else
  if isdir("sm") 
    rm("sm") 
  end # if

  url = "https://raw.githubusercontent.com/zhaoli-IHEP/FeAmGen_artifacts/main/Models/sm.tar.bz2"
  Downloads.download( url, "./sm.tar.bz2" )
  @assert (bytes2hex∘sha256)("sm.tar.bz2") == 
      "4ee796b534e9fd595ac4e12c52ae7847bc938f2c91d3a7daa4ff4c27de3d1204"
  run( `tar xjvf sm.tar.bz2` )
  println( "sm.tar.bz2 has been downloaded and decompressed." )
end # if

#-------------------------------------------
# Start running
for nloop in [0,1]

  open( "seed_eeHZ_proc_$(nloop)Loop.yaml", "w" ) do infile
    write( infile, generic_eeHZ_seed_proc_yaml_str(nloop=nloop) )
  end # close

  digest_seed_proc( "seed_eeHZ_proc_$(nloop)Loop.yaml" )

  generate_amp( "eminus_eplus_TO_H_Z_$(nloop)Loop/eminus_eplus_TO_H_Z.yaml" )

end # for nloop

@info "eeHZ_Test ends @ $(now())"



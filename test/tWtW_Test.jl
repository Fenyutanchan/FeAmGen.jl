using Dates, Downloads, FeAmGen, SHA

start = now()
@info "tWtW_Test starts @ $(start)"

#----------------------------------------------------------------------------
# t->b+W 1-loop, 2-loop including tbW vertices tests
#----------------------------------------------------------------------------
seed_proc_yaml_str( ; nloop::Int64 = 2::Int64 ) = """
# input file for calculation details
# model related information

# model name
model_name: "sm_tbW"

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
n_loop: $(nloop)
# order of QCD counter-term vertices
QCDCT_order: 0   

# order of QCD coupling gs in the amplitude
Amp_QCD_order: $(2*nloop) 
# order of QED coupling ee in the amplitude
Amp_QED_order: 0  
# order of special coupling in the amplitude
Amp_SPC_order: 2  

# min ep power in the amplitude
Amp_Min_Ep_Xpt: $(-2*nloop)
# max ep power in the amplitude
Amp_Max_Ep_Xpt: 0

# incoming and outgoing information
incoming: [ "Wplus", "t" ]          # incoming particles
outgoing: [ "Wplus", "t" ]          # outgoing particles 

# whether to check the consistency between two versions of amplitudes
check_consistency: false

"""

#-------------------------------
# Fetch the Model files.
if isdir("sm_tbW") && 
  (bytes2hex∘sha256)("sm_tbW") == 
      "fb20f86211878fa31a41cb0f18e9aa02db76cef6b61c27ebc8896ad40b943020"
  println( "sm_tbW has been found." )
else
  if isdir("sm_tbW") 
    rm("sm_tbW") 
  end # if

  url = "https://raw.githubusercontent.com/zhaoli-IHEP/FeAmGen_artifacts/main/Models/sm_tbW.tar.bz2"
  Downloads.download( url, "./sm_tbW.tar.bz2" )
  @assert (bytes2hex∘sha256)("sm_tbW.tar.bz2") == 
      "fe06b669291b1a9f7255e576bea1705ec24a6d238430c0d73a46b5909bab41d2"
  run( `tar xjvf sm_tbW.tar.bz2` )
  println( "sm_tbW.tar.bz2 has been downloaded and decompressed." )
end # if

#-------------------------------------------
# Start running
for nloop in [3]

  open( "seed_tWtW_proc_$(nloop)Loop.yaml", "w" ) do infile
    write( infile, seed_proc_yaml_str(nloop=nloop) )
  end # close

  digest_seed_proc( "seed_tWtW_proc_$(nloop)Loop.yaml" )

  generate_amp( "Wplus_t_TO_Wplus_t_$(nloop)Loop/Wplus_t_TO_Wplus_t.yaml" )

end # for nloop

@info "tWtW_Test ends @ $(now()) started from $(start)"



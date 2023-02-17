using SymEngine, FeAmGen, Test, YAML, JLD2, Dates, Pipe, Downloads, SHA 

start = now()
@info "tWb_Test starts @ $(start)"

#----------------------------------------------------------------------------
# t->b+W 1-loop, 2-loop including tbW vertices tests
#----------------------------------------------------------------------------
seed_proc_yaml_str( ; nloop::Int64 = 2::Int64 ) = """
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
n_loop: $(nloop)
# order of QCD counter-term vertices
QCDCT_order: 0   

# order of QCD coupling gs in the amplitude
Amp_QCD_order: $(2*nloop) 
# order of QED coupling ee in the amplitude
Amp_QED_order: 1  
# order of special coupling in the amplitude
Amp_SPC_order: 0  

# min ep power in the amplitude
Amp_Min_Ep_Xpt: $(-2*nloop)
# max ep power in the amplitude
Amp_Max_Ep_Xpt: 0

# incoming and outgoing information
incoming: [ "t" ]          # incoming particles
outgoing: [ "Wplus", "b" ]               # outgoing particles 

# whether to check the consistency between two versions of amplitudes
check_consistency: false

"""

#-------------------------------
# Fetch the Model files.
if isdir("sm_CKMdiag_Haa") && 
  (bytes2hex∘sha256)("sm_CKMdiag_Haa") == 
      "42890871318263aabd93bba66cfe8aa6013587e5ee6f8af607fa03f62652ee2e"
  println( "sm_CKMdiag_Haa has been found." )
else
  if isdir("sm_CKMdiag_Haa") 
    rm("sm_CKMdiag_Haa") 
  end # if

  url = "https://raw.githubusercontent.com/zhaoli-IHEP/FeAmGen_artifacts/main/Models/sm_CKMdiag_Haa.tar.bz2"
  Downloads.download( url, "./sm_CKMdiag_Haa.tar.bz2" )
  @assert (bytes2hex∘sha256)("sm_CKMdiag_Haa.tar.bz2") == 
      "926b05268a906a803cdb738349fd5b78351c8619996f1a16564f29a169543be9"
  run( `tar xjvf sm_CKMdiag_Haa.tar.bz2` )
  println( "sm_CKMdiag_Haa.tar.bz2 has been downloaded and decompressed." )
end # if

#-------------------------------------------
# Start running
for nloop in [3]

  open( "seed_tWb_proc_$(nloop)Loop.yaml", "w" ) do infile
    write( infile, seed_proc_yaml_str(nloop=nloop) )
  end # close

  digest_seed_proc( "seed_tWb_proc_$(nloop)Loop.yaml" )

  generate_amp( "t_TO_Wplus_b_$(nloop)Loop/t_TO_Wplus_b.yaml" )

end # for nloop

@info "tWb_Test ends @ $(now()) started from $(start)"



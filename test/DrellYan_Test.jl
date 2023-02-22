using Dates, Downloads, FeAmGen, SymEngineExt

@info "DrellYan_Test starts @ $(now())"

#----------------------------------------------------------------------------
# Drell-Yan 0-loop, 1-loop, 2-loop tests
#----------------------------------------------------------------------------
generic_dy_seed_proc_yaml_str( ; nloop::Int64 = 2::Int64 ) = """
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
n_loop: $nloop    
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
incoming: [ "dbar", "u" ]          # incoming particles
outgoing: [ "Wplus" ]               # outgoing particles 

# whether to check the consistency between two versions of amplitudes
check_consistency: true

"""

#-------------------------------
# Fetch the Model files.
if isdir("sm_CKMdiag_Haa") && 
  calc_sha256( filter( file_name->(last∘splitext)(file_name)==".py", readdir("sm_CKMdiag_Haa",join=true) ) ) == 
      "cfadd77f9c1383d50fbedada430174db68871f0a365a93fd6fa7ddfde6869c47"
  println( "sm_CKMdiag_Haa has been found." )
else
  if isdir("sm_CKMdiag_Haa") 
    rm("sm_CKMdiag_Haa") 
  end # if

  url = "https://raw.githubusercontent.com/zhaoli-IHEP/FeAmGen_artifacts/main/Models/sm_CKMdiag_Haa.tar.bz2"
  Downloads.download( url, "./sm_CKMdiag_Haa.tar.bz2" )
  @assert calc_sha256("sm_CKMdiag_Haa.tar.bz2") == 
      "a004f1e79ce4cfefb165f2560a8f10670ce9d420245e10ba4e37719ce51b7d3c"
  run( `tar xjvf sm_CKMdiag_Haa.tar.bz2` )
  println( "sm_CKMdiag_Haa.tar.bz2 has been downloaded and decompressed." )
end # if

#-------------------------------------------
# Start running
for nloop in [0,1,2]

  open( "seed_dy_proc_$(nloop)Loop.yaml", "w" ) do infile
    write( infile, generic_dy_seed_proc_yaml_str(nloop=nloop) )
  end # close

  digest_seed_proc( "seed_dy_proc_$(nloop)Loop.yaml" )

  generate_amp( "dbar_u_TO_Wplus_$(nloop)Loop/dbar_u_TO_Wplus.yaml" )

end # for nloop

@info "DrellYan_Test ends @ $(now())"



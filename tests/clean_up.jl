run( `rm dy_seed_proc_0Loop.yaml`, wait=false )
run( `rm -r dbar_u_TO_Wplus_0Loop`, wait=false )
run( `rm -r dbar_u_TO_Wplus_0Loop_amplitudes`, wait=false )
run( `rm -r dbar_u_TO_Wplus_0Loop_visuals`, wait=false )

run( `rm dy_seed_proc_1Loop.yaml`, wait=false )
run( `rm -r dbar_u_TO_Wplus_1Loop`, wait=false )
run( `rm -r dbar_u_TO_Wplus_1Loop_amplitudes`, wait=false )
run( `rm -r dbar_u_TO_Wplus_1Loop_visuals`, wait=false )

run( `rm dy_seed_proc_2Loop.yaml`, wait=false )
run( `rm -r dbar_u_TO_Wplus_2Loop`, wait=false )
run( `rm -r dbar_u_TO_Wplus_2Loop_amplitudes`, wait=false )
run( `rm -r dbar_u_TO_Wplus_2Loop_visuals`, wait=false )

run( `rm qgraf_out.dat`, wait=false )
run( `rm sm.log`, wait=false )
run( `rm sm_CKMdiag_Haa.log`, wait=false )


run( `rm gbtw_seed_proc_0Loop.yaml`, wait=false )
run( `rm -r g_b_TO_t_Wminus_0Loop`, wait=false )
run( `rm -r b_g_TO_Wminus_t_0Loop_amplitudes`, wait=false )
run( `rm -r b_g_TO_Wminus_t_0Loop_visuals`, wait=false )

run( `rm gbtw_seed_proc_1Loop.yaml`, wait=false )
run( `rm -r g_b_TO_t_Wminus_1Loop`, wait=false )
run( `rm -r b_g_TO_Wminus_t_1Loop_amplitudes`, wait=false )
run( `rm -r b_g_TO_Wminus_t_1Loop_visuals`, wait=false )

run( `rm gbtw_seed_proc_2Loop.yaml`, wait=false )
run( `rm -r g_b_TO_t_Wminus_2Loop`, wait=false )
run( `rm -r b_g_TO_Wminus_t_2Loop_amplitudes`, wait=false )
run( `rm -r b_g_TO_Wminus_t_2Loop_visuals`, wait=false )

run( `rm eeHZ_seed_proc_0Loop.yaml`, wait=false )
run( `rm -r eplus_eminus_TO_H_Z_0Loop`, wait=false )
run( `rm -r eminus_eplus_TO_H_Z_0Loop_amplitudes`, wait=false )
run( `rm -r eminus_eplus_TO_H_Z_0Loop_visuals`, wait=false )

run( `rm eeHZ_seed_proc_1Loop.yaml`, wait=false )
run( `rm -r eplus_eminus_TO_H_Z_1Loop`, wait=false )
run( `rm -r eminus_eplus_TO_H_Z_1Loop_amplitudes`, wait=false )
run( `rm -r eminus_eplus_TO_H_Z_1Loop_visuals`, wait=false )


run( `rm ggttbar_seed_proc_0Loop.yaml`, wait=false )
run( `rm -r g_g_TO_t_tbar_0Loop`, wait=false )
run( `rm -r g_g_TO_t_tbar_0Loop_amplitudes`, wait=false )
run( `rm -r g_g_TO_t_tbar_0Loop_visuals`, wait=false )

run( `rm ggttbar_seed_proc_1Loop.yaml`, wait=false )
run( `rm -r g_g_TO_t_tbar_1Loop`, wait=false )
run( `rm -r g_g_TO_t_tbar_1Loop_amplitudes`, wait=false )
run( `rm -r g_g_TO_t_tbar_1Loop_visuals`, wait=false )

rm( "TSI_integrals", recursive=true )
rm( "shiftUP_TSI_integrals", recursive=true )
rm( "IRD_integrals", recursive=true )



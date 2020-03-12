


########################################
function show_welcome_message()::Nothing
########################################

  MIRACLE_str = """
  ######################################################################################
  #                                                                                    #
  #      ____    ____  _____  _______          _        ______  _____     ________     #
  #     |_   \\  /   _||_   _||_   __ \\        / \\     .' ___  ||_   _|   |_   __  |    #
  #       |   \\/   |    | |    | |__) |      / _ \\   / .'   \\_|  | |       | |_ \\_|    #
  #       | |\\  /| |    | |    |  __ /      / ___ \\  | |         | |   _   |  _| _     #
  #      _| |_\\/_| |_  _| |_  _| |  \\ \\_  _/ /   \\ \\_\\ `.___.'\\ _| |__/ | _| |__/ |    #
  #     |_____||_____||_____||____| |___||____| |____|`.____ .'|________||________|    #
  #                                                                                    #
  #                                                                                    #
  #                  MultIple RAdiation Correction aLmighty gEnerator                  #
  #                                                                                    #
  #                               Version 0225.2017                                    #
  #                                                                                    #
  #                                   Zhao  Li                                         #
  #                                   IHEP-CAS                                         #
  #                               zhaoli@ihep.ac.cn                                    #
  #                                                                                    #
  ######################################################################################
  """
  
  # https://onlineasciitools.com/convert-text-to-ascii-art
  # Using "kban" "Fitted" font 
  SEDE_str = """
  ######################################################################################
  #                                                                                    #
  #                      .|'''.|'||''''| '||''|. '||''''|                              #
  #                      ||..  ' ||  .    ||   || ||  .                                #
  #                       ''|||. ||''|    ||    ||||''|                                #
  #                     .     '||||       ||    ||||                                   #
  #                     |'....|'.||.....|.||...|'.||.....|                             #
  #                                                                                    #
  #         SeDe: Julia program for Sector Decomposition @ project MIRACLE             #
  #                                                                                    #
  #                                   Zhao  Li                                         #
  #                                   IHEP-CAS                                         #
  #                               zhaoli@ihep.ac.cn                                    #
  #                                                                                    #
  ######################################################################################
  """
  FEAMGEN_str = """
  ######################################################################################
  #                                                                                    #
  #       '||''''|'||''''|     |    '||    ||'..|'''.|'||''''| '|.   '|'               #
  #        ||  .   ||  .      |||    |||  |||.|'     ' ||  .    |'|   |                #
  #        ||''|   ||''|     |  ||   |'|..'||||    ....||''|    | '|. |                #
  #        ||      ||       .''''|.  | '|' ||'|.    || ||       |   |||                #
  #       .||.    .||.....|.|.  .||..|. | .||.''|...'|.||.....|.|.   '|                #
  #                                                                                    #
  #     FeAmGen: Julia program for Feynman Amplitude Generation @ project MIRACLE      #
  #                                                                                    #
  #                                   Zhao  Li                                         #
  #                                   IHEP-CAS                                         #
  #                               zhaoli@ihep.ac.cn                                    #
  #                                                                                    #
  ######################################################################################
  """

  println(FEAMGEN_str)

  println("""

  Two functions are avaiable for user access.

  digest_seed_proc( "seed_proc.yaml", "../Models" )

  generate_amp( "parton_parton_TO_parton_t/b_u_TO_d_t.yaml", "../Models" )

  """)

  return nothing
end # function show_welcome_message





#####################################################################
function green_message( statement::String, message::String )::Nothing
#####################################################################

  print( statement )
  printstyled( message, "\n", color=:green )

  return nothing
end # function green_message



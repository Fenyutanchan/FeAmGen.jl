


########################################
function show_welcome_message()::Nothing
########################################

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



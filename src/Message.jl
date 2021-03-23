


########################################
"""
    show_welcome_message()

This function gives information for users.
Return nothing.
"""
function show_welcome_message()::Nothing
########################################

  println("""

  Two functions are avaiable for user access.

  digest_seed_proc( "seed_proc.yaml", "../Models" )

  generate_amp( "parton_parton_TO_parton_t/b_u_TO_d_t.yaml", "../Models" )

  """)

  return nothing

end # function show_welcome_message





#####################################################################
"""
    green_message( statement::String, message::String )::Nothing

Print the `statement` in normal color and `message` in green.
Return nothing.
"""
function green_message( statement::String, message::String )::Nothing
#####################################################################

  print( statement )
  printstyled( message, "\n", color=:green )

  return nothing
end # function green_message



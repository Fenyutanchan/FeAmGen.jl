







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



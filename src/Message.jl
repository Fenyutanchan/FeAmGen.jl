

#####################################################################
"""
    box_message( message::String; color=:light_cyan )::Nothing

Return nothing.
"""
function box_message( message::String; color=:light_cyan )::Nothing
#####################################################################

  len = length( message )
  bar = "-"^len
  printstyled( """

$bar
$message
$bar

""", color=color )

  return nothing
end # function box_message











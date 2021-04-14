diagram_index_list = [1, 2]
for diagram_index in diagram_index_list
  run( `lualatex visual_diagram$(diagram_index)` )
end

diagram_index_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
for diagram_index in diagram_index_list
  run( `lualatex visual_diagram$(diagram_index)` )
end

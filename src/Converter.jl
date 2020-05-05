
############################################################################################
function convert_couplingfactor( diagram_index::Int64, couplingfactor::Basic )::String
############################################################################################

  symbol_list = free_symbols( couplingfactor )
  n_symbol = length(symbol_list)

  replace_str_list = Vector{String}( undef, n_symbol )
  for symbol_index ∈ 1:n_symbol
    one_symbol = symbol_list[symbol_index]
    symbol_str = string(one_symbol)
    len = length(symbol_str)
    if symbol_str == "shat"
      replace_str_list[symbol_index] = "shat -> s"
    elseif symbol_str[1:5] == "gcsub"
      replace_str_list[symbol_index] = "$(symbol_str) -> Subscript[g,$(symbol_str[6:len])]"
    elseif symbol_str[1:3] == "ver"
      replace_str_list[symbol_index] = "$(symbol_str) -> Subscript[v,$(symbol_str[4:len])]"
    else
      error( "Exception: $(symbol_str)" )
    end # if
  end # for symbol_index

  file_name = "convert_couplingfactor_diagram$(diagram_index)"

  file = open( file_name*".m", "w" )
  write( file, """
expr = $(couplingfactor);
stream=OpenWrite["$(file_name).out"];
WriteString[ stream, expr//.{$(join( replace_str_list, "," ))}//TeXForm ];
Close[stream];
""" )
  close(file)

  #printstyled( "  run MathKernel -script $(file_name).m ... \n", color=:green )
  #println( "  Start @", Dates.now() )
  run( pipeline( `MathKernel -script $(file_name).m`, file_name*".log" ) )
  #println( "  Done @", Dates.now() )
  printstyled( "  Done MathKernel -script $(file_name).m in thread #$(Threads.threadid()) \n", color=:green )

  file = open( file_name*".out", "r" )
  result_str = replace( read( file, String ), r"\s"=>"" )
  close(file)

  rm( file_name*".m" )
  rm( file_name*".log" )
  rm( file_name*".out" )

  return result_str

end # function convert_couplingfactor










########################################################################
function convert_color_list( diagram_index::Int64, color_list::Vector{Basic} )::Vector{String}
########################################################################

  n_color = length(color_list)

  color_str_list = Vector{String}( undef, n_color )
  for color_index ∈ 1:n_color
    one_color = color_list[color_index]

    symbol_list = filter( s_ -> string(s_)[1:2] == "cl", free_symbols( one_color ) )
    symbol_str_list = map( string, symbol_list )
    replace_str_list = map( s_ -> "$(s_) -> Subscript[$(s_[3]),$(s_[4:length(s_)])]", symbol_str_list )

    color_mma_str = gen_mma_str(one_color)

    file_name = "convert_color_diagram$(diagram_index)_color$(color_index)"

    file = open( file_name*".m", "w" )
    write( file, """
expr = $(color_mma_str);
expr = expr//.{$(join(replace_str_list,","))};
expr = expr//.{im -> I, ca -> Subscript[C,A], cf -> Subscript[C,F], SUNT[x__] -> Subscript[T,x], SUNF[x__] -> Subscript[f,x] };
stream=OpenWrite["$(file_name).out"];
WriteString[ stream, expr//TeXForm ];
Close[stream];
""" )
    close(file)

    #printstyled( "  run MathKernel -script $(file_name).m ... \n", color=:green )
    #println( "  Start @", Dates.now() )
    run( pipeline( `MathKernel -script $(file_name).m`, file_name*".log" ) )
    #println( "  Done @", Dates.now() )
    printstyled( "  Done MathKernel -script $(file_name).m in thread #$(Threads.threadid()) \n", color=:green )
  
    file = open( file_name*".out", "r" )
    result_str = replace( read( file, String ), r"\s"=>"" )
    close(file)
  
    rm( file_name*".m" )
    rm( file_name*".log" )
    rm( file_name*".out" )

    color_str_list[color_index] = result_str

  end # for color_index

  return color_str_list

end # function convert_color_list




#####################################################################################################################################
function convert_lorentz_list( diagram_index::Int64, lorentz_list::Vector{Basic}, ext_mom_list::Vector{Basic}, scale2_list::Vector{Basic} )::Vector{String}
#####################################################################################################################################


  loop_mom_replace_str_list = String[ "q1 -> Subscript[q,1]", "q2 -> Subscript[q,2]", "q3 -> Subscript[q,3]" ]
  ext_mom_str_list = map( string, ext_mom_list )
  ext_mom_replace_str_list = map( s_ -> "$(s_) -> Subscript[$(s_[1]),$(s_[2])]", ext_mom_str_list )

  n_lorentz = length( lorentz_list )

  lorentz_str_list = Vector{String}( undef, n_lorentz )
  for lorentz_index ∈ 1:n_lorentz
    one_lorentz = lorentz_list[lorentz_index]

    dummy_symbol_list = filter( s_ -> length(string(s_)) > 7 && string(s_)[1:7] == "dummyMU", free_symbols(one_lorentz) )
    dummy_symbol_str_list = map( string, dummy_symbol_list )
    dummy_symbol_replace_str_list = map( s_ -> "$(s_) -> Subscript[ \\[Mu], $(s_[8:length(s_)]) ]", dummy_symbol_str_list )

    epsMu_symbol_list = filter( s_ -> length(string(s_)) > 5 && string(s_)[1:5] == "epsMU", free_symbols(one_lorentz) )
    epsMu_symbol_str_list = map( string, epsMu_symbol_list )
    epsMu_symbol_replace_str_list = map( s_ -> "$(s_) -> Subscript[ \\[Nu], $(s_[6:length(s_)]) ]", epsMu_symbol_str_list )

    gcsub_list = filter( s_ -> length(string(s_)) > 5 && string(s_)[1:5] == "gcsub", free_symbols(one_lorentz) )
    gcsub_str_list = map( string, gcsub_list )
    gcsub_replace_str_list = map( s_ -> "$(s_) -> Subscript[ g, $(s_[6:length(s_)]) ]", gcsub_str_list )

    scale_list = free_symbols( sum( scale2_list ) )
    scale_str_list = map( string, scale_list )
    n_scale = length( scale_str_list )
    scale_replace_str_list = Vector{String}( undef, n_scale )
    for scale_index ∈ 1:n_scale
      one_scale_str = scale_str_list[scale_index]
      if length(one_scale_str) > 3 && one_scale_str[1:3] == "ver"
        scale_replace_str_list[scale_index] = "$(one_scale_str) -> Subscript[s,$(one_scale_str[4:length(one_scale_str)])]"
      elseif one_scale_str[1:1] == "m"
        scale_replace_str_list[scale_index] = "$(one_scale_str) -> Subscript[M,$(one_scale_str[2:length(one_scale_str)])]"
      elseif one_scale_str == "shat"
        scale_replace_str_list[scale_index] = "shat -> s"
      else
        error( "Exception: $(one_scale_str)" )
      end # if
    end # for scale_index

    ver_replace_str_list = Vector{String}( undef, 10 )
    for ver_index ∈ 1:10
      ver_replace_str_list[ver_index] = "ver$(ver_index) -> Subscript[s,$(ver_index)]"
    end # for ver_index


    lorentz_mma_str = gen_mma_str(one_lorentz)

    file_name = "convert_lorentz_diagram$(diagram_index)_lorentz$(lorentz_index)"

    file = open( file_name*".m", "w" )
    write( file, """

expr = $(lorentz_mma_str);

MomList = {q1, q2, q3, $(join(map(string,ext_mom_list),","))};
vanishing = Map[# :> 0 &, MomList];

expr = expr //. FermionChain[ x1__, GA[mu_], x2__ ] * FV[mom_,mu_] :> FermionChain[ x1, GA[mom], x2 ];
expr = expr //. FermionChain[ x1__, GA[mom_/;Coefficient[mom,unity]=!=0], x2__ ] :> FermionChain[ x1, \\[Gamma][(mom/.unity:>0)]+Coefficient[mom,unity], x2 ];
expr = expr //. DiracTrace[ x1___, GA[mom_/;Coefficient[mom,unity]=!=0], x2___ ] :> DiracTrace[ x1, \\[Gamma][(mom/.unity:>0)]+Coefficient[mom,unity], x2 ];
expr = expr //. FermionChain[ x1__, GA[mom_/; (mom /. vanishing) == 0], x2__ ] :> FermionChain[ x1, \\[Gamma][mom], x2 ];

expr = expr //. FermionChain[x__] :> Dot[x] //. { PR -> Subscript[P, R], PL -> Subscript[P, L], im -> I, 
                  UB[idx_, x__] -> Subscript[U, idx], U[idx_, x__] -> Subscript[u, idx], 
                  VB[idx_, x__] -> Subscript[V, idx], V[idx_, x__] -> Subscript[v, idx],
                  VecEps[idx_,mu_,x__] -> Superscript[Subscript[\\[Epsilon],idx],mu], VecEpsC[idx_,mu_,x__] -> Superscript[Subscript[\\[Epsilon],idx,c],mu],
                  FV[k1_,mu1_]*FV[k2_,mu2_]*LMT[mu1_,mu2_] -> SP[k1,k2] };

expr = expr //.{$(join(loop_mom_replace_str_list,","))} //.{$(join(ext_mom_replace_str_list,","))};

expr = expr //. {$(join(dummy_symbol_replace_str_list,","))};
expr = expr //. {$(join(epsMu_symbol_replace_str_list,","))};
expr = expr //. {$(join(gcsub_replace_str_list,","))};
expr = expr //. {$(join(scale_replace_str_list,","))};
expr = expr //. {$(join(ver_replace_str_list,","))};
expr = expr //. GA[x_ /; MemberQ[MomList, x]] -> \\[Gamma].x;
expr = expr //. GA[x_ /; ! MemberQ[MomList, x]] -> Subscript[\\[Gamma], x];
expr = expr //. diim -> d //. SP[x1_,x2_] -> Subscript[S,P][x1,x2] //. DiracTrace[x__] -> Subscript[D,T][x] //. FV[mom_,mu_] -> Subscript[F,4][mom,mu];

stream=OpenWrite["$(file_name).out"];
WriteString[ stream, expr//TraditionalForm//TeXForm ];
Close[stream];
""" )
    close(file)

    #printstyled( "  run MathKernel -script $(file_name).m ... \n", color=:green )
    #println( "  Start @", Dates.now() )
    run( pipeline( `MathKernel -script $(file_name).m`, file_name*".log" ) )
    #println( "  Done @", Dates.now() )
    printstyled( "  Done MathKernel -script $(file_name).m in thread #$(Threads.threadid()) \n", color=:green )
  
    file = open( file_name*".out", "r" )
    result_str = replace( read( file, String ), r"\s"=>"" )
    close(file)
  
    rm( file_name*".m" )
    rm( file_name*".out" )
    rm( file_name*".log" )

    result_str = replace( result_str, "U"=>"\\bar{u}" )
    result_str = replace( result_str, "V"=>"\\bar{v}" )

    lorentz_str_list[lorentz_index] = result_str

  end # for lorentz_index

  return lorentz_str_list

end # function convert_lorentz_list















using AmpTools
using FORM_jll
using JLD2
using Pkg
using SymEngine

##############################
function calc_color_square(
    one_color::Basic,
    conj_color::Basic 
)::Basic
##############################

  file = open( "calc.frm", "w" )
  write( file, """

  #-
  Off Statistics;

  format nospaces;
  format maple;

  symbols nc;

  #include color.frm

  Local colorConj = $(conj_color);

  *** make conjugate for colorConj only
  id SUNT = SUNTConj;
  id sunTrace = sunTraceConj;
  .sort

  Local colorFactor = $(one_color);
  .sort

  Local colorSquare = colorConj*colorFactor;
  .sort
  drop colorConj;
  drop colorFactor;
  .sort

  #call calc1_CF();
  .sort 

  #call calc2_CF();
  .sort 

  id ca = nc;
  id cf = (nc^2-1)/(2*nc);
  id ca^(-1) = nc^(-1);
  id ca^(-2) = nc^(-2);
  .sort

  #write <calc.out> "%E", colorSquare
  #close <calc.out>
  .sort

  .end

  """ )
  close( file )

  run( pipeline( `$(form()) calc.frm`, "calc.log" ) )

  file = open( "calc.out", "r" )
  result_expr = (Basic∘read)( file, String )
  close( file )

  rm( "calc.frm" )
  rm( "calc.out" )
  rm( "calc.log" )
  
  return result_expr

end # function calc_color_square





#########################
function main()::Nothing
#########################

  art_dir = Pkg.Artifacts.artifact"FeAmGen"
  cp( "$(art_dir)/scripts/color.frm", "color.frm", force=true )

  amp_dir = "./Wplus_t_TO_Wplus_t_3Loop_amplitudes"

  root, dirs, files = (first∘collect∘walkdir)(amp_dir)
  jld_list = filter( s->endswith(s,".jld2"), files )

  @vars ca cf nc

  conj_color = Basic("DeltaFun(cla4, clb2)")
  unique_color_list = Vector{Basic}()
  for file_name in jld_list
    file = jldopen( "$(amp_dir)/$(file_name)", "r" ) 
    color_list = (to_Basic∘read)( file, "amp_color_list" )
    close( file )

    color_square_list = map( x->calc_color_square(x,conj_color), color_list )
    color_square_list = map( x->subs(x,Basic("im")=>im), color_square_list )
    #println( "[ $(file_name) ]" )
    for one_square in color_square_list
      nc4_coeff = coeff(one_square,nc,Basic(4))
      if !iszero(nc4_coeff) 
        println( "[ $(file_name) ]" )
        println( "  $(one_square)" )
      end # if
      term_list = get_add_vector_expand(one_square)
      filtered_term_list = filter( x->SymEngine.get_symengine_class(x)∉[:Integer,:Rational,:Complex],
                                   (unique∘map)( drop_coeff, term_list ) )
      union!( unique_color_list, filtered_term_list )
    end # for one_square
  end # for file_name

  @show unique_color_list

  rm( "color.frm" )

  return nothing

end # function main


########
main()
########

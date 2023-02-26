
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

trace5_index_list = [2, 3, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 31, 46, 47, 50, 51, 52, 53, 62, 68, 69, 72, 73, 76, 77, 78, 79, 329, 331, 349, 350, 357, 358, 371, 372, 379, 380, 387, 388, 395, 396, 417, 419, 421, 423, 437, 439, 441, 443, 488, 490, 522, 523, 936, 938, 939, 953, 954, 955, 956, 988, 989, 996, 997, 1004, 1005, 1012, 1013, 1024, 1025, 1026, 1027, 1040, 1041, 1042, 1043, 1098, 1102, 1106, 1110, 1150, 1152, 1154, 1156, 1204, 1213, 1270, 1271, 1278, 1279, 1292, 1293, 1300, 1301, 1332, 1333, 1334, 1335, 1336, 1337, 1338, 1339, 1341, 1354, 1355, 1356, 1357, 1358, 1359, 1360, 1361, 1363, 1378, 1379, 1386, 1387, 1442, 1443, 1450, 1451, 1794, 1795, 1797, 1798, 1804, 1805, 1811, 1812, 1813, 1814, 1818, 1819, 1822, 1823, 1825, 1826, 1827, 1828, 1829, 1830, 1831, 1832, 1841, 1842, 1846, 1847, 1848, 1850, 1855, 1856, 1857, 1858, 1859, 1861, 1862, 1863, 1873, 1874, 1876, 1877, 1878, 1879, 1882, 1884, 1893, 1894, 1895, 1896, 1923, 1924, 1925, 1926, 1927, 1928, 1929, 1930, 1939, 1940, 1942, 1943, 1944, 1945, 1946, 1947, 1950, 1951, 1952, 1953, 1961, 1962, 1963, 1964, 1975, 1976, 1978, 1979, 1989, 1990, 1991, 1992, 2005, 2006, 2008, 2009, 2010, 2011, 2012, 2013, 2023, 2024, 2026, 2027, 2034, 2035, 2036, 2037, 2044, 2045, 2046, 2047, 2058, 2059, 2060, 2061, 2062, 2063, 2064, 2065, 2066, 2067, 2068, 2069, 2070, 2071, 2072, 2073, 2084, 2085, 2086, 2087, 2088, 2089, 2090, 2091, 2092, 2093, 2094, 2095, 2096, 2097, 2098, 2099, 2110, 2111, 2112, 2113, 2134, 2135, 2136, 2137, 2158, 2159, 2160, 2161, 2162, 2163, 2165, 2166, 2178, 2179, 2180, 2181, 2182, 2183, 2184, 2185, 2186, 2188, 2189, 2190, 2195, 2196, 2197, 2198, 2199, 2200, 2203, 2204, 2226, 2227, 2228, 2229, 2230, 2231, 2232, 2233, 2236, 2237, 2238, 2239, 2240, 2241, 2242, 2243, 2244, 2245, 2246, 2247, 2248, 2249, 2250, 2251, 2252, 2253, 2254, 2255, 2256, 2257, 2258, 2259, 2261, 2262, 2263, 2264, 2265, 2266, 2267, 2268, 2287, 2288, 2290, 2291, 2292, 2293, 2294, 2295, 2296, 2297, 2298, 2299, 2301, 2302, 2303, 2304, 2305, 2306, 2307, 2308, 2309, 2310, 2311, 2312, 2313, 2314, 2315, 2316, 2317, 2318, 2319, 2320, 2321, 2322, 2323, 2324, 2344, 2345, 2346, 2347, 2348, 2349, 2350, 2359, 2370, 2374, 2378, 2382, 2433, 2434, 2435, 2450, 2451, 2452, 2470, 2539, 2540, 2547, 2548, 2555, 2556, 2563, 2564]


  art_dir = Pkg.Artifacts.artifact"FeAmGen"
  cp( "$(art_dir)/scripts/color.frm", "color.frm", force=true )

  amp_dir = "./Wplus_t_TO_Wplus_t_3Loop_amplitudes"

  root, dirs, files = (first∘collect∘walkdir)(amp_dir)
  #jld_list = filter( s->endswith(s,".jld2"), files )
  jld_list = map( x->"amplitude_diagram$(x).jld2", trace5_index_list )

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
      nc4_coeff = coeff(one_square,nc,Basic(3))
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

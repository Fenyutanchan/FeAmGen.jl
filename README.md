# FeAmGen.jl

FeAmGen: Julia program for Feynman Amplitude Generation      

-----------------------------------------------

Install this package.

```julia
(@v1.8) pkg> add https://github.com/zhaoli-IHEP/FeAmGen.jl.git
```

-----------------------------------------------

```julia
digest_seed_proc( "seed_proc.yaml", "Models" )
```
`digest_seed_proc( <seed process file>, <model directory> )`
This is used to generate process input files according to generic seed process file and model directory.

-----------------------------------------------

```julia
generate_amp( "parton_parton_TO_parton_t_0Loop/b_u_TO_d_t.yaml", "Models" )
```
`generate_amp( <process file>, <model directory> )`
This is used to generate diagrams and amplitudes by using the previously generated process input file and model directory.

-----------------------------------------------

`generate_integral( <YAML file> )`
This is used to generate expression for given YAML file, which has specific format shown in the test examples.

-----------------------------------------------


Here we use directly the UFO model files stored in the `<model directory>`.

***NOTE:***
***If the UFO model files are based on Python2, one can convert them into Python3 scripts by using the script `2to3`, which may be found in for example `Python-3.8.1/Tools/scripts`***
***And PyCall can be build with relevant version of Python. The method can be found in PyCall.jl homepage.***

And we need `QGRAF` and `FORM` packages installed.

The results contains amplitude in file `amplitude_diagram<number>.out` and the Feynman diagrams in file `visual_diagram<number>.tex`.
Explicitly one could use `lualatex visual_diagram<number>.tex` to generate PDF file with `tikz-feynman.sty`.

------------------------------------------------

Documentation can be generated by using `julia make.jl` in `docs`.




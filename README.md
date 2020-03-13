# FeAmGen.jl

FeAmGen: Julia program for Feynman Amplitude Generation @ project MIRACLE      

-----------------------------------------------

```julia
digest_seed_proc( "seed_proc.yaml", "Models" )
```
`digest_seed_proc( <seed process file>, <model directory> )`
This is used to generate process input files according to generic seed process file and model directory.

-----------------------------------------------

```julia
generate_amp( "parton_parton_TO_parton_t/b_u_TO_d_t.yaml", "Models" )
```
`generate_amp( <process file>, <model directory> )`
This is used to generate diagrams and amplitudes by using the previously generated process input file and model directory.

-----------------------------------------------

Here we use directly the UFO model files stored in `<model directory>`.

***NOTE: Since we are using UFO format, PyCall.jl needs to be compiled with Python2 instead of Python3.***

And we need `QGRAF` and `FORM` packages installed.

The results contains amplitude in file `amplitude.out` and the Feynman diagrams in file `visual_graphs.tex`.
Explicitly one could use `lualatex visual_graphs.tex` to generate PDF file with `tikz-feynman.sty`.


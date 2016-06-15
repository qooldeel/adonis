#!/bin/bash

# $1 file name for reduced ODE data
# $2 extract extract column from $1
# $3 file name for full ODE data
# $4 extract extract column from $3
# $5 output name

LAB="\$Y_\{\\\\ce\{O\}\}\$"
#LAB="\$T\$"
KEY="set key center right"
#"set xtics (0, .02, .04, .06, .08, .1)"

#draw x-axis tic every 0.02
XTICS="set xtics 0.02" 

COLOR1="#5F9EA0"
COLOR2="#DC143C"

echo "Create .eps and .tex files..."
gnuplot -persist <<PLOT
set terminal epslatex size 8.5 cm, 7.0 cm color colortext 
set output '$5.tex'
set format xy "$%g$"
${KEY}
${XTICS}

set mxtics 2 
set mytics 2
set grid xtics ytics mxtics mytics
set xlabel "\$t\$"
plot '$1' using 1:$2 with lines lw 2 lc rgbcolor "${COLOR1}" title "${LAB} (red)", '$3' using 1:$4 with lines lw 2 lc rgbcolor "${COLOR2}" title "${LAB} (full)"

quit
PLOT
echo "Done ..."


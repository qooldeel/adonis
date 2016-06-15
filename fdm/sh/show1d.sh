#!/bin/bash 

# $1     #filename RED 
#  $2   #column from $1 to be printed 
# $3   # filename FULL
# $4   # column from $2 to be printed

# $5  # only when plotting to epslatex is demanded

LAB="\$T\$"

COLOR1="#5F9EA0"
COLOR2="#DC143C"

if [ $2 -eq 1 ]; then
    echo -e "This is the x-coordinate. Plot from col 2 onwards only \nEXIT"
    exit 113
fi

if [ $4 -eq 1 ]; then
    echo -e "This is the x-coordinate. Plot from col 2 onwards only \nEXIT"
    exit 113
fi

gnuplot -persist <<PLOT
set terminal epslatex size 8.5 cm, 7.0 cm color colortext 
set output '$5.tex'
set format xy "$%g$"

set mxtics 2 
set mytics 2
set grid xtics ytics mxtics mytics

set xlabel "\$x\$"
plot '$1' using 1:$2 with lines lw 2 lc rgbcolor "${COLOR1}" title "${LAB} (red)", '$3' using 1:$4 with lines lw 2 lc rgbcolor "${COLOR2}" title "${LAB} (full)"
quit
PLOT
echo "Done ..."
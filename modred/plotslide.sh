#!/bin/bash 

FILE=$1

gnuplot -persist <<PLOT
set xlabel "x"
set ylabel "concentrations"
set mxtics 4
set mytics 4
set grid xtics ytics mxtics mytics
plot '${FILE}' using 1:2 title "H2" with lines lw 2, '${FILE}' using 1:3 with lines lw 2 title "H2O"

quit
PLOT
echo "Done ..."
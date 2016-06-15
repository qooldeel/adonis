#!/bin/bash

TIME=1   #by default this column contains time points

# $1  file to plot
# $2  column to be plot against time

gnuplot -persist <<PLOT

set mxtics 2 
set mytics 2
set grid xtics ytics mxtics mytics

plot '$1' using ${TIME}:$2 with lines lw 2 

quit
PLOT
echo "Done ..."

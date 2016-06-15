#!/bin/bash

## assume the matrix has been previously output in the usual C-style, i.e.
## looping first over rows then over columns.
## Then this scripts plots the matrix with GNUplot
## Taken from 
##    http://stackoverflow.com/questions/8299103/plotting-a-correlation-matrix-with-gnuplot

FILE=$1

gnuplot -persist <<PLOT
set terminal postscript eps enhanced color
set output "SparsitypatternOfJacobian2D_FDM.eps"


## this seems not to work properly
 # set pm3d map
 # splot "$FILE" matrix

##alternative:
set xlabel "# columns"
set ylabel "# rows"
set title "Jacobian pattern of 2D MOL"
plot "$FILE" matrix with image


quit
PLOT
echo "Done ..."
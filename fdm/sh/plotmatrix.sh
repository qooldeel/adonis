#!/bin/bash

## assume the matrix has been previously output in the usual C-style, i.e.
## looping first over rows then over columns.
## Then this scripts plots the matrix with GNUplot
## Taken from 
##    http://stackoverflow.com/questions/8299103/plotting-a-correlation-matrix-with-gnuplot

FILE=$1
FONTSIZE=23

FNAME=$(basename $FILE)
EPSNAME=${FNAME%.*}
echo -e "\n EPS file created: $EPSNAME.eps "
gnuplot -persist <<PLOT
##places the origin at the top-left corner. This produces images for 2D FDM
## like in [Y. SAAD, "Iterative Methods for Sparse Linear Systems", p. 55]
set yrange [:] reverse
## no extensions to next tics
set autoscale fix
set terminal postscript eps enhanced color font "Times-Roman, $FONTSIZE"
#set output "SparsitypatternOfJacobian2D_FDM.eps"
set output "$EPSNAME.eps"

## this seems not to work properly
 # set pm3d map
 # splot "$FILE" matrix 

##alternative:
# show xtics and xlabel on top, don't show anything at bottom
unset xtics
set x2tics 
set x2label "# columns"

set ylabel "# rows"
set title "Jacobian pattern of 2D MOL"
plot "$FILE" matrix with image notitle


quit
PLOT
echo " Done ..."
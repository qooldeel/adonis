#!/bin/bash 

FILE=$1     #filename 
COLUMN=$2   #column to be printed against columns 1 and 2 (the x and y-axes) 
gnuplot -persist <<PLOT

# set title "temperature profile"

# wanna print to postscript??? ...
# set terminal postscript eps enhanced color
# set output "messdat.eps"

set pm3d map          #at b
set palette rgbformulae 33,13,10
# set samples 50, 50
# set isosamples 50, 50
splot "$FILE" using 1:2:$2
unset pm3d

quit
PLOT
echo "Done ..."
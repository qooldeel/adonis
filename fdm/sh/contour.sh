#!/bin/bash 

## For interpolation see: 
## http://www.gnuplotting.org/tag/pm3d/

FILE=$1     #filename 
COLUMN=$2   #column to be printed against columns 1 and 2 (the x and y-axes) 

echo "---------------------------------"
echo "CONTOUR PLOT FOR 2D NAVIER-STOKES"
echo "---------------------------------"

if [ $2 -eq 1 ]; then
    echo -e "This is the x-coordinate. Plot from col 3 onwards only \nEXIT"
    exit 113
elif [ $2 -eq 2 ]; then 
    echo -e "This is the y-coordinate. Plot from col 3 onwards only \nEXIT"
    exit 113
else
    echo "column o.k."
    # if [ $2 -eq 3 ];  then
    # 	echo -e "plot rho \n"
    # elif [ $2 -eq 4 ]; then
    # 	echo -e "plot v1 \n"
    # elif [ $2 -eq 5 ]; then
    # 	echo -e "plot v2 \n"
    # elif [ $2 -eq 6 ]; then
    # 	echo -e "plot T \n"
    # elif [ $2 -ge 7 ]; then
    # 	res=$(($2 - 7))
    # 	echo -e "plot species Y_$res \n"
    # else
    # 	echo "unspecified column"
    # fi
fi

#gnuplot 
gnuplot -persist <<PLOT

# set title "temperature profile"

#scientific format of x- and y-axis tics, respectively
#set format x "%.1e"
#set format y "%.1e"


# set xtics 0,.0003125,0.02
# set ytics  0,.0003125,0.0025

#set grid

set xlabel "x"
set ylabel "y"

## tell gnuplot to disable extensions of the axis range to the next tic:
## see http://gnuplot.sourceforge.net/docs_4.2/node157.html
set autoscale fix

# wanna print to postscript??? ...
# set terminal postscript eps enhanced color
# set output "messdat.eps"

set title "${FILE}"

## for interpolation: n,m are the additional points along the x- and y-axis
## setting n = m = 0 forces gnuplot to choose the correct number of interpol.
## points by itself
###TRY this:            width,height  
#set size ratio 0.1 #make graph 1/10 as high as wide
set pm3d map interpolate 0,0       #at b
set palette rgbformulae 33,13,10
# set samples 50, 50
 set isosamples 100, 100
splot "$FILE" using 1:2:$2 notitle
unset pm3d

quit
PLOT
echo "Done ..."

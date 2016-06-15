#!/bin/bash

FULL=$1            #$1 FULL Solution 
RED1st=$2          #$2 1st REDUCTION
OUTFILE=$3         #$3 Outfile  (WITHOUT suffix)

#define some fancy colors here ;)
# if you're not sure, visit 
#  http://www.uni-hamburg.de/Wiss/FB/15/Sustainability/schneider/gnuplot/colors.htm
 COLOR1="#5F9EA0"
 COLOR2="#DC143C"
 COLOR3="#8A2BE2"
 COLOR4="#FF7F50"


## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## CAUTION $ is a reserved character. If you want to use it, then type \$
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
LABH2="\$\\\\conc\{\\\\ce\{H2\}\}\$"
LABH2O="\$\\\\conc\{\\\\ce\{H2O\}\}\$"

#KEY="set key top right"
KEY="set key center right"

echo "Create .eps and .tex files..."
gnuplot -persist <<PLOT
set terminal epslatex size 8.5 cm, 7.0 cm color colortext 
set output '${OUTFILE}.tex'
set format xy "$%g$"
set mxtics 2 
set mytics 2
set grid xtics ytics mxtics mytics
set xlabel "\$x\$"
${KEY}
plot '${FULL}' using 1:2  with lines lw 2 lc rgbcolor "${COLOR1}" title "${LABH2} (full)", '${FULL}' using 1:3 with lines lw 2 lc rgbcolor "${COLOR2}" title "${LABH2O} (full)", '${RED1st}' using 1:2  with lines lw 2 lc rgbcolor "${COLOR3}"  title "${LABH2} (red)",'${RED1st}' using 1:3  with lines lw 2  lc rgbcolor "${COLOR4}" title "${LABH2O} (red)"

quit
PLOT
echo "DONE!"
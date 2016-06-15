#!/bin/bash 

## Takes 2 arguments: $1: number of slides and $2: file name addition (FULL,
##                    RED_1st, ...
TMPDAT=movie.data
echo "set yrange [0.1:0.8]" > $TMPDAT
for((i=1;i<$1;i++)) do 
  SOL=Sol1D_$2_$i.gnu
  echo "plot \"$SOL\" using 1:2 w l lw 2, \"$SOL\" using 1:3 w l lw 2" >> $TMPDAT
  echo "pause 0.15" >> $TMPDAT
done  

cat $TMPDAT | gnuplot 

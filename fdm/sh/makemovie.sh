#!/bin/bash 

## Takes 2 arguments: $1: number of slides and $2: file name 
TMPDAT=movie.data
echo "set yrange [0.1:0.8]" > $TMPDAT
for((i=1;i<$1;i++)) do 
  SOL=$2$i.gnu
  echo "splot \"$SOL\" using 1:2:$3" >> $TMPDAT
  echo "pause 0.15" >> $TMPDAT
done  

cat $TMPDAT | gnuplot 

#!/bin/bash 

## Takes 2 arguments: $1: number of slides 
##                    $2: filename1 WITHOUT number
##                    $3: filename2 WITHOUT number   
##                        and suffix, e.g. 'data13.fig' then write 'data'
##                    $4: moviename (without suffix)

SUFFIX=avi
MOVIENAME=$4
TYPE=jpg
FRAMEPERSEC=2


SOL=
LOS=

echo "Removing existing '.${TYPE}' files..."
rm *.${TYPE}

for((i=1;i<$1;i++)) do 
## write into right ordering format, i.e. 0000, 0001, 0002, etc. needed for 
## mencoder
## check out "http://www.astrophysik.uni-kiel.de/index.php?option=com_content&view=article&id=80%3Agnuplot&Itemid=88&lang=de"
    if [ $i -lt 10 ]; then
        ext="000$i"
    elif [ $i -lt 100 ]; then
        ext="00$i"
    elif [ $i -lt 1000 ]; then
        ext="0$i"
    else
        ext="$i"
    fi
   SOL=$2$i.gnu
   LOS=$3$i.gnu
  gnuplot -e "set terminal jpeg; set mxtics 2; set mytics 2; set grid xtics ytics mxtics mytics; set title 'Timestep: ${ext}'; plot '$SOL' using 1:2 w l lw 2 title \"[H2] (full)\", '$SOL' using 1:3 w l lw 2 title \"[H2O] (full)\",'$LOS' using 1:2 w l lw 2 title \"[H2] (red)\", '$LOS' using 1:3 w l lw 2 title \"[H2O] (red)\"" > pic_$ext.${TYPE}
done  

#ffmpeg -i pic_%d.$TYPE ${MOVIENAME}.${SUFFIX}


## run mencoder from mplayer package:
mencoder mf://*${TYPE} -mf fps=${FRAMEPERSEC}:type=${TYPE} -ovc lavc -lavcopts vcodec=msmpeg4v2 -oac copy -o ${MOVIENAME}.${SUFFIX}
#!/bin/bash 

## Takes 2 arguments: $1: number of slides 
##                    $2: file name WITHOUT number and suffix, i.e.
##                        instead of 'data13.fig', write 'data'
##                    $3: moviename (wo suffix)
##                    $4: column to be plotted 

SUFFIX=avi
MOVIENAME=$3
TYPE=jpg
FRAMEPERSEC=5

rm *jpg   #erase all jpg - otherwise you'll spoil your results


if [ $4 -eq 1 ]; then
    echo -e "This is the x-coordinate. Plot from col 2 onwards only \nEXIT"
    exit 113
fi

SOL=
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
  gnuplot -e "set terminal jpeg; plot '$SOL' using 1:$4 w l lw 2" > pic_$ext.${TYPE}
done  

#ffmpeg -i pic_%d.$TYPE ${MOVIENAME}.${SUFFIX}


## run mencoder from mplayer package:
mencoder mf://*${TYPE} -mf fps=${FRAMEPERSEC}:type=${TYPE} -ovc lavc -lavcopts vcodec=msmpeg4v2 -oac copy -o ${MOVIENAME}.${SUFFIX}
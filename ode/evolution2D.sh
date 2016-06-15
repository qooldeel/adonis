#!/bin/bash 
## Creates movie of evolution on a rectangular 2D domain.  
## Conservation variables are stored columnwise, for each (x,y) value
##
## Takes 2 arguments: $1: number of slides 
##                    $2: file name WITHOUT number 
##                        and suffix, e.g. 'data13.fig' then write 'data'
##                    $3: moviename
##                    $4: column to be extracted. In 2D:
##                           1   2   3    4    5    6  || 7  ....  
##                           x   y   rho  v1   v2   T  || species
##                         

SUFFIX=avi
MOVIENAME=$3
TYPE=jpg
FRAMEPERSEC=5
COLUMN=$4

XVAL=1
YVAL=2

echo "Remove all '${TYPE}' files first..."
rm *.${TYPE}

if [ -z "$4" ]; then   # test if variable is set
    printf "Number of input arguments incomplete, pal!!\n"
    exit 1
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

## insert this to get contour lines
## set contour base; set cntrparam cubicspline; set cntrparam levels auto 20; unset clabel
    gnuplot -e "set terminal jpeg; set pm3d map; set palette rgbformulae 33,13,10; set title 'Timestep: ${ext}'; set contour base; set cntrparam levels auto 20; unset clabel; splot '$SOL' using ${XVAL}:${YVAL}:${COLUMN}; unset pm3d" > pic_$ext.${TYPE}
done  

#ffmpeg -i pic_%d.$TYPE ${MOVIENAME}.${SUFFIX}


## run mencoder from mplayer package:
mencoder mf://*${TYPE} -mf fps=${FRAMEPERSEC}:type=${TYPE} -ovc lavc -lavcopts vcodec=msmpeg4v2 -oac copy -o ${MOVIENAME}.${SUFFIX}
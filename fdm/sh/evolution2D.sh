#!/bin/bash 
## Creates movie of evolution on a rectangular 2D domain.  
## Conservation variables are stored columnwise, for each (x,y) value
##
## Takes 4 arguments: $1: number of slides 
##                    $2: file name WITHOUT number 
##                        and suffix, e.g. 'data13.fig' then write 'data'
##                    $3: moviename
##                    $4: column to be extracted. In 2D:
##                           1   2   3    4    5    6  || 7  ....  
##                           x   y   rho  v1   v2   T  || species
##                    $5: (optional) if you start from slide n != 0
## optionally, a fifth argument can be provided, especially when
## the slides are not counted from index 1.                        

SUFFIX=avi
MOVIENAME=$3
TYPE=jpg
COLUMN=$4

FRAMEPERSEC=5
SMOOTH=bspline
ORDER=10  #the bigger the number the smoother the contours
LEVEL=50  #20

XVAL=1  #x-axis values
YVAL=2  #y-axis values

# e.g. if you are not starting from slide 1 but, say, from slide 51
OPTIONALARG=

FROMSLIDE=0

# use a fifth argument by checking the number of arguments via $#
if [ $# -ge 4 ]; then
    echo "O.k. you are using more than 4 arguments right now. 5th arg = '$5'"
    OPTIONALARG=$5
    FROMSLIDE=${OPTIONALARG}
fi


if [ $4 -eq 1 ]; then
    echo -e "This is the x-coordinate. Plot from col 3 onwards only \nEXIT"
    exit 113
elif [ $4 -eq 2 ]; then 
    echo -e "This is the y-coordinate. Plot from col 3 onwards only \nEXIT"
    exit 113
fi

echo "Remove all '${TYPE}' files first..."
rm *.${TYPE}

if [ -z "$4" ]; then   # test if variable is set
    printf "Number of input arguments incomplete, pal!!\n"
    exit 1
fi

IDX=
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
    
    IDX=$i
    if [[ ${FROMSLIDE} -ne 0 ]]; then # mind the double brackets here
	IDX=$((IDX + $((FROMSLIDE - 1))))  #perform addition
	#echo " IDX = $IDX"
    fi

    SOL=$2${IDX}.gnu

## insert this to get contour lines
## set contour base; set cntrparam cubicspline; set cntrparam levels auto 20; unset clabel
    gnuplot -e "set terminal jpeg; set autoscale fix; set pm3d map interpolate 0,0; set palette rgbformulae 33,13,10; set title 'Timestep: ${ext}'; set contour base; set cntrparam order ${ORDER};set cntrparam ${SMOOTH}; set cntrparam levels auto ${LEVEL}; unset clabel; splot '$SOL' using ${XVAL}:${YVAL}:${COLUMN}; unset pm3d" > pic_$ext.${TYPE}
done  

#ffmpeg -i pic_%d.$TYPE ${MOVIENAME}.${SUFFIX}


## run mencoder from mplayer package:
mencoder mf://*${TYPE} -mf fps=${FRAMEPERSEC}:type=${TYPE} -ovc lavc -lavcopts vcodec=msmpeg4v2 -oac copy -o ${MOVIENAME}.${SUFFIX}

echo "Done! :)"

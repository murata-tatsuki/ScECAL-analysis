#！/bin/sh
path="$PWD"


datafile=20201216_20210308
#datafile=20201216_1924_cosmicRayTst

inputfile="/mnt2/ScECAL_CR/results/trackFit/"
#inputfile="/mnt2/ScECAL_CR/results/wPedestal/"
caliFile="/mnt2/ScECAL_CR/results/calibration/"

DIR_IN=${inputfile}${datafile}".root"
DIR_OUT1=${caliFile}

./main ${DIR_IN} ${DIR_OUT1}
echo "all over"

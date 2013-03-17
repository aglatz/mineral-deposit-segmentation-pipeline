#!/bin/bash

usage()
{
	echo
	echo "Standalone execution of the N4 bias field correction from Slicer3."
	echo
	echo "Usage: $0 <Slicer3 basedir> <N4 arguments>"
	echo
	exit 1
}

if [ $# -lt 2 ]
then
	usage;
fi

LIB_PATH1=${1}/lib/Slicer3/Plugins/;
LIB_PATH2=${1}/lib/InsightToolkit/;
shift
env LD_LIBRARY_PATH=${LIB_PATH1}:${LIB_PATH2} ${LIB_PATH1}/N4ITKBiasFieldCorrection $@ 

exit 0;

#!/bin/bash

usage()
{
	echo
	echo "Usage: $0 <fsl path> <volume to check>"
	echo
	echo "This script returns the number of non-zero "
	echo "voxels of the input volume."
	echo
	exit 1
}

if [ $# -ne 2 ]
then
	usage;
fi

VOL=`${1}fslstats ${2} -V 2> /dev/zero`;
if [ $? -eq 0 ] # check if input is ok
then
	echo $VOL | cut -d ' ' -f 1;
else
	echo 0;
fi

exit 0

#!/bin/bash

SLC_PFX=/usr/local/Slicer3/
FSL_PFX=/usr/local/fsl/4.1.9/

# Return environment variable name
VAR_NAME=$1
VAR_VAL=`set | grep -w "${VAR_NAME}=" | sed s/"${VAR_NAME}="//g`

if [ -z $VAR_VAL ]
then
	echo "";
	exit 1
else
	echo $VAR_VAL
	exit 0
fi

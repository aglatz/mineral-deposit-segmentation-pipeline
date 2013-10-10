#!/bin/bash

SDIR=$1 # Script directory
FSL_PFX=`${SDIR}/setup_environment.sh FSL_PFX` # Load variables
IN=$2 # Input filename
OUT=$3 # Target filename

# Simple conversion to get file format right
#${FSL_PFX}bin/fslmaths $IN $OUT
${SDIR}/exec_matlab.sh ${SDIR}/prepareinput.mt $IN $OUT

#!/bin/bash

SDIR=$1 # Script directory
FSL_PFX=`${SDIR}/setup_environment.sh FSL_PFX` # Load variables
IN=$2 # Input filename
OUT=$3 # Target filename

# Simple conversion to get file format right
${FSL_PFX}bin/fslmaths $IN $OUT

SFORMCODE=`${FSL_PFX}bin/fslorient -getsformcode $OUT`
QFORMCODE=`${FSL_PFX}bin/fslorient -getqformcode $OUT`
CMDARG=
if [ $SFORMCODE -eq 0 ]
then
	CMDARG="-setsformcode 2"
fi
if [ $QFORMCODE -eq 0 ]
then
	CMDARG="$CMDARG -setqformcode 2"
fi
if [ "$CMDARG" != "" ]
then
	# In case of ANALYZE files we have to set the
	# q and s form code. Otherwise N4 doesn't work
	# correctrly and flips the volumes left-right.
	${FSL_PFX}bin/fslorient $CMDARG $OUT
fi


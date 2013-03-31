#!/bin/bash
# This script invokes the MATLAB template script (first argument) after
# substituting the bash variables ARG_0, ARG_1, ... in the script
# with the remaining arguments.
# The script is executed in the directory where it resides (e.g.
# /my/script/hu.mt is executed in /my/script).
#
# Author: Andreas Glatz <a.glatz@sms.ed.ac.uk>
#

usage() 
{
	echo
	echo "Usage: $0 <MATLAB script> [ARG_0 ARG_1 ... ]"
	echo
	echo "       <MATLAB template script> MATLAB template script to execute."
	echo "       [ARG_0 ARG_1 ... ]       Arguments for the MATLAB script."
	echo
	echo "Note: The template script must have the ending .mt"
	exit 1
}

function canonize_path()
{
	local IN_PATH="${1}"
	cd -P -- ${IN_PATH} &> /dev/null && \
	echo $(pwd -P)
}


# MAIN ######################################################################
SCRIPT_EXT=".mt"
SCRIPT_PATH=$1
if [ ! -f $SCRIPT_PATH ]
then
	echo "Error: Could not read MATLAB script ${SCRIPT_PATH}!"
	usage
fi
SCRIPT_NAME=`basename $SCRIPT_PATH $SCRIPT_EXT`
SCRIPT_NAME_TMP=`basename $SCRIPT_PATH`
if [ "$SCRIPT_NAME" = "$SCRIPT_NAME_TMP" ]
then
	echo "Error: Matlab template script does end in ${SCRIPT_EXT}!"
	usage
fi
SCRIPT_DIR=$(canonize_path $(dirname $SCRIPT_PATH));

# Build sed command with all arguments
ARG_CNT=2
SC='+'
START_MSG="START_${SCRIPT_NAME}"
ERROR_MSG="ERROR_${SCRIPT_NAME}"
SEDEXP="s${SC}\${*SCRIPT_DIR}*${SC}${SCRIPT_DIR}${SC}g;"
SEDEXP="s${SC}\${*ERROR_MSG}*${SC}${ERROR_MSG}${SC}g;$SEDEXP"
SEDEXP="s${SC}\${*START_MSG}*${SC}${START_MSG}${SC}g;$SEDEXP"
while [ $ARG_CNT -le $# ]
do
	ARG_NAME="ARG_$((ARG_CNT-2))";
	eval ARG_VAL=$(echo \"\$${ARG_CNT}\");
	SEDEXP="s${SC}\${*${ARG_NAME}}*${SC}${ARG_VAL}${SC}g;$SEDEXP"
	ARG_CNT=$((ARG_CNT+1))
done

# Use mktemp - it fails if the same file already exists.
# Add two random numbers because previous versions of mktemp
# do not replace X's if they are not behind a '.' and at the end
COOKIE=$RANDOM
TMPSCRIPT_PATH=tmp${COOKIE}XXX_${SCRIPT_NAME}.m
TMPSCRIPT_PATH=`mktemp $TMPSCRIPT_PATH`
if [ $? -ne 0 ]
then
	echo "Error: Could not create $TMPSCRIPT_PATH!"
	usage
fi

# Substitude arguments and add to file
SEDEXP="sed -e '${SEDEXP}' $SCRIPT_PATH >> $TMPSCRIPT_PATH"
eval $SEDEXP
if [ $? -ne 0 ]
then
	echo "Error: Could not substitute arguments!"
	usage
fi

# Check if we have matlab
FOUND=0
matlab -e &> /dev/zero && FOUND=1
if [ $FOUND -eq 0 ]
then
	echo "Error: Command 'matlab' not found!"
	usage
fi

# Run script
ERROR=0
OUT_FILE=`echo $TMPSCRIPT_PATH | sed s/tmp/out/g`
ERR_FILE=`echo $TMPSCRIPT_PATH | sed s/tmp/err/g`
# echo $OUT_FILE $ERR_FILE
# Just collect output for standart error (2)
matlab -nodesktop -nosplash -r `basename ${TMPSCRIPT_PATH} .m` \
&> $ERR_FILE > $OUT_FILE
cat $ERR_FILE | grep ${ERROR_MSG} > /dev/zero && ERROR=1
cat $OUT_FILE | sed 1,/${START_MSG}/d
rm -f $ERR_FILE $OUT_FILE
if [ $ERROR -eq 1 ]
then
	echo "Error: Cannot excecute '${TMPSCRIPT_PATH}'!"
	rm -f ${IMG_OUT}
	usage
else
	rm -f ${TMPSCRIPT_PATH}
fi
exit 0
#!/bin/bash

get_volume()
{
	TMP=`fslstats $1 -V`
	TMP=`echo $TMP | cut -d ' ' -f 2`
	echo $TMP | cut -d '.' -f 1
}

get_sum()
{
	VOL=`get_volume $1`
	MEAN=`fslstats $1 -M`
	SUM=`echo "$VOL*$MEAN" | bc`
	echo $SUM | cut -d '.' -f 1
}

### MAIN ######################################################################
if [ $# -lt 2 ]
then
	echo "Usage: $0 <subject csv file>"
	exit 1
fi

FSLOUTPUTTYPE=NIFTI_GZ
NAME=${2}_bin_notreg_mni_mask
IN=${NAME}.nii.gz
OUT=./${NAME}_comb.nii.gz
SUBJECTS=$1

rm -rf $OUT
VOL_MNI=0
VOL_ORI=0
CNT=0;
while read LINE
do
	MASK=`echo $LINE | cut -d ',' -f 1`;
	MASK="../$MASK/$IN";
	echo "Adding $MASK ..."
	VOL1=`get_volume $MASK`
	if [ $VOL1 -gt 0 ]
	then
		CNT=$((CNT+1));
		fslmaths $MASK -div $VOL1 tmp_mask -odt float;
		if [ -f $OUT ];
		then
			fslmaths $OUT -add tmp_mask $OUT -odt float;
		else
			fslmaths tmp_mask $OUT -odt float;
		fi;
		VOL_MNI=$((VOL_MNI+VOL1))
		TMP="`dirname $MASK`/T2swHypo_mask"
		VOL=`get_volume $TMP`
		VOL_ORI=$((VOL_ORI+VOL))
		VOL2=`get_sum $OUT`
		echo "$VOL1 $VOL $VOL_MNI $VOL_ORI $VOL2" 
	fi
done < $SUBJECTS
fslmaths $OUT -div $((CNT*8)) $OUT -odt float;

exit 0

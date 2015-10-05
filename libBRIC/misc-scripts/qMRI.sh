#!/bin/bash
# This scrip converts DICOM to NIFTI, then fits mono-exponentials to voxel
# signals of multi-echo sequences.
#

me2r2()
{
	STR=$1
	OUT=$2
	SKIP=$3

	cd $OUTDIR
	TMP=`grep -rwn $STR\] * | cut -d '/' -f 1 | uniq | xargs`
	echo "Found $STR in the following directories: $TMP"
	for DIR in $TMP
	do
		ECHOES=`ls $DIR/Input/S_*.nii.gz | sed 's+.*_\([0-9][0-9]*\)\..*+\1+g' | sort -g | xargs`
		echo "All Echoes: $ECHOES"
		if [ $SKIP -gt 0 ]
		then
			NUM=0
			CNT=0
			RES=
			for E in $ECHOES
			do
				if [ $CNT -ge $SKIP ]
				then
					CNT=0
					if [ $NUM -le 0 ]
					then
						RES=$E
					else
						RES="$RES $E"
					fi
					NUM=$((NUM+1))
				else
					CNT=$((CNT+1))
				fi
			done
			ECHOES=$RES
			echo "Echoes after skipping $SKIP: $ECHOES"
		fi
		
		echo "convert_genii: $OUTDIR/$DIR"
		$TOOLSDIR/libBRIC/misc-scripts/exec_matlab.sh \
			$TOOLSDIR/libBRIC/qMRI/convert_genii.mt \
			$DIR/Input/S \
			`echo $ECHOES | wc -w` \
			`echo $ECHOES | cut -d ' ' -f 1,2` \
			0 0 \
			$DIR/S \
			1
		rm -f $DIR/S_*.img $DIR/S_mag_ln.nii.gz \
			  $DIR/S_rea.nii.gz $DIR/S_ima.nii.gz

		echo "me2R2: $OUTDIR/$DIR"
		$TOOLSDIR/libBRIC/misc-scripts/exec_matlab.sh \
			$TOOLSDIR/libBRIC/qMRI/me2R2.mt \
			$DIR/S_mag \
			`echo $ECHOES | wc -w` \
			`echo $ECHOES | cut -d ' ' -f 1,2` \
			100 \
			0 \
			$DIR/$OUT
	done
}

usage() 
{
	echo "usage: $0 <tools dir> <dcmtk dir> <dicom dir> <output dir>"
    exit 1
}

checkdir()
{
	local _resultvar=$1
	local _indir=$2

	cd $_indir &> /dev/zero
	if [ $? -ne 0 ]
	then
		echo "ERROR: $_resultval directory does not exist: $_indir. Exiting."
		exit 1
	fi
	eval $_resultvar="`pwd`"
	cd - &> /dev/zero
}

### MAIN ######################################################################
if [ $# -lt 4 ]
then
	usage;
fi

# /home/s1063233/mineral-deposit-segmentation-pipeline 
checkdir TOOLSDIR $1

# /home/s1063233/dcmtk/install
checkdir DCMTKTOOL $2
export PATH=$PATH:$DCMTKTOOL/bin
DCMDICTPATH=$DCMTKTOOL/share/dcmtk/dicom.dic:$DCMTKTOOL/share/dcmtk/private.dic
export DCMDICTPATH=$DCMDICTPATH

# 1 2 3 4 ...
checkdir INDIR $3

# output dir
checkdir OUTDIR $4

# Convert DICOMs
if true; then
echo "DCM2NII: $IN -> $OUT"
for DIR in `ls $INDIR | xargs`
do
	echo "DCM2NII: $INDIR/$DIR"
	mkdir -p $OUTDIR/$DIR/Input;
	$TOOLSDIR/libBRIC/misc-scripts/dcm2nii.pl \
		$INDIR/$DIR $OUTDIR/$DIR/Input &> $OUTDIR/$DIR/Input/dcm2nii.pl.log
done
fi

# Reconstruct R2s maps
echo ""
echo "------------"
read -p "Do you want to reconstruct the T2* map from MEGE? [y|n]" yn
if [ x$yn = xy -o x$yn = xY ]
then
	echo "me2r2 MEGE R2s 0"
	me2r2 MEGE R2s 0
fi

# Reconstruct R2 maps
echo ""
echo "------------"
read -p "Do you want to reconstruct the T2 map from MESE? [y|n]" yn
if [ x$yn = xy -o x$yn = xY ]
then
	echo "me2r2 2DMESE R2 0"
	me2r2 2DMESE R2 0 #1
fi

echo ""
echo "DONE."

exit 0;

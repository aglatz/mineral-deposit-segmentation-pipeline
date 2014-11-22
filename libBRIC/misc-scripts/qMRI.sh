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
	TMP=`grep -rwn $STR * | cut -d '/' -f 1 | uniq | xargs`
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
			10 \
			0 \
			$DIR/$OUT
	done
}

### MAIN ######################################################################
if [ $# -lt 2 ]
then
	echo "usage: $0 <indir> <outdir>"
	exit 1
fi

IN=$1
OUT=$2

# Setup the environment
cd /home/s1063233/mineral-deposit-segmentation-pipeline || exit 1;
TOOLSDIR=`pwd`;

cd /home/s1063233/dcmtk/install || exit 1;
PATH=$PATH:`pwd`/bin;
DCMTKTOOL=`pwd`
DCMDICTPATH=$DCMTKTOOL/share/dcmtk/dicom.dic:$DCMTKTOOL/share/dcmtk/private.dic
export DCMDICTPATH=$DCMDICTPATH

# Check the input
cd $IN || exit 1; export INDIR=`pwd`
cd $OUT || exit 1; export OUTDIR=`pwd`

# Convert DICOMs
echo "DCM2NII: $IN -> $OUT"
for DIR in `ls $INDIR | xargs`
do
	echo "DCM2NII: $INDIR/$DIR"
	mkdir -p $OUTDIR/$DIR/Input;
	$TOOLSDIR/libBRIC/misc-scripts/dcm2nii.pl \
		$INDIR/$DIR $OUTDIR/$DIR/Input &> $OUTDIR/$DIR/Input/dcm2nii.pl.log
done

# Reconstruct R2s maps
echo ""
echo "------------"
me2r2 3DMEGE R2s 0

# Reconstruct R2 maps
echo ""
echo "------------"
me2r2 2DMESE R2 1

echo ""
echo "DONE."

exit 0;

#!/bin/bash

if [ $# -lt 2 ]
then
	echo "Usage $0 <archive version> <archive.tar.gz> <basedir>"
	exit 1
fi

VER=$1
AR=$2
BDIR=$3

ARBDIR=`basename $AR .tar.gz`
SNO=`echo $ARBDIR | rev | cut -d '_' -f 1 | rev`
if [ $VER -gt 1 ]
then
	NUM=2
else
	NUM=3
fi
ARRDIRS=`tar tfz $AR | cut -d '/' -f $NUM | uniq | xargs`
echo $ARRDIRS

for ARRDIR in $ARRDIRS
do
	if [ $VER -gt 1 ]
	then
		ARDIR="$ARBDIR/$ARRDIR"
	else
		ARDIR="$ARBDIR/$SNO/$ARRDIR"
	fi
	echo -n "Working on: $ARDIR ..."
	tar -zxvf $AR $ARDIR
	DIR="$BDIR/$ARRDIR/Input"
	echo "$DIR ..."
	mkdir -p $DIR
	`dirname $0`/dcm2nii.pl $ARDIR $DIR
	rm -r $ARBDIR
done

exit 0

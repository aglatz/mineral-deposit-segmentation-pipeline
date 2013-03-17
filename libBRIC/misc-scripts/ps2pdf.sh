#!/bin/bash
# This script converts *.ps files into *.pdf files using the Linux
# script ps2pdf. This script also implements a workaround so that
# the images of the ps file are saved in the same way in the pdf
# (no change in color depth or resolution) and the
# colormap of the images are embedded in the final pdf.
#

INDIR=`dirname "$1"`
INFILE=`basename "$1" .ps`

if [ "$INFILE" = "" ]
then
	echo
	echo "Usage $0 <input.ps>"
	echo
	exit 1
fi

IN=${INDIR}/${INFILE}
OUT=${INDIR}/${INFILE}
echo "Info: Writing output to ${OUT}.pdf"

# Append line to ps file to tell ps2pdf leave
# the images as they are intended to be...
TMPFILE=${IN}_tmp
sed '/%%EndComments/ a\
% AGL\
/setdistillerparams where {pop}{userdict /setdistillerparams {pop} put}ifelse\
<< /AutoFilterColorImages false\
/AutoFilterGrayImages false\
/ColorImageFilter /FlateEncode\
/GrayImageFilter /FlateEncode\
/UseCIEColor true\
/ColorACSImageDict .printerACSImageDict >>\
setdistillerparams\
%AGL\
' ${IN}.ps > ${TMPFILE}.ps

# Convert ps to pdf
ps2pdf ${TMPFILE}.ps ${TMPFILE}.pdf
rm -f ${TMPFILE}.ps

# Rename outfile
mv ${TMPFILE}.pdf ${OUT}.pdf

exit 0

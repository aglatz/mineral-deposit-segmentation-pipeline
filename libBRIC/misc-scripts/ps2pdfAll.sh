#!/bin/sh
# This script searches *.ps files in all subdirectories and (i)converts them
# into pdf format using the script ps2pdf.sh, and (ii) compresses the original
# *.ps file with gzip to save disk space.
# This script has to reside in the same path as ps2pdf.sh.
#
PS2PDF_PATH=`dirname $0`; PS2PDF_PATH=`readlink -f $PS2PDF_PATH`; 

find -L . -name \*.ps \
	-execdir sh -c "cd `dirname {}`;FILE=`basename {}`;env PATH=$PATH:$PS2PDF_PATH ps2pdf.sh \$FILE; gzip -fv \$FILE" \;

#!/bin/bash
# Executes the pipline that is specified in the Makefile from the current directory.
# Author: Andreas Glatz <a.glatz@sms.ed.ac.uk>
#

usage()
{
	echo
	echo "Usage: $0 -s <node script dir> -N <#cpus> -k -v -r [all|<target1,target2,...>] -i <CSV file>"
	echo
	echo "Options:"
	echo "   -s <node script dir>           Directory with node scripts"
	echo "   -N <#cpus>                     Max. number of CPU cores available to distribute work"
	echo "   -k                             Keep intermediate files"
	echo "   -v                             Verbose checking of input file"
	echo "   -r [all|<target1,target2,...>] Remove specified targets from intermediate directories"
	echo "   -i                             Skip checking paths"
	echo "Operants:"
	echo "   <CSV file>                     CSV file containing paths for pipeline"
	echo
	exit 1
}

get_cannonical_path()
{
	cd $1 &>/dev/zero;
	if [ $? -eq 0 ];
	then
		pwd; 
		cd - &> /dev/zero;
	else
		# This will give an error further down
		# as the directory does not exist.
		echo "/$1";
	fi
}

strip_whitespace()
{
	echo "$1" | tr -cd '[:print:]' | sed 's/^ *//; s/ *$//'
}

# Getopts ######################################################################
echo; # Pretty output

SCRIPT_DIR=
MAX_JOBS=1
INP_FILE=
VERBOSE_FLAG=0
KEEP_FLAG=0
REMOVE=
SKIPCHECK_FLAG=0
while getopts s:f:N:r:kvi OPT
do
	case $OPT in
		s) SCRIPT_DIR=`get_cannonical_path $OPTARG`;;
		N) MAX_JOBS=$OPTARG;;
		k) KEEP_FLAG=1;;
		v) VERBOSE_FLAG=1;;
		r) REMOVE=$OPTARG;;
		i) SKIPCHECK_FLAG=1;;
		*) usage;;
	esac
done
shift $((OPTIND-1));
INP_FILE=$*;

if [ -z $INP_FILE ] || [ ! -f $INP_FILE ];
then
	echo "Error: Cannot read CSV file $INP_FILE";
	usage;
fi
if [ ! -z $SCRIPT_DIR ] && [ ! -d $SCRIPT_DIR ];
then
	echo "Error: argument to -s not a directory!"
	usage;
fi;
MAX_CPU=`cat /proc/cpuinfo | grep processor | tail -1 | cut -d ':' -f 2`
[ $MAX_JOBS -eq 1 ] &> /dev/zero;
if [ $? -eq 2 ]
then
	IS_NUM=0;
else
	IS_NUM=1;
fi;
if [ ! -z $MAX_JOBS ];
then
	if [ $IS_NUM -eq 1 ];
	then
		if [ $MAX_JOBS -lt 1 -o $MAX_JOBS -gt $((MAX_CPU+1)) ];
		then
			echo "Error: argument to -N is less than 1 or greater than $((MAX_CPU+1))!"
			usage;
		fi;
	else
		echo "Error: argument to -N is not a number!"
		usage;
	fi;
fi;

# MAIN starts here ############################################################
BASE_DIR=`pwd`
JOBS_RUNNING=0
COOKIE=$RANDOM
OUT_DIR_TEMP="`basename $0 .sh`_${COOKIE}XXX"
OUT_DIR=`mktemp -d $OUT_DIR_TEMP`
if [ $? -ne 0 ]
then
	echo "Error: Could not create output directory!"
	exit 1
fi
chmod 0755 ${OUT_DIR}
SUFFIX=".csv"
MAKE_FILE="${OUT_DIR}/Makefile"
LOG_FILE="${MAKE_FILE}.log"
DATA_FILE="${OUT_DIR}/all${SUFFIX}"
SANE_FILE="${OUT_DIR}/sane${SUFFIX}"
ERR_FILE="${OUT_DIR}/error${SUFFIX}"
STAT_FILE_NAME="stat${SUFFIX}"
STAT_FILE="${OUT_DIR}/${STAT_FILE_NAME}"
INPM_FILE="Input${SUFFIX}"
OUTM_FILE="Output${SUFFIX}"

# Create static part of the Makefile
EXTRA_ARGS="KEEP_FLAG=$KEEP_FLAG";
if [ ! -z $SCRIPT_DIR ];
then
	EXTRA_ARGS="SCRIPT_DIR=${SCRIPT_DIR}/ $EXTRA_ARGS";
fi
cat << EOF > $MAKE_FILE
.DEFAULT:
	+ \$(MAKE)	-C \$@ -f ${BASE_DIR}/Makefile \\
				CUR_DIR=\$@ PATH=\$(PATH):${BASE_DIR} \\
				$EXTRA_ARGS

EOF

# Create subject directories and add subject targets
ALL_DIR=
LINE_NO=1
echo "Info: Checking paths in file $INP_FILE..."
while read LINE
do
	CUR_DIR=`echo $LINE | cut -d ',' -f 1`
	
	# Remove all leading '/' (we don't allow absolute paths)
	# and record result in the data file
	CUR_DIR=`echo $CUR_DIR | sed -e 's@^/*\(.*\)@\1@g'`
	echo $CUR_DIR >> $DATA_FILE

	# Create dir if it doesn't exists
	mkdir -p ${CUR_DIR}
	if [ $? -ne 0 ]
	then
		echo -n "Error in line $LINE_NO of file $INP_FILE: ";
		echo "First column does not contain a valid subject path!";
		exit 1;
	fi
	chmod 0755 ${CUR_DIR}

	if [ "$REMOVE" != "" ]
	then
		if [ "$REMOVE" = "all" ]
		then
			if [ $VERBOSE_FLAG -eq 1 ]
			then
				echo "Info: Removing all targets from ${CUR_DIR}..."
			fi
			rm -rf ${CUR_DIR}/*
		else
			for TGT in ${REMOVE//,/ }
			do
				if [ $VERBOSE_FLAG -eq 1 ]
				then
					echo "Info: Removing ${CUR_DIR}/${TGT}.* ..."
				fi
				rm -f ${CUR_DIR}/${TGT}.*
			done
			# Remove missing targets files
			make -C ${CUR_DIR} -f ${BASE_DIR}/Makefile \
				$EXTRA_ARGS \
				clean-mtg &> $LOG_FILE
		fi
	else
		# Remove missing targets files
		make -C ${CUR_DIR} -f ${BASE_DIR}/Makefile \
			$EXTRA_ARGS \
			clean-mtg &> $LOG_FILE
	fi
	
	# Add current directory to makefile as a target
	ALL_DIR="$ALL_DIR \"${CUR_DIR}\""

	# Check all column entries with input data
	# if they are valid and add them to Input.csv
	rm -f ${CUR_DIR}/${INPM_FILE}; 
	COL_NO=2; # Current column number
	NO_COLS=8; # Total number of colums
	NO_REQ_COLS=3; # Number of required columns
	while [ $COL_NO -le $NO_COLS ]
	do
		EMPTY=0;
		PROBLEM=0;
		# Get right column
		STR=`echo $LINE | cut -d ',' -f $COL_NO`;
		# Remove non print chars and white space at the start/end
		STR=`strip_whitespace "$STR"`;
		if [ $SKIPCHECK_FLAG -eq 0 ]
		then
			if [ "$STR" != "" ]
			then
				RET=`sh -c "rsync -az --list-only $STR 2>/dev/zero"`
				if [ "$RET" = "" ]
				then
					PROBLEM=1;
				fi
			else
				EMPTY=1;
			fi
		fi
		if [ \( $EMPTY -eq 1 -o $PROBLEM -eq 1 \) -a $COL_NO -le $NO_REQ_COLS ]
		then
			echo -n "Error: Required path in column:$COL_NO, "
			echo "line:$LINE_NO is invalid! Exiting.";
			exit 1;
		fi
		if [ $EMPTY -eq 1 -a $VERBOSE_FLAG -eq 1 ]
		then
			echo -n "Info: Optional path in column:$COL_NO, "
			echo "line:$LINE_NO is missing!";
		fi
		if [ $PROBLEM -eq 1 -a $VERBOSE_FLAG -eq 1 ]
		then
			echo -n "Info: Optional path in column:$COL_NO, "
			echo "line:$LINE_NO is invalid and will be ignored!";
		fi
		echo -n "$STR," >> ${CUR_DIR}/${INPM_FILE}
		COL_NO=$((COL_NO+1));
	done
	echo "" >> ${CUR_DIR}/${INPM_FILE};

	# Check if there is an output directory
	rm -f ${CUR_DIR}/${OUTM_FILE};
	COL_NO=$((NO_COLS+1));
	STR=`echo $LINE | cut -d ',' -f $COL_NO`;
	# Remove non print chars and white space at the start/end
	STR=`strip_whitespace "$STR"`;
	if [ "$STR" != "" ]
	then
		if [ ! -d "$STR" ]
		then
			if [ $VERBOSE_FLAG -eq 1 ]
			then
				echo -n "Info: Creating output path $STR from "
				echo "column: $COL_NO, line:$LINE_NO.";
			fi
			mkdir -p "$STR"
		fi
		STR=`get_cannonical_path "$STR"`
		if [ "$STR" = "" ]
		then
			echo -n "Error: Could not create $STR from "
			echo "column: $COL_NO, line:$LINE_NO! Exiting.";
			exit 1;
		fi
	else
		if [ $VERBOSE_FLAG -eq 1 ]
		then
			echo -n "Info: Optional output path in column: $COL_NO, "
			echo "line:$LINE_NO is missing!"; 
		fi
	fi
	echo "$STR" >> ${CUR_DIR}/${OUTM_FILE}
	
	LINE_NO=$((LINE_NO+1));
done < $INP_FILE
echo "Info: Checked $((LINE_NO-1)) lines in file $INP_FILE."
	
# Add 'all' target and finalize makefile
cat << EOF >> $MAKE_FILE
.PHONY: all
all: ${ALL_DIR}
# END OF MAKEFILE
EOF

# Start processing ...
echo -n "Info: Starting processing on $MAX_JOBS CPU core(s) "
echo -n "(logging to $LOG_FILE)... "
make -k -f $MAKE_FILE -j $MAX_JOBS &> $LOG_FILE
echo "done."

# Sequentially gather stats of generated images
# and missing targets to avoid cluttered output.
FIRST_SUBJECT=1;
for CUR_DIR in `cat $DATA_FILE | xargs`
do
	# We explicitly change dir so that
	# make doesn't output that it entered/left
	# the directory and because print_image_stats()
	# needs to be executed in the subject dir.
	cd ${BASE_DIR}/${CUR_DIR}

	# Gather info about missing targets
	make -f ${BASE_DIR}/Makefile \
		$EXTRA_ARGS \
		print-mtg >> ${BASE_DIR}/${ERR_FILE}.tmp

	# Gather image stats
	if [ -f ${STAT_FILE_NAME} ]
	then
		# If it's the first subject then we also copy
		# the header.. otherwise just the data.
		if [ $FIRST_SUBJECT -eq 1 ]
		then
			cat ${STAT_FILE_NAME} >> ${BASE_DIR}/$STAT_FILE;
			FIRST_SUBJECT=0;
		else
			sed -e '1d' ${STAT_FILE_NAME} >> ${BASE_DIR}/$STAT_FILE;
		fi
	fi

	# Done ...
	cd ${BASE_DIR}
done

# Create list of targets with and without errors
if [ -f ${ERR_FILE}.tmp ]
then
	# Print some stats
	echo ----------------------------------------------------------------------
	echo Missing files:
	cat ${ERR_FILE}.tmp
	echo ----------------------------------------------------------------------

	# Compile the sorted file path list of complete directories
	cat ${ERR_FILE}.tmp \
	| sed -n -e 's/\(.*\)\/.*$/\1/p' | sort | uniq \
	> ${ERR_FILE}

	cat ${DATA_FILE} \
	| sort | uniq \
	> ${DATA_FILE}.tmp
	mv ${DATA_FILE}.tmp ${DATA_FILE}

	comm -23 ${DATA_FILE} ${ERR_FILE} \
	> ${SANE_FILE}
fi

# Done ...
exit 0

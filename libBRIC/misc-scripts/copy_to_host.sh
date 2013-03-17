#!/bin/bash

THOST_DEF="brain@dcn060062.dcn.ed.ac.uk"
TPORT_DEF=2222

usage()
{
	echo
	echo "Usage: $0 [-s <remote ssh host>[:<ssh port>]] <remote path> <items>"
	echo
	echo "Additionally, this application expects either a list of directories"
	echo "or files at stdin. If the application just finds files at stdin"
	echo "then the <items> argument can be obmitted."
	echo
	echo "Options:"
	echo "   -s             Remote ssh host and port, if necessary."
	echo
	echo "Operants:"
	echo "   <remote path>  The remote path where to copy the items."
	echo "   <items>        The items which should be copied from the local"
    echo "                  path to the remote path (could be redundant... "
	echo "                  see explanation above)."
	echo
	exit 1
}

check_dir()
{
	# Check if path is ok
	DIR_OK=y
	cd $1 || DIR_OK=n;
	if [ $DIR_OK == n ]
	then
		echo "Error: Cannot cd to $1!"
		usage;
	fi
	cd -
}


# Getopts #####################################################################
echo; # Pretty output

while getopts s:f:N:k OPT
do
	case $OPT in
		s) THOST=`echo $OPTARG | cut -d ':' -f 1`; 
		   TPORT=`echo $OPTARG | cut -d ':' -f 2`;
		   if [ "$THOST" != "" -a "$TPORT" != "" -a "$THOST" == "$TPORT" ]
		   then
				TPORT=
		   fi;;
		*) usage;;
	esac
done
shift $((OPTIND-1));

# Host
if [ -z $THOST ]
then
	echo "Info: Using default host: $THOST_DEF"
	THOST=$THOST_DEF
fi

# Port
if [ -z $TPORT ]
then
	echo "Info: Using default port: $TPORT_DEF"
	TPORT=$TPORT_DEF
fi

# Remote path
if [ $# -lt 1 ] # We need at least one operant (remote path).
then
	echo "Error: Too less operants!"
	echo
	usage;
fi
TPATH=$1;
shift;

# Items to copy in local directory
TOCOPY1=$*;

echo "Remote host, port and path: ----------"
echo "$THOST:$TPORT:$TPATH"
echo "Items to copy: -----------------------"
echo "$TOCOPY1"
echo "--------------------------------------"
echo 

CURPATH=`pwd`

# Let's do something...
for ITEM in `xargs`
do
	TOCOPY="" # Reset what we have to copy
	cd $CURPATH # Reset current path

	if [ -f $ITEM ]
	then
		if [ "$CURPATH" = "`dirname $ITEM`" ]
		then
			RELDIR="."
		else
			# Remove curpath portion (if $ITEM contains an absolute path)
			RELDIR=`dirname $ITEM | sed -e 's+'$CURPATH\/'++'`
		fi
		check_dir $RELDIR

		TOCOPY=$ITEM
	else
		if [ -d $ITEM ]
		then
			# Remove curpath portion (if $ITEM contains an absolute path)
			RELDIR=`echo $ITEM | sed -e 's+'$CURPATH\/'++'`
			check_dir $RELDIR

			if [ "$TOCOPY1" = "" ]
			then
				echo "Error: Specify items to copy!"
				usage;
			fi
			for TMP in $TOCOPY1
			do
				if [ -f "$RELDIR/$TMP" ]
				then
					TOCOPY="$TOCOPY $RELDIR/$TMP"
				fi
			done
			if [ "$TOCOPY" = "" ]
			then
				echo "Error: Cannot access any item to copy!"
				usage;
			fi
		else
			echo "Error: Unsupported file type $ITEM!"
			usage;
		fi
	fi

	# Create target dir
	TPATHFULL="${TPATH}/${RELDIR}"
	echo "ssh -p$TPORT $THOST mkdir -p $TPATHFULL"
	ssh -p$TPORT $THOST "mkdir -p $TPATHFULL" || exit 1

	# Copy
	echo "scp -P$TPORT $TOCOPY $THOST:$TPATHFULL" 
	scp -P$TPORT $TOCOPY $THOST:$TPATHFULL || exit 1
done

exit 0

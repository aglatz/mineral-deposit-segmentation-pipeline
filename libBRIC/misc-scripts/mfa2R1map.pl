#!/usr/bin/perl

use strict;
use warnings;

die("\nUsage: mfa2R1map.pl <despot> <noisescale> <reg_flag> <dir1> <dir2> ... <dirN>\n") if ($#ARGV < 4);

printf("Parameters: @ARGV\n");

my $despot = shift;
my $noisescale = shift;
my $do_reg = shift;
my @arglist;
my $argcnt = 0;
my $dir;
my $niinames = " ";
for $dir (@ARGV) {
	unless ( opendir(DH, "$dir/Input") ) {
		printf("Warning: Cannot read $dir! Wrong structure? Skipping...\n");
		next;
	}
	my @files = readdir(DH);
	closedir(DH);

	my %arg = (
		"idx" => $argcnt,
		"txtname" => "",
		"niiname" => "",
		"imgname" => "",
		"fa" => "",
		"tr" => "",
		"ti" => "",
		"nslices" => "",
		"txgain" => "",
		"arxgain" => "",
		"drxgain" => ""
	);

	my $file;
	my $cnt = 0;
	for $file (@files) {
		if ($file =~ /.txt/) {
			$arg{"txtname"} = "$dir/Input/$file";
			$cnt = $cnt + 1;
		}
		if ($file =~ /.nii/) {
			system("FSLOUTPUTTYPE=NIFTI_GZ; fslsplit $dir/Input/$file in_$argcnt" . "_ -t; " .
				   "fslcpgeom $dir/Input/$file in_$argcnt" . "_0000 -d");
			$arg{"niiname"} = "in_$argcnt" . "_0000";
			$cnt = $cnt + 1;
		}
	}

	unless ( $cnt == 2 ) {
		printf("Warning: Too many files in $dir/Input!\n");
		next;
	}

	unless ( open(FH, $arg{"txtname"}) ) {
		printf("Warning: Cannot open txt file" . $arg{"txtname"} . "! Skipping...\n");
		next;
	}
	$cnt = 0;
	my $slicethick;
	my $sliceposstart;
	my $sliceposend;
	while (my $line = <FH>) {
		if ( $line =~ /^\(0018,1314\) DS \[(.*)\]/ ) { # FlipAngle
			$arg{"fa"} = $1;
			$cnt = $cnt + 1;
		}
		if ( $line =~ /^\(0018,0080\) DS \[(.*)\]/ ) { # RepetitionTime
			$arg{"tr"} = $1;
			$cnt = $cnt + 1;
		}
		if ( $line =~ /^\(0018,0082\) DS \[(.*)\]/ ) { # InversionTime
			$arg{"ti"} = $1;
			$cnt = $cnt + 1;
		}
		if ( $line =~ /^\(0018,0050\) DS \[(.*)\]/ ) { # SliceThickness
			$slicethick = $1;
			$cnt = $cnt + 1;
		}
		if ( $line =~ /^\(0019,1019\) DS \[(.*)\]/ ) { # FirstScanLocation
			$sliceposstart = $1;
			$cnt = $cnt + 1;
		}
		if ( $line =~ /^\(0019,101b\) DS \[(.*)\]/ ) { # LastScanLocation
			$sliceposend = $1;
			$cnt = $cnt + 1;
		}
		if ( $line =~ /^\(0019,1094\) SS (\d*) / ) { # TransmitGain
			$arg{"txgain"} = $1;
			$cnt = $cnt + 1;
		}
		if ( $line =~ /^\(0019,1095\) SS (\d*) / ) { # AnalogReceiverGain
			$arg{"arxgain"} = $1;
			$cnt = $cnt + 1;
		}
		if ( $line =~ /^\(0019,1096\) SS (\d*) / ) { # DigitalReceiverGain
			$arg{"drxgain"} = $1;
			$cnt = $cnt + 1;
		}
	}
	close(FH);
	unless ( $cnt == 9 ) {
		printf("$cnt Warning: Not all fields valid in file $dir/Input" . $arg{"txtname"} . "!\n");
		next;
	}
	$arg{"nslices"} = abs($sliceposend - $sliceposstart) / $slicethick + 1;

	# same niinames for registration code below
	$niinames = $niinames . " " . $arg{'niiname'};

	$arglist[$argcnt] = \%arg;
	$argcnt = $argcnt + 1;
}

system("FSLOUTPUTTYPE=NIFTI_GZ; fslmerge -t S $niinames;");
if ( $do_reg ) {
	system("FSLOUTPUTTYPE=NIFTI_GZ; mv S.nii.gz S_tmp.nii.gz;" . 
		   "mcflirt -in S_tmp -out S -refvol 0 -plots -report");
}
system("FSLOUTPUTTYPE=NIFTI_GZ; fslsplit S S_ -t;");
$argcnt = 0;
for my $argref (@arglist) {
	$argref->{"imgname"} = "S_$argref->{'ti'}_$argref->{'tr'}_$argref->{'fa'}";	
	system(sprintf("FSLOUTPUTTYPE=NIFTI_PAIR; fslmaths S_%04d $argref->{'imgname'} -odt float", $argcnt));
	$argcnt = $argcnt + 1;
}
system("rm S_*.nii.gz");

my $arg1_name = "";
my $arg1_fa = "";
my $arg1cnt = 0;
my $arg2_name = "";
my $arg2_ti = "";
my $arg2cnt = 0;
my $tr = -1;
my $irtr = -1;
my $irfa = -1;
my $nslices = -1;
my $txgain = -1;
my $arxgain = -1;
my $drxgain = -1;
for my $argref (@arglist) {
	if ( $nslices < 0 ) {
		$nslices = $argref->{'nslices'};
	} else {
		if ( $argref->{'nslices'} != $nslices ) {
			printf("Warning: Nslices of $argref->{'niiname'} is different!\n");
		}
	}
	if ( $txgain < 0 ) {
		$txgain = $argref->{'txgain'};
	} else {
		if ( $argref->{'txgain'} != $txgain ) {
			printf("Warning: txgain of $argref->{'niiname'} is different!\n");
		}
	}
	if ( $arxgain < 0 ) {
		$arxgain = $argref->{'arxgain'};
	} else {
		if ( $argref->{'arxgain'} != $arxgain ) {
			printf("Warning: arxgain of $argref->{'niiname'} is different!\n");
		}
	}
	if ( $drxgain < 0 ) {
		$drxgain = $argref->{'drxgain'};
	} else {
		if ( $argref->{'drxgain'} != $drxgain ) {
			printf("Warning: drxgain of $argref->{'niiname'} is different!\n");
		}
	}

	if ( !$argref->{'ti'} ) {
		if ( $tr < 0 ) {
			# Set tr
			$tr = $argref->{'tr'};
		} else {
			# Check tr
			if ( $argref->{'tr'} != $tr ) {
				printf("Warning: TR of $argref->{'niiname'} is different!\n");
			}
		}
		$arg1_name = $arg1_name . " " . $argref->{'imgname'};
		$arg1_fa = $arg1_fa . " " . $argref->{'fa'};
		$arg1cnt = $arg1cnt + 1;
	} else {
		if ( $irtr < 0 ) {
			# Set tr
			$irtr = $argref->{'tr'};
		} else {
			# Check tr
			if ( $argref->{'tr'} != $irtr ) {
				printf("Warning: IRTR of $argref->{'niiname'} is different!\n");
			}
		}
		if ( $irfa < 0 ) {
			# Set fa
			$irfa = $argref->{'fa'};
		} else {
			# Check fa
			if ( $argref->{'fa'} != $irfa ) {
				printf("Warning: IRFA of $argref->{'niiname'} is different!\n");
			}
		}
		$arg2_name = $arg2_name . " " . $argref->{'imgname'};
		$arg2_ti = $arg2_ti . " " . $argref->{'ti'};
		$arg2cnt = $arg2cnt + 1;
	}
}

my $cmd;
my @output;
if ( !$arg2cnt ) {
	$cmd = "$despot 1 $arg1cnt $tr $arg1_name $arg1_fa ./ $noisescale";
	@output = ("DESPOT1_T1Map", "DESPOT1_MoMap");
} else {
	$nslices = $nslices + 2; # another 2 are added by despot
	# Use field strength "30" because it shortens the TI times to the actual
	# values... if not R1 values are too low.
	$cmd = "$despot 2 $arg1cnt $tr $arg1_name $arg1_fa $arg2cnt $arg2_name " .
		   "$arg2_ti $irtr $irfa $nslices 30 1 ./ $noisescale 1";
	@output = ("DESPOT1HIFI_B1Map", "DESPOT1HIFI_MoMap", "DESPOT1HIFI_T1Map");
}
printf("Calling: $cmd\n");
system($cmd);

for my $file (@output) {
	system("rm -f $file.nii.gz; fslchfiletype NIFTI_GZ $file");
}

system("rm -f *.hdr *.img");

for my $file (@output) {
	system("FSLOUTPUTTYPE=NIFTI_GZ; fslcpgeom $arglist[0]->{'niiname'} $file -d");
}

system("rm -f in_*.nii.gz");

exit 0;

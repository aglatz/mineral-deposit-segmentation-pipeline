#!/usr/bin/perl

use strict;
use warnings;

die("\nUsage: getacqparams.pl <dir1> ... <dirN>\n") if ($#ARGV < 0);

#printf("Parameters: @ARGV\n");

my $dir;
my $argcnt = 0;
printf("Dir,SliceThickness,RepetitionTime,EchoTime,InversionTime,NumberOfAverages," .
	   "PixelBandwidth,Mx,My,FlipAngle,FOV,NumberOfEchos,TransmitGain,AnalogReceiverGain," .
	   "DigitalReceiverGain,PulseSequenceName,Px,Py\n");
for $dir (@ARGV) {
	unless ( opendir(DH, "$dir/Input") ) {
		printf("Warning: Cannot read $dir! Wrong structure? Skipping...\n");
		next;
	}
	my @files = readdir(DH);
	closedir(DH);

	my $file;
	my $txtname;
	my $niiname;
	my $cnt = 0;
	for $file (@files) {
		if ($file =~ /.txt$/) {
			$txtname = "$dir/Input/$file";
			$cnt = $cnt + 1;
		}
	}

	unless ( $cnt == 1 ) {
		printf("Warning: Too many txt files in $dir/Input!\n");
		next;
	}

	unless ( open(FH, $txtname) ) {
		printf("Warning: Cannot open txt file" . $txtname . "! Skipping...\n");
		next;
	}
	$cnt = 0;
	printf("$dir,");
	while (my $line = <FH>) {
		if ( $line =~ /^\(0018,1314\) DS \[(.*)\]/ ) { # FlipAngle
			printf("$1,");
			$cnt = $cnt + 1;
		}
		if ( $line =~ /^\(0019,109c\) LO \[(.*)\]/ ) { # PulseSequenceName
			printf("$1,");
			$cnt = $cnt + 1;
		}
		if ( $line =~ /^\(0018,0080\) DS \[(.*)\]/ ) { # RepetitionTime
			printf("$1,");
			$cnt = $cnt + 1;
		}
		if ( $line =~ /^\(0018,0082\) DS \[(.*)\]/ ) { # InversionTime
			printf("$1,");
			$cnt = $cnt + 1;
		}
		if ( $line =~ /^\(0019,1094\) SS (\d*) / ) { # TransmitGain
			printf("$1,");
			$cnt = $cnt + 1;
		}
		if ( $line =~ /^\(0019,1095\) SS (\d*) / ) { # AnalogReceiverGain
			printf("$1,");
			$cnt = $cnt + 1;
		}
		if ( $line =~ /^\(0019,1096\) SS (\d*) / ) { # DigitalReceiverGain
			printf("$1,");
			$cnt = $cnt + 1;
		}
		if ( $line =~ /^\(0018,0081\) DS \[(.*)\] / ) { # EchoTime
			printf("$1,");
			$cnt = $cnt + 1;
		}
		if ( $line =~ /^\(0018,0083\) DS \[(.*)\] / ) { # NumberOfAverages
			printf("$1,");
			$cnt = $cnt + 1;
		}
		if ( $line =~ /^\(0019,107e\) SS (\d*) / ) { # NumberOfEchos
			printf("$1,");
			$cnt = $cnt + 1;
		}
		if ( $line =~ /^\(0018,0095\) DS \[(.*)\] / ) { # PixelBandwidth
			printf("$1,");
			$cnt = $cnt + 1;
		}
		if ( $line =~ /^\(0018,0050\) DS \[(.*)\] / ) { # SliceThickness
			printf("$1,");
			$cnt = $cnt + 1;
		}
		if ( $line =~ /^\(0018,1310\) US 0\\(\d*)\\(\d*)\\0 / ) { # AcquisitionMatrix
			printf("$1,$2,");
			$cnt = $cnt + 1;
		}
		if ( $line =~ /^\(0028,0030\) DS \[(.*)\\(.*)\] / ) { # PixelSpacing
			printf("$1,$2,");
			$cnt = $cnt + 1;
		}
		if ( $line =~ /^\(0019,101e\) DS \[(.*)\] / ) { # DisplayFieldOfView
			printf("$1,");
			$cnt = $cnt + 1;
		}
	}
	close(FH);
	printf("\n");
	unless ( $cnt == 15 ) {
		printf("$cnt Warning: Not all fields valid in file $dir/Input" . $txtname . "!\n");
		next;
	}
	$argcnt = $argcnt + 1;
}

exit 0;

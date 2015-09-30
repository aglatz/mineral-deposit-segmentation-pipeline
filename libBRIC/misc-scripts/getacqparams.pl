#!/usr/bin/perl

use strict;
use warnings;

die("\nUsage: getacqparams.pl <dir1> ... <dirN>\n") if ($#ARGV < 0);

#printf("Parameters: @ARGV\n");

# Print HTML header and header row
print("<!DOCTYPE HTML>\n");
print("<html>\n");
print("<body>\n");
print("<table style=\"width:100%\">\n");
println_tr();
my @hdr = (	"Dir",
			"FlipAngle",
			"PulseSequenceName",
			"RepetitionTime",
			"InversionTime",
			"TransmitGain",
			"AnalogReceiverGain",
			"DigitalReceiverGain",
			"EchoTime",
			"NumberOfAverages",
			"NumberOfEchos",
			"PixelBandwidth",
			"SliceThickness" ,
			"MatrixSize",
			"PixelSize",
			"FOV"
		  );
foreach my $col (@hdr) {
	println_td_ntd($col);
}
println_ntr();

for my $dir (@ARGV) {
	unless ( opendir(DH, "$dir/Input") ) {
		warn("Warning: Cannot read $dir! Wrong structure? Skipping...\n");
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
		warn("Warning: Too many txt files in $dir/Input!\n");
		next;
	}

	unless ( open(FH, $txtname) ) {
		warn("Warning: Cannot open txt file" . $txtname . "! Skipping...\n");
		next;
	}

	my @arr = ();
	while (my $line = <FH>) {
		if ( $line =~ /^\(0018,1314\) DS \[(.*)\]/ ) { # FlipAngle
			$arr[0] = $1;
		}
		if ( $line =~ /^\(0019,109c\) LO \[(.*)\]/ ) { # PulseSequenceName
			$arr[1] = $1;
		}
		if ( $line =~ /^\(0018,0080\) DS \[(.*)\]/ ) { # RepetitionTime
			$arr[2] = $1;
		}
		if ( $line =~ /^\(0018,0082\) DS \[(.*)\]/ ) { # InversionTime
			$arr[3] = $1;
		}
		if ( $line =~ /^\(0019,1094\) SS (\d*) / ) { # TransmitGain
			$arr[4] = $1;
		}
		if ( $line =~ /^\(0019,1095\) SS (\d*) / ) { # AnalogReceiverGain
			$arr[5] = $1;
		}
		if ( $line =~ /^\(0019,1096\) SS (\d*) / ) { # DigitalReceiverGain
			$arr[6] = $1;
		}
		if ( $line =~ /^\(0018,0081\) DS \[(.*)\] / ) { # EchoTime
			$arr[7] = $1;
		}
		if ( $line =~ /^\(0018,0083\) DS \[(.*)\] / ) { # NumberOfAverages
			$arr[8] = $1;
		}
		if ( $line =~ /^\(0019,107e\) SS (\d*) / ) { # NumberOfEchos
			$arr[9] = $1;
		}
		if ( $line =~ /^\(0018,0095\) DS \[(.*)\] / ) { # PixelBandwidth
			$arr[10] = $1;
		}
		if ( $line =~ /^\(0018,0050\) DS \[(.*)\] / ) { # SliceThickness
			$arr[11] = $1;
		}
		if ( $line =~ /^\(0018,1310\) US 0\\(\d*)\\(\d*)\\0 / ) { # AcquisitionMatrix
			$arr[12] = "$1, $2";
		}
		if ( $line =~ /^\(0028,0030\) DS \[(.*)\\(.*)\] / ) { # PixelSpacing
			$arr[13] = "$1, $2";
		}
		if ( $line =~ /^\(0019,101e\) DS \[(.*)\] / ) { # DisplayFieldOfView
			$arr[14] = $1;
		}
	}
	close(FH);

	# Print the current row
	println_tr();
	println_td_ntd($dir);
	for (my $i=0; $i<=$#arr; $i++) {
		if ($arr[$i]) {
			println_td_ntd($arr[$i]);
		} else {
			println_td_ntd_nan(1);
		}
	}
	println_ntr();
}

# Print end of HTML
print("</table>\n");
print("</body>\n");
print("</html>\n");

exit 0;

#--- Subs ---------------------------------------------------------------------
sub println_tr
{
	print("<tr>\n");
}

sub println_ntr
{
	print("</tr>\n");
}

sub println_td_ntd
{
	my $str = shift;
	print("\t<td>$str</td>\n");
}

sub println_td_ntd_nan
{
	my $i = 0;
	my $cnt = shift;
	
	while ($i < $cnt) {
		println_td_ntd("N/A");
		$i = $i + 1;
	}
}

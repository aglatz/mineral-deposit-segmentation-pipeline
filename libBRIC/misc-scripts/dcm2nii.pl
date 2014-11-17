#!/usr/bin/perl
# Requires: dcmtk and fsl utils, as well as dcm2nii

use strict;
use warnings;
use POSIX;

### Subroutines ###############################################################

sub uniq
{
	my %seen;
	return grep {!$seen{$_}++} @_;
}

sub generate_nifti
{
	my $type = shift;
	my $echotime = shift;
	my $indir = shift;
	my $outdir = shift;
	my $file;
	my @files;
	for $file (@_) {
		push @files, join('/', $indir, $file);
	}

	# Sanity checking
	my @echotimes;
	my @types;
	my $cmd = "dcmdump --search 0043,102f --search EchoTime @files |";
	open(ES, $cmd) || die("Error: checking dicom files!");
	while(my $line = <ES>) {
		if( $line =~ /.*\[(.*)\].*EchoTime/ ) {
			push @echotimes, $1;
		}
		if( $line =~ /.*SS (\d).*RawDataType/ ) {
			push @types, $1;
		}
	}
	close(ES);
	@echotimes = uniq @echotimes;
	@types = uniq @types;
	#print "@echotimes\n";
	#print "$#echotimes $echotimes[0] $echotime\n";
	#print "@types\n";
	#print "$#types $types[0] $type\n";
	if ($#echotimes != 0 || 
	    $echotimes[0]*1000 != $echotime ||
		$#types != 0) {
		die("Error: Inconsitent dicom files!");
	}

	# Generate shell script
	my $script_name = $outdir . sprintf("/dcm2nii_%d_%d.sh", $type, $echotime);
	printf("Info: Generating nifti %s with %d slices...", 
			$script_name, $#files+1);
	open(FH, ">" . $script_name) 
		|| die("Error: Cannot create file $script_name");
	print FH "dcm2nii -o $outdir -a y -v n -e n -d n -i n -r n @files";
	print FH " &> $script_name.log\n";
	print FH "NAME=\`ls -alrt $outdir/*.nii.gz | tail -n 1 | ",
				"sed -n \'s/.* \\\(.*\\\)\.nii\.gz/\\1/p\'\`\n";
	print FH "OUT=S_$type\_$echotime.nii.gz\n";
	print FH "sleep 1\n";
	print FH "mv \$NAME.nii.gz $outdir/\$OUT\n";
	print FH "echo $outdir/\$OUT\nexit 0\n";
	close(FH);

	# Execute script
	chmod 0755, $script_name;
	my $out = `$script_name`; chomp($out);
	die("Error: Could not create nifti file $out!") if (!-e $out || -z $out);
	system("rm $script_name $script_name.log");
	print("done.\n");
	return $out;
}


### Main #####################################################################

die("\nUsage: dcm2nii.pl <dcm dir> <dir>\n") if( $#ARGV < 1 );

my $dcmdir = shift;
my $dcmout = shift;
my $ftypes = shift;

# Get echos and locations
my @echotimes;
my @slicelocations;
my $seriesdesc;
my $patientname;
my $cmd = "dcmdump --search EchoTime --search SliceLocation "
		  . "--search SeriesDescription --search PatientName "
		  . "--scan-directories $dcmdir |";
open(ES, $cmd) || die("Error: Cannot open $dcmdir!");
while(my $line = <ES>) {
	if( $line =~ /.*\[(.*)\].*EchoTime/ ) {
		push @echotimes, $1;
	}
	if( $line =~ /.*\[(.*)\].*SliceLocation/ ) {
		push @slicelocations, $1;
	}
	if( $line =~ /.*\[(.*)\].*SeriesDescription/ ) {
		$seriesdesc = $1;
	}
	if( $line =~ /.*\[(.*)\].*PatientName/ ) {
		$patientname = $1;
	}
}
close(ES);

print "Info: $seriesdesc - $patientname\n\n";

@echotimes = sort {$a <=> $b} uniq @echotimes;
printf "Info: Found %d echotimes: @echotimes\n\n", ($#echotimes+1);

@slicelocations = sort {$a <=> $b} uniq @slicelocations; 
printf "Info: Found %d slicelocations: @slicelocations\n\n",
		($#slicelocations+1);

# Get filenames and figure out which one is a dicom file
my $DH;
opendir($DH, $dcmdir) || die("Error: Cannot open $dcmdir!");
my @files = readdir($DH);
closedir($DH);

my @dcmfiles;
my @dcmidx;
my $file;
for $file (@files) {
	my $cmd = "dcmdump --search InstanceNumber $dcmdir/$file |";
	unless( open(ES, $cmd) ) {
		printf("Warning: Cannot open $file. Skipping.\n");
		next;
	}
	my $line = <ES>;
	if( defined($line) && $line =~ /.*\[(.*)\].*InstanceNumber/ ) {
		push @dcmfiles, $file;
		push @dcmidx, $1;
	} else {
		printf("Warning: InstanceNumber tag missing in $file. Skipping.\n");
	}
	close(ES);
}

#my @hu=uniq @dcmidx;
#print "$#dcmidx $#hu\n";

my @sortidx = sort {
	$dcmidx[$a] <=> $dcmidx[$b]
} 0..$#dcmfiles;
#print "@dcmidx[@sortidx]\n";

@dcmfiles = @dcmfiles[@sortidx]; # Sorted list of dicom images
#print "@dcmfiles[0..32]\n";
printf("Info: Found %d images\n\n", ($#dcmfiles+1));

$cmd = "dcmdump $dcmdir/$dcmfiles[0] > $dcmout/$dcmfiles[0].txt";
system($cmd);

# Calculate number of filetypes
my $types;
my $types_float;
if ( defined($ftypes) ) {
	$types = $ftypes;
} else {
	$types_float = ($#dcmfiles+1) / (($#echotimes+1) * ($#slicelocations+1));
	$types = sprintf("%0.0f", $types_float);
	if ( abs($types_float - $types) > 0.1 ) {
		print "$#echotimes $#slicelocations $types_float $types\n";
		die("Error: Found an inconsitent number of files!");
	}
}
printf "Info: Found %0.1f image types\n\n", $types_float;

# Generate Nifti files	
for(my $echo_idx = 0; $echo_idx <= $#echotimes; $echo_idx = $echo_idx+1) {
	my @flist;
	my $echotime = $echotimes[$echo_idx]*1000;

	for (my $type_idx = 0; $type_idx < $types; $type_idx = $type_idx+1) {
		my @selidx;
		my $off = $type_idx + $echo_idx*$types;
		my $inc = $types * ($#echotimes+1);
		for(my $idx = $off; $idx <= $#dcmfiles; $idx=$idx + $inc) {
			push @selidx, $idx;
		}
		#print "@selidx\n";
		push @flist, generate_nifti($type_idx, $echotime,
									$dcmdir, $dcmout, @dcmfiles[@selidx]);
	}
	my $cmd = "FILE=\"";
	for $file (@flist) {
		$cmd = $cmd . "$file ";
	}
	$cmd = $cmd . "\"; fslmerge -t $dcmout/S_$echotime \$FILE; rm \$FILE;";
	system($cmd);
}

exit 0;



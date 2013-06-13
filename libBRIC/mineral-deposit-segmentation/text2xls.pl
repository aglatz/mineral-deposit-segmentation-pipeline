#!/usr/bin/perl -w

# text2xls.pl
# pod at tail

use strict;
use Spreadsheet::WriteExcel;
use Getopt::Long;

# Crank it up
print("\nStarting $0\n");

# Get command line switches and arguments
my %option = ('defaultdelim' => ',');
GetOptions(
    'infile=s'  => \$option{csvin},
    'outfile=s' => \$option{xlsout},
    'delim=s'   => \$option{delim},
    'header!'   => \$option{header},
    );
unless (defined ($option{csvin} && $option{xlsout})) {
    &USAGE();
    exit;
    }
unless (defined $option{delim}) {
    $option{delim} = $option{defaultdelim};
    }
unless (defined $option{header}) {
    $option{header} = 0;
    }

# Do the dirty work
open (CSVFILE, $option{csvin})    or die "Error opening $option{csvin}
+: $!";
my $workbook  = Spreadsheet::WriteExcel -> new($option{xlsout});
my $worksheet = $workbook -> addworksheet();
my $row       = 0;
while (<CSVFILE>) {
    chomp;
    my @field = split("$option{delim}", $_);
    my $column = 0;
    foreach my $token (@field) {
        $worksheet -> write($row, $column, $token);
        $column++;
        }
    $row++;
    }

# Prettify row as header
if($option{header} == 1) {
    my $header = $workbook -> addformat();
        $header -> set_bold();
        $header -> set_align('center');
        $header -> set_bg_color('tan');
        $header -> set_border();
    $worksheet -> set_row(0, undef, $header);
    }

# Optional, unless doing some external action against xls file
$workbook  -> close()         or die "Error closing $workbook: $!";

# Go tell it on the mountain
print("  infile  = $option{csvin}\n");
print("  outfile = $option{xlsout}\n");
print("Finished $0\n\n");

########################################################
sub USAGE {
print <<EOF

Required switches:
  -i or --infile csv_infilename
  -o or --output xls_outfilename

Optional switches:
-d or --delim alt_delimiter
  default delimiter is ','
-h or --header
  Accepts no arguments.
  Bolds and background colors contents of topmost row.

Example:
  csv2xls.pl -i file1.csv -o file2.xls -d * -h

EOF
;
}
########################################################

=head1 Name

 text2xls.pl

=head1 Description

 Convert text file in csv format to Excel binary format

=head1 Usage

 Required switches:
   -i or --infile csv_infilename
   -o or --output xls_outfilename

 Optional switches:
 -d or --delim alt_delimiter
   default delimiter is ','
 -h or --header
   Accepts no arguments.
   Bolds and background colors contents of topmost row.

 Example:
   text2xls.pl -i file1.csv -o file2.xls -d * -h

=head1 Tested

 Spreadsheet::WriteExcel 0.31
 Getopt::Long            2.19
 Perl                    5.00503
 Debian                  2.2r2 "Espy"
 Excel                   9.0.3821 SR-1
                         97 SR-2
 Win32                   5.0.2195 SP-1
                         4.0 SP-6a

=head1 Updated

 2001-04-17   10:45  Renamed as text2xls - avoid conflict with existin
+g csv2xls.
                     Tweaked todos re feedback from jmcnamara (S::WE a
+uthor).
 2001-04-16          Initial working code & posted to PerlMonks

=head1 Todos

 Use Text::xSV by tilly or Text::CSV_XS for parsing text infile
 Review and fix "prettify row header" which somewhat works, but not qu
+ite right.
   Support for "freeze pane" anticipated in future S::WE
 Test with different delimiters.
 Use Pod::Usage to instead of &USAGE(). But appears to only be in Perl
+ 5.6+ 8^(
 Use @headerparms and foreach loop in headers section to reduce duplic
+ation.
 Accept multiple infiles.
 Save as "file.xls" (provided infile is "file") with -r switch.
 Add support for "page delimiter" character for multiple worksheets.
 Add support for $VERSION.

=head1 Author

 ybiC

=head1 Credits

 Thanks to jmcnamara, Petruchio and zdog for suggestions,
 and to jmcnamara for excellent pod in Spreadsheet::WriteExcel,
 and vroom, of course for PM.

=cut

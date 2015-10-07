#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

#"Generates a key file to lookup original protein names in the output file of NetMHC 3.4 from the original 21-mer FASTA file for wildtype(WT) and mutant(MT) proteins"

use vars qw($opt_i $opt_o);
getopts('i:o:');

if(! $opt_i || ! $opt_o){

   print "Usage: $0\n";
   print "-i < input FASTA file with variant sequences for wildtype(WT) and mutant(MT) proteins generated using 'GenerateVariantSequences.pl' and filtered using 'FilterSeq.pl'>\n";
   print "-o < output Key file for lookup>\n";
	die "ERROR : Please make sure all parameters are defined\n";
}

my $input_file = $opt_i;
my $output_file = $opt_o;


if(! -e $input_file){
   die "File $input_file doesn't exist -- fatal.\n";
}

open(IFH,$input_file) || die "Can't open $input_file.\n";
open(OUT,">$output_file") || die "Can't open $output_file for writing.\n";


    my $i = 1;
    while (my $line = <IFH>) {
        chomp $line;
        if ( $line =~ /^>/ ) {
            my $original_name = $line;
            my $new_name      = "Entry_" . $i;
            print OUT join( "\t", $new_name, $original_name ) . "\n";
            $i++;
        }
    }
close(IFH);
close (OUT);

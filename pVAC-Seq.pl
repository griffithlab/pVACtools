#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

sub usage {
    "Usage: $0
\t --input-file FILENAME ; Input TSV File with variants (please provide complete path)
\t --sample-name Name of Sample ; will be used as prefix for output files
\t --netmhc-path Path to local NetMHC3.4 installation
\t --allele Allele name to predict epitope prediction. Multiple alleles can be specified using a comma-separated list. For a list of available alleles, use: netMHC -A
\t --varpeptide-length length of the peptide sequences in the input FASTA file ; default 21
\t --epitope-length length of subpeptides(epitopes) to predict ; Multiple lengths can be specified using a comma-separated list. Typical epitope lengths vary between 8-11.
\t --binding-cutoff ; report only epitopes where the mutant allele has ic50 binding scores below this value ; default 500
\t --min-fc ; Minimum fold change between mutant binding score and wild-type score. The default is 0, which filters no results, but 1 is often a sensible default (requiring that binding is better to the MT than WT)
\t --output-dir DIRNAME ; Output directory for writing all result files
";
}

my $input_file;
my $output_dir;
my $netmhc_path;
my $binding_threshold   = 500;
my $minimum_fold_change = 0;
my $sample_name;
my @allele;
my @epitope_len;
my $peptide_sequence_length = 21;

GetOptions(
    "input-file|i=s"        => \$input_file,
    "sample-name|s=s"       => \$sample_name,
    "netmhc-path|n=s"       => \$netmhc_path,
    "allele|a=s"            => \@allele,
    "varpeptide-length|l:f" => \$peptide_sequence_length,
    "epitope-length|e=s"    => \@epitope_len,
    "binding-cutoff|c:f"    => \$binding_threshold,
    "min-fc|c:f"            => \$minimum_fold_change,
    "output-dir|o=s"        => \$output_dir
) or die usage;

if ( !-e $input_file ) {
    die "ERROR : Input File $input_file doesn't exist.\n" . usage;
}

if ( !$sample_name ) {
    die "ERROR : Sample name not provided.\n" . usage;
}

if ( !-e $netmhc_path ) {
    die "ERROR : NetMHC3.4 installation path not provided or incorrect.\n"
      . usage;
}

if ( !@allele ) {
    die "ERROR : Please specify allele for NetMHC epitope predictions.\n" . usage;
}

if ( !@epitope_len ) {
    die "ERROR : Please specify length for predicting NetMHC epitopes.\n" . usage;
}

if ( !-d $output_dir || !-w $output_dir || !-e $output_dir ) {
    die "ERROR : Please specify output dir (writeable) for results.\n" . usage;
}

my $fasta_file     = $sample_name . '_' . $peptide_sequence_length . '.fa';
my $fasta_key_file = $sample_name . '_' . $peptide_sequence_length . '.key';

### Generate Variant FASTA Seq ###
my $fa_cmd =
"perl bin/GenerateVariantSequences.pl -i $input_file  -o  $output_dir/$fasta_file -l $peptide_sequence_length";
print "\n#GENERATING VARIANT PEPTIDE FASTA FILE: ";
my $stderr_a = `$fa_cmd`;
if ( $stderr_a =~ /^Usage/ ) {
    die "ERROR: Generating variant peptide FASTA file";
}
else { print "COMPLETED"; }

### Generate FASTA Key File ###
my $fa_key_cmd =
"perl bin/GenerateFastaKey.pl -i $output_dir/$fasta_file -o $output_dir/$fasta_key_file";
print "\n#GENERATING FASTA KEY FILE: ";
my $stderr_b = `$fa_key_cmd`;
if   ( $stderr_b =~ /^Usage/ ) { die "ERROR: Generating FASTA key file"; }
else                           { print "COMPLETED"; }

### Valid NETMHC Allele list ###

my $netmhc_a = `$netmhc_path -A`;
my @allele_arr = split( "\n", $netmhc_a );

@allele      = split( /,/, join( ',', @allele ) );
@epitope_len = split( /,/, join( ',', @epitope_len ) );

### Run NETMHC allele-length combinations ###
foreach my $epl (@epitope_len) {
    foreach my $a (@allele) {
        if ( grep { $_ eq $a } @allele_arr ) {
            my $net_out = $sample_name . '.' . $a . '.' . $epl . '.netmhc.xls';
            print "\n#RUNNING NETMHC ON ALLELE $a AND EPITOPE LENGTH $epl\n";
            my $netmhc_cmd =
'$netmhc_path -a $a -l $epl $output_dir/$fasta_file -x $output_dir/$net_out';
            my $tmp =
`$netmhc_path -a $a -l $epl $output_dir/$fasta_file -x $output_dir/$net_out`;

        }
        else {
            print "NetMHC allele not valid. Please check using"
              . $netmhc_path . '-A';
        }
        print "COMPLETED";
    }

}

### Parse NETMHC results and create FOF with file names ##
my $fof = $output_dir . '/' . $sample_name . '.fof';
open( my $fh, '>', $fof ) or die "Could not open file '$fof'";

foreach my $epl (@epitope_len) {
    foreach my $a (@allele) {
        my $net_out =
            $output_dir . '/'
          . $sample_name . '.'
          . $a . '.'
          . $epl
          . '.netmhc.xls';
        my $net_parsed =
            $output_dir . '/'
          . $sample_name . '.'
          . $a . '.'
          . $epl
          . '.netmhc.parsed';
        my $parse_cmd =
"perl bin/ParseOutputNetmhc.pl -i $net_out  -o $net_parsed  -k $output_dir/$fasta_key_file";
        print "\n#PARSING NETMHC OUTPUT:";
        my $stderr_c = `$parse_cmd`;
        if   ( $stderr_c =~ /^Usage/ ) { die "ERROR: Parsing NETMHC output"; }
        else                           { print "COMPLETED"; }
        print $fh $net_parsed . "\n";

    }

}
close $fh;

## Binding Filters ##
my $filt_out = $sample_name . "_filtered.xls";
my $b_cmd =
"perl bin/BindingFilter.pl -i $input_file -f $fof -o $output_dir/$filt_out -c $minimum_fold_change -b $binding_threshold";
print "\n#RUNNING BINDING FILTERS: ";
system($b_cmd);
my $stderr_d = `$b_cmd`;
if   ( $stderr_d =~ /^Usage/ ) { die "ERROR: Running Binding Filters"; }
else                           { print "COMPLETED"; }

print "\n\nDONE: pVac-Seq.pl has completed. File $filt_out contains list of binding-filtered putative neoantigens.\n"
. "We recommend appending coverage information and running bin/CoverageFilers.pl to filter based on sequencing coverage information as well.\n";

__END__

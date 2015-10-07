#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper ;
use Getopt::Long;



sub usage { "Usage: $0
\t --input-file FILENAME ; Input TSV File with variants
\t --sample-name Name of Sample ; will be used as prefix for output files
\t --netmhc-path Path to local NetMHC3.4 installation
\t --allele Allele name to predict epitope prediction. Multiple alleles can be specified using a comma-separated list. For a list of available alleles, use: netMHC -A 
\t --fasta-length length of the peptide sequences in FASTA file ; default 21
\t --epitope-length length of subpeptides(epitopes) to predict ; Multiple lengths can be specified using a comma-separated list
\t --binding-cutoff ; report only epitopes where the mutant allele has ic50 binding scores below this value ; default 500
\t --min-fc ; Minimum fold change between mutant binding score and wild-type score. The default is 0, which filters no results, but 1 is often a sensible default 
   (requiring that binding is better to the MT than WT)
\t --output-file FILENAME ; Filtered File
"}

my $input_file='';
my $output_file;
my $netmhc_path='';
my $binding_threshold = 500;
my $minimum_fold_change = 0;
my $sample_name ='';
my @allele;
my @epitope_len;
my $peptide_sequence_length = 21;



GetOptions ("input-file|i=s" => \$input_file,
			"sample-name|s=s" => \$sample_name,
			"netmhc-path|n=s"=> \$netmhc_path,
			"allele|a=s"	=> \@allele,	
			"fasta-length|l:f" => \$peptide_sequence_length,
			"epitope-length|e=s" => \@epitope_len,
			"binding-cutoff|c:f" => \$binding_threshold,
			"min-fc|c:f" => \$minimum_fold_change,			
            "output-file|o=s" => \$output_file
            )
or die usage;

if(! -e $input_file){
   die  "ERROR : Input File $input_file doesn't exist.\n".
   usage ;
}

if(!$sample_name){
   die  "ERROR : Sample name not provided.\n".
   usage ;
}


if(! -e $netmhc_path){
   die  "ERROR : NetMHC3.4 installation path not provided or incorrect.\n".
   usage ;
}


if(!@allele){
   die  "ERROR : Please specify allele for NetMHC epitope predictions.\n".
   usage ;
}

if(!@epitope_len){
   die  "ERROR : Please specify length for predicting NetMHC epitopes.\n".
   usage ;
}

if(!$output_file){
   die  "ERROR : Please specify output file for filtered results.\n".
   usage ;
}


my $fasta_file = $sample_name.'_'.$peptide_sequence_length.'.fa';
my $fasta_key_file = $sample_name.'_'.$peptide_sequence_length.'.key';



### Generate Variant FASTA Seq ###
system ("perl GenerateVariantSequences.pl -i ".$input_file." -o ".$fasta_file." -l ". $peptide_sequence_length);


### Generate FASTA Key File ###
system ("perl GenerateFastaKey.pl -i ".$fasta_file." -o ".$fasta_key_file);



### Valid NETMHC Allele list ###

my $netmhc_a = `$netmhc_path -A`;
my @allele_arr = split("\n", $netmhc_a);

@allele = split(/,/,join(',',@allele));
@epitope_len = split(/,/,join(',',@epitope_len));



### Run NETMHC allele-length combinations ###
foreach my $epl (@epitope_len)
{
	foreach my $a (@allele)
	{
		if (grep {$_ eq $a} @allele_arr ) {
		my $net_out = $sample_name.'.'.$a.'.'.$epl.'.netmhc.xls';
	 	my $tmp = `$netmhc_path -a $a -l $epl $fasta_file -x $net_out`;
		}	
		else
		{
		print "NetMHC allele not valid. Please check using".$netmhc_path.'-A';
		}
	}	
	
}


### Parse NETMHC results and create FOF with file names ## 
my $fof = $sample_name.'.fof';
open(my $fh, '>', $fof) or die "Could not open file '$fof'";

foreach my $epl (@epitope_len)
{
	foreach my $a (@allele)
	{
		my $net_out = $sample_name.'.'.$a.'.'.$epl.'.netmhc.xls';
		my $net_parsed = $sample_name.'.'.$a.'.'.$epl.'.netmhc.parsed';
		system ("perl ParseOutputNetmhc.pl -i ".$net_out." -o ".$net_parsed." -k ".$fasta_key_file);
		print $fh $net_parsed ."\n";
	}

}	
close $fh;	

## Binding Filters ##
my $b_cmd = "perl BindingFilter.pl -i $input_file -f $fof -o $output_file -c $minimum_fold_change -b $binding_threshold";
system ($b_cmd);







__END__

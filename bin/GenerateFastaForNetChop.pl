#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper ;
use Getopt::Long;


sub usage { "Usage: $0
\t --input-file FILENAME ; Input Filtered File with predicted epitopes (please provide complete path)
\t --output-file DIRNAME ; Output FASTA filename for putative neoepitopes
"}

my $input_file='';
my $output_file;



GetOptions ("input-file|i=s" => \$input_file,
            "output-file|o=s" => \$output_file
            )
or die usage;


if(! -e $input_file){
   die  "ERROR : Input File $input_file doesn't exist.\n".
   usage ;
}

open(IFH,$input_file) || die "Can't open $input_file.\n";

#1 chromosome_name
#2 start	
#3 stop	
#4 reference	
#5 variant	
#6 gene_name	
#7 transcript_name	
#8 amino_acid_change	
#9 ensembl_gene_id	
#10 wildtype_amino_acid_sequence	
#11 GeneName	
#12 HLAallele	
#13 PeptideLength	
#14 SubPeptidePosition	
#15 MTScore	
#16 WTScore	
#17 MTEpitopeSeq	
#18 WTEpitopeSeq	
#19 FoldChange

open(OUT,">$output_file") || die "Can't open $output_file for writing.\n";



while (my $line = <IFH>) {
        chomp $line;
        $line =~ s/[*]$//g;
      #  print $line."\n";
       if($line=~/^chromosome_name/){
       	next;}
        my @arr =  split(/\t/, $line);
        my $gene_name = $arr[10];
        my $hla = $arr[11];
        my $aachange = $arr[7];
        
        my $header = $gene_name.'.'.$aachange;
        print OUT ">".$header."\n";
        print OUT $arr[16]."\n";
        
            }

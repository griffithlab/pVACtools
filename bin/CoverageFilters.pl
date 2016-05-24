#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

sub usage {
    "Usage: $0
\t --input-file FILENAME ; Read count info appended to Binding Filtered File
\t --output-file FILENAME ; Filtered File
\t --norm Set flag if Normal is provided ; --nonorm if no normal sample [DEFAULT]
\t --tdna Set flag if Tumor DNA is provided ; --notdna if no tumor DNA [DEFAULT]
\t --trna Set flag if Tumor RNA is provided ; --notrna if no tumor RNA [DEFAULT]
\t --expn Set flag if Gene Expression is provided ; --noexpn if no expression data [DEFAULT]
\t --normal-cov <VALUE> Normal Coverage Cutoff; sites above this cutoff will be considered ; default 5
\t --tdna-cov <VALUE> Tumor DNA Coverage Cutoff; sites above this cutoff will be considered; default 10 
\t --trna-cov <VALUE> Tumor RNA Coverage Cutoff; sites above this cutoff will be considered; default 10 
\t --normal-vaf <VALUE> Normal VAF Cutoff; sites BELOW this cutoff in normal will be considered ; default 2 
\t --tdna-vaf <VALUE> Tumor DNA VAF Cutoff; sites above this cutoff will be considered ; default 40 
\t --trna-vaf <VALUE> Tumor RNA VAF Cutoff; sites above this cutoff will be considered ; default 40
\t --expn-val <VALUE> Gene Expression (FPKM) Cutoff ; default 1";
}

my $input_file  = '';
my $output_file = '';
my $norm_cov    = 5;
my $norm_vaf    = 2;
my $dna_tum_cov = 10;
my $dna_tum_vaf = 40;
my $rna_tum_cov = 10;
my $rna_tum_vaf = 40;
my $expn_val    = 1;

my $norm = 0;    # option variable with default value (true)
my $tdna = 0;
my $trna = 0;
my $expn = 0;

GetOptions(
    "input-file|i=s"  => \$input_file,
    "normal-cov:f"    => \$norm_cov,
    "normal-vaf:f"    => \$norm_vaf,
    "tdna-cov:f"      => \$dna_tum_cov,
    "tdna-vaf:f"      => \$dna_tum_vaf,
    "trna-cov:f"      => \$rna_tum_cov,
    "trna-vaf:f"      => \$rna_tum_vaf,
    "expn-val:f"      => \$expn_val,
    "norm!"           => \$norm,
    "tdna!"           => \$tdna,
    "trna!"           => \$trna,
    "expn!"           => \$expn,
    "output-file|o=s" => \$output_file
) or die usage;

if ( !-e $input_file ) {
    die "Input File $input_file doesn't exist.\n";
}

if ( !$output_file ) {
    die "Output File not specified.\n";
}

if ( ( defined $norm ) && ( $norm eq 1 ) ) {
    print "Normal DNA is provided. Normal Coverage Cutoff is "
      . $norm_cov
      . " and VAF cutoff is "
      . $norm_vaf . "\n";
}
else {
    print "No Normal DNA will be used for filtering\n";

}

if ( ( defined $tdna ) && ( $tdna eq 1 ) ) {
    print "Tumor DNA is provided. Tumor DNA Coverage Cutoff is "
      . $dna_tum_cov
      . " and VAF cutoff is "
      . $dna_tum_vaf . "\n";
}
else {
    print "No Tumor DNA will be used for filtering.\n";

}

if ( ( defined $trna ) && ( $trna eq 1 ) ) {
    print "Tumor RNA is provided. Tumor RNA Coverage Cutoff is "
      . $rna_tum_cov
      . " and VAF cutoff is "
      . $rna_tum_vaf . "\n";
}
else {
    print "No Tumor RNA will be used for filtering\n";

}

if ( ( defined $expn_val ) && ( $expn eq 1 ) ) {
    print "Gene expression values will be filtered using cutoff of "
      . $expn_val . "\n";
}
else {
    print "No Gene expression values provided.\n";

}

#### INPUT AND OUTPUT FILE FORMAT ##
#1	chromosome_name
#2	start
#3	stop
#4	reference
#5	variant
#6	gene_name
#7	transcript_name
#8	amino_acid_change
#9	ensembl_gene_id
#10	wildtype_amino_acid_sequence
#11	GeneName
#12	HLAallele
#13	PeptideLength
#14	SubPeptidePosition
#15	MTScore
#16	WTScore
#17	MTEpitopeSeq
#18	WTEpitopeSeq
#19	FoldChange
#20	NormalRefCount
#21 NormalVarCount
#22 TumorDNARefCount
#23 TumorDNAVarCount
#24 TumorRNARefCount
#25 TumorRNAVarCount
#26 GeneExpFPKM

## TODO : add functionality for gene tracking file#

open( IFH, $input_file )     || die "Can't open $input_file.\n";
open( OUT, ">$output_file" ) || die "Can't open $output_file for writing.\n";

while ( my $line = <IFH> ) {
    chomp $line;
    my @arr = split( "\t", $line );
    if ( $line =~ /^chr/ ) {
        print OUT $line . "\n";
        next;
    }
    else {

        my $nref   = $arr[19];
        my $nvar   = $arr[20];
        my $nc     = $nref + $nvar;
        my $nvaf   = ( $nvar / ( $nc + 0.00001 ) ) * 100;
        my $nvaf_r = sprintf( "%.3f", $nvaf );

        my $tdref = $arr[21];
        my $tdvar = $arr[22];
        my $tdc   = $tdref + $tdvar;
        my $tdvaf = ( $tdvar / ( $tdc + 0.00001 ) ) * 100;

        my $tdvaf_r = sprintf( "%.3f", $tdvaf );

        my $trref   = $arr[23];
        my $trvar   = $arr[24];
        my $trc     = $trref + $trvar;
        my $trvaf   = ( $trvar / ( $trc + 0.00001 ) ) * 100;
        my $trvaf_r = sprintf( "%.3f", $trvaf );

        my $fpkmval = $arr[25];

        my ( $flag_n, $flag_dna, $flag_rna, $flag_expn );
        my @flag_arr;

        #################### NORMAL FILTER ####################

        if ( $norm eq 1 ) {

            if ( $nc >= $norm_cov && $nvaf_r <= $norm_vaf ) {
                $flag_n = 1;
            }
            else {
                $flag_n = 0;

      #				print $line."\t".$nc."<". $norm_cov."\t".$nvaf_r.">=".$norm_vaf."\n";
            }

            push( @flag_arr, $flag_n );
        }
        #################### DNA FILTER ####################

        if ( $tdna eq 1 ) {
            if ( $tdc >= $dna_tum_cov && $tdvaf_r >= $dna_tum_vaf ) {

                $flag_dna = 1;
            }

            else

            {
                $flag_dna = 0;

            }

            push( @flag_arr, $flag_dna );
        }
        #################### RNA FILTER ####################

        if ( $trna eq 1 ) {
            if ( $trc >= $rna_tum_cov && $trvaf_r >= $rna_tum_vaf )

            {

                $flag_rna = 1;
            }

            else {

                $flag_rna = 0;
            }

            push( @flag_arr, $flag_rna );
        }

        #################### GENE EXPRESSION FILTER ####################

        if ( $expn eq 1 )

        {
            if ( $fpkmval >= $expn_val ) {
                $flag_expn = 1;

            }
            else {
                $flag_expn = 0;

            }
            push( @flag_arr, $flag_expn );
        }

        ####################

        if ( grep { $_ eq '0' } @flag_arr ) {

            #		print Dumper \@flag_arr;
            next;
        }
        else {
            print OUT $line . "\n";

        }

    }
}

__END__
				
	

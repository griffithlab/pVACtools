# Based on myParseReportsTopGenesTopAlleleStringent_4.pl by TWYLIE (Feb 2012)
#Takes in a FOF with parsed NetMHC files for different allele-length combinations and outputs best candidates per gene based on binding affinity.
use strict;
use warnings;
use File::Basename qw(basename);
use Getopt::Std;


use vars qw($opt_i $opt_f $opt_o $opt_c $opt_b);
getopts('i:f:o:c:b:');


if(! $opt_f || ! $opt_i || ! $opt_o || ! $opt_c || ! $opt_b){

   print "Usage: $0\n";
   print "-i <input list of variants>\n";
   print "-f <FOF containing list of parsed epitope files for different allele-length combinations (same sample)>\n";
   print "-o <Output .xls file containing list of filtered epitopes based on binding affinity for each allele-length combination per gene>\n";
   print "-c <Minimum fold change between mutant binding score and wild-type score. The default is 0, which filters no results, but 1 is often a sensible default 
   (requiring that binding is better to the MT than WT)>\n";
   print "-b <report only epitopes where the mutant allele has ic50 binding scores below this value ; default 500>\n";

	die "ERROR : Please make sure all parameters are defined\n";
}

my $variant_file = $opt_i;
my $fof_file = $opt_f;
my $output_file = $opt_o;
my $binding_threshold = $opt_b || 500;
my $minimum_fold_change = $opt_c || 0;


#$binding_threshold = 500;
#$minimum_fold_change = 0;

if(! -e $variant_file){
   die "Variant File $variant_file doesn't exist -- fatal.\n";
}


    my %variants;
    open(IFH,$variant_file) || die "Can't open $variant_file.\n";
    while (my $line = <IFH>) {
            chomp $line;
            #header?
            if($line=~/^chromosome_name/){
                $variants{"header"} = $line;
                next;
            }
            my @F = split("\t",$line); 
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
                       
            my $gene = $F[5];
            my $amino_acid_change = $F[7];
            $variants{join("\t",($gene,$amino_acid_change))} = $line;
        }
       if(defined($variants{"header"})){
		open(OUT,">$output_file") || die "Can't open $output_file for writing.\n";	
            print OUT join("\t",($variants{"header"},
                                     "GeneName","MutationEffect","HLAallele","PeptideLength",
                                     'SubPeptidePosition','MT Score',
                                     'WT Score','MT Epitope Seq',
                                     'WT Epitope Seq','Fold Change')) . "\n";
        }

    # read in the netMHC predictions and filter
        open(FOF,$fof_file) || die "Can't open $fof_file.\n";


    my %prediction;
    my $threshold = $binding_threshold;
    while (my $file = <FOF>) {
        chomp $file;
        my $basename = basename( $file);
        # TEST.HLA-A01:01.9.netmhcM.parsed
        my @f      = split( /\./, $basename );
        my $sample = $f[0];
       $sample =~ s/_netmhc//g;
        my $allele = $f[1];
        my $length = $f[2];

        my $mode = 'filtered';
        my $i = 0;
        
        open(FILE,$file)
  		or die "Could not open file '$file' $!";
 		FILELINE:
   		while (<FILE>) {
        chomp;
        my $line = $_;
        $i++;
        next if ($i == 1);  # skips headers....
        my (
            $gene_name,
            $point_mutation,
            $sub_peptide_position,
            $mt_score,
            $wt_score,
            $mt_epitope_seq,
            $wt_epitope_seq,
            $fold_change,
           ) = split( "\t", $line );


        $gene_name = {
                      gene_name            => $gene_name,
                      allele               => $allele,
                      point_mutation       => $point_mutation,
                      sub_peptide_position => $sub_peptide_position,
                      mt_score             => $mt_score,
                      wt_score             => $wt_score,
                      mt_epitope_seq       => $mt_epitope_seq,
                      wt_epitope_seq       => $wt_epitope_seq,
                      fold_change          => $fold_change,
                     };
        push( @{ $prediction{$mode}->{$sample}->{$length}->{genes} }, \$gene_name );
    }
    

}
close (FOF);


##

my %best;

MODE:
foreach my $mode (sort keys %prediction) {
    SAMPLE:
    foreach my $sample (sort keys %{ $prediction{$mode} }) {
       LENGTH:
        foreach my $length (sort keys %{ $prediction{$mode}->{$sample} }) {
            GENES:
            foreach my $gene (sort @{ $prediction{$mode}->{$sample}->{$length}->{genes} }) {
                # BEST
                unless( !$best{$sample}->{$$gene->{gene_name}}->{SCORE} ) {
                    if ($$gene->{mt_score} < $best{$sample}->{$$gene->{gene_name}}->{SCORE}) {
                        $best{$sample}->{$$gene->{gene_name}}->{SCORE} = $$gene->{mt_score};
                        $best{$sample}->{$$gene->{gene_name}}->{GENES} = [];
                        $$gene->{sample} = $sample;
                        $$gene->{length} = $length;
                        $$gene->{mode}   = $mode;
                        push( @{ $best{$sample}->{$$gene->{gene_name}}->{GENES} }, $gene );
                    }
                    elsif ($$gene->{mt_score} == $best{$sample}->{$$gene->{gene_name}}->{SCORE}) {
                        $best{$sample}->{$$gene->{gene_name}}->{SCORE} = $$gene->{mt_score};
                        $$gene->{sample} = $sample;
                        $$gene->{length} = $length;
                        $$gene->{mode}   = $mode;
                        push( @{ $best{$sample}->{$$gene->{gene_name}}->{GENES} }, $gene );
                    }
                }
                else {
                    $best{$sample}->{$$gene->{gene_name}}->{SCORE} = $$gene->{mt_score};
                    $$gene->{sample} = $sample;
                    $$gene->{length} = $length;
                    $$gene->{mode}   = $mode;
                    push( @{ $best{$sample}->{$$gene->{gene_name}}->{GENES} }, $gene );
                }

            }
        }
    }
}

# REPORTING
    foreach my $sample (sort keys %best) {
        foreach my $gene (sort keys %{ $best{$sample} }) {
            foreach my $entry (@{ $best{$sample}->{$gene}->{GENES} }) {
                if (($$entry->{mt_score} < $threshold) && ($$entry->{fold_change} > $minimum_fold_change)){
                    #write out epiptopes
    #                print OUT $entry;
                    #also write out variants with epitopes appended
                    if(defined($variant_file)){
                        my $key = join("\t",$gene,$$entry->{point_mutation});
                        if(defined($variants{$key})){                            
                            print OUT join("\t", ($variants{$key},
                                                      $gene,$$entry->{point_mutation},$$entry->{allele},
                                                      $$entry->{length},$$entry->{sub_peptide_position},$$entry->{mt_score},
                                                      $$entry->{wt_score},$$entry->{mt_epitope_seq},
                                                      $$entry->{wt_epitope_seq},$$entry->{fold_change})) . "\n";

                        } else {
                            
                            print "Couldn't find variant for " . $gene . " " . $$entry->{point_mutation} . "in variant_file";
                        }
                    }
                }
            }
        }
    }

close (OUT);
close (IFH);


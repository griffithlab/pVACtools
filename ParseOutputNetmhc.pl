# FOR NETMHC3.4 : Parses output from NetMHC3.4 for MHC Class I epitope prediction; Uses a special key file that could be generated using GenerateFastaKey.pl --help.The parsed TSV file contains predictions for only those epitopes that contain the mutant SNV
use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

#"Generates a key file to lookup original protein names in the output file of NetMHC 3.4 from the original 21-mer FASTA file for wildtype(WT) and mutant(MT) proteins"

use vars qw($opt_i $opt_o $opt_k);
getopts('i:o:k:');

if(! $opt_i || ! $opt_o || ! $opt_k) {

   print "Usage: $0\n";
   print "-i < Raw output file from Netmhc>\n";
   print "-o  < Parsed output file >\n";
   print "-k < Key file for lookup of FASTA IDs>\n";
  
	die "ERROR : Please make sure all parameters are defined\n";
}

my $input_file = $opt_i;
my $output_file = $opt_o;
my $key_file = $opt_k ; 

if(! -e $input_file){
   die "File $input_file doesn't exist -- fatal.\n";
}

open(IFH,$input_file) || die "Can't open $input_file.\n";
open(OUT,">$output_file") || die "Can't open $output_file for writing.\n";

   my $key_hash = key_hash();
   my %position_score;
   my ($netmhc_results, $epitope_seq) = make_hashes_from_input($key_hash);

 	print OUT join("\t",
        "Gene Name",
        "Point Mutation",
        "Sub-peptide Position",
        "MT score",
        "WT score",
        "MT epitope seq",
        "WT epitope seq",
        "Fold change"
    ) . "\n";
    
    
    my $protein_type = 'MT';
    for my $protein_name (sort keys %{$netmhc_results->{$protein_type}}) {
        for my $variant_aa (sort keys %{$netmhc_results->{$protein_type}->{$protein_name}}) {
            %position_score = %{$netmhc_results->{$protein_type}{$protein_name}{$variant_aa}};
            my @positions = sort {$position_score{$a} <=> $position_score{$b}} keys %position_score;
            my $total_positions = scalar @positions;

                for (my $i = 0; $i < $total_positions; $i++) {

                    if ($epitope_seq->{'MT'}->{$protein_name}->{$variant_aa}->{$positions[$i]} ne
                        $epitope_seq->{'WT'}->{$protein_name}->{$variant_aa}->{$positions[$i]})
                        # Filtering if mutant amino acid present
                    {
                        my $position = $positions[$i];
                        
                print OUT join("\t", $protein_name, $variant_aa, $position, $position_score{$position}) . "\t";
   			    print OUT $netmhc_results->{'WT'}->{$protein_name}->{$variant_aa}->{$position} . "\t";
			    print OUT $epitope_seq->{'MT'}->{$protein_name}->{$variant_aa}->{$position} . "\t";
   				print OUT $epitope_seq->{'WT'}->{$protein_name}->{$variant_aa}->{$position} . "\t";
   				my $fold_change =
    		    $netmhc_results->{'WT'}->{$protein_name}->{$variant_aa}->{$position} / $position_score{$position};
  				my $rounded_FC = sprintf("%.3f", $fold_change);
 		   		print OUT $rounded_FC . "\n";
                        
                    }
                }
            }

            }


sub key_hash {
	open(KEY,$key_file) || die "Can't open $key_file.\n";
    my %key_hash;
    
    while (my $keyline = <KEY>) {
        chomp $keyline;
        #Entry_1	>WT.GSTP1.R187W
        my ($new_name, $original_name) = split(/\t/, $keyline);
        $original_name =~ s/>//g;
        $key_hash{$new_name} = ();
        $key_hash{$new_name}{'name'} = $original_name;

    }

    close(KEY);

    return \%key_hash;
}

sub make_hashes_from_input {
    my $key_hash = shift;
    my (%netmhc_results, %epitope_seq);
    while (my $line = <IFH>) {
    chomp $line;
    	if ($line =~ /^Entry/) {
        my @result_arr = split(/\t/, $line);
        my $position         = $result_arr[1];
        my $score            = $result_arr[3];
        my $epitope          = $result_arr[2];
        my $protein_new_name = $result_arr[0];
		
        my (@protein_arr);
        if ($line =~ /^Entry/) {
            if (exists($key_hash->{$protein_new_name})) {
                my $protein = $key_hash->{$protein_new_name}{'name'};
                @protein_arr = split(/\./, $protein);
            }
        }
    
        else {
            next;
        }
        my $protein_type = $protein_arr[0];
        my $protein_name = $protein_arr[1];
        my $variant_aa =  $protein_arr[2];

        $netmhc_results{$protein_type}{$protein_name}{$variant_aa}{$position} = $score;
        $epitope_seq{$protein_type}{$protein_name}{$variant_aa}{$position} = $epitope;
    } }

    return \%netmhc_results, \%epitope_seq;
}
close (IFH);
close (OUT);
##### END ####



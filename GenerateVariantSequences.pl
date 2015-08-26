use strict;
use warnings;
use Data::Dumper ;
use Getopt::Std;


use vars qw($opt_i $opt_l $opt_o);
getopts('i:l:o:');
if(! $opt_i || ! $opt_l || ! $opt_o){

   print "Usage: $0\n";
   print "-i <input list of variants>\n";
   print "-l <length of the peptide sequences>\n";
   print "-o <output FASTA file>\n";
	die "ERROR : Please make sure all parameters are defined\n";
}

my $input_file = $opt_i;
my $peptide_sequence_length = $opt_l;
my $output_file = $opt_o;


if(! -e $input_file){
   die "File $input_file doesn't exist -- fatal.\n";
}

my $var_type = 'missense';
open(IFH,$input_file) || die "Can't open $input_file.\n";

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

open(OUT,">$output_file") || die "Can't open $output_file for writing.\n";

while (my $line = <IFH>) {
        chomp $line;
        $line =~ s/[*]$//g;
      #  print $line."\n";
        my @protein_arr =  split(/\t/, $line);
        if ( $protein_arr[7] =~ /^([A-Z])(\d+)([A-Z])/ ) {
            my ($wildtype_amino_acid, $position, $mutant_amino_acid) = ($1, $2 - 1, $3);
       		my $wildtype_sequence = $protein_arr[9];

            my @arr_wildtype_sequence = split('',$wildtype_sequence);
            if ($wildtype_amino_acid ne $arr_wildtype_sequence[$position]) {
                next;
            }
            else {
                my ($mutation_position, @wildtype_arr) = get_wildtype_subsequence_for_printing($position, \@arr_wildtype_sequence, \@protein_arr,$peptide_sequence_length);
                my @mutant_arr = @wildtype_arr;
         
                $mutant_arr[$mutation_position] = $mutant_amino_acid;
                print_to_output(\@wildtype_arr, \@mutant_arr, \@protein_arr, $position);
            }
        }
    }

    close(IFH);
	close (OUT);
  #  return 1;


sub get_wildtype_subsequence_for_printing {
    my ($position, $arr_wildtype_sequence_ref, $protein_arr,$peptide_sequence_length) = @_;
    my @arr_wildtype_sequence = @$arr_wildtype_sequence_ref;
    #If the wildtype sequence is shorter than the desired peptide sequence
    #length we use the wildtype sequence length instead so that the extraction
    #algorithm below works correctly
    if (scalar(@arr_wildtype_sequence) < $peptide_sequence_length) {
        $peptide_sequence_length = scalar(@arr_wildtype_sequence);
    }
    # We want to extract a subset from @arr_wildtype_sequence that is
    # $peptide_sequence_length long so that the $position ends
    # up in the middle of the extracted sequence.
    # If the $position is too far toward the beginning or end of
    # @arr_wildtype_sequence there aren't enough amino acids on one side
    # to achieve this.
    my (@wildtype_arr, $mutation_position);
    my $one_flanking_sequence_length = ($peptide_sequence_length - 1) / 2;
    if (distance_from_start($position, @arr_wildtype_sequence) < $one_flanking_sequence_length) {
        @wildtype_arr = @arr_wildtype_sequence[ 0 ... ($peptide_sequence_length - 1) ];
        $mutation_position = $position;
    }
    elsif (distance_from_end($position, @arr_wildtype_sequence) < $one_flanking_sequence_length) {
        @wildtype_arr = @arr_wildtype_sequence[ ($#arr_wildtype_sequence - $peptide_sequence_length + 1) ... $#arr_wildtype_sequence];
        $mutation_position = $peptide_sequence_length - distance_from_end($position, @arr_wildtype_sequence) - 1;
    }
    elsif (
        (distance_from_start($position, @arr_wildtype_sequence) >= $one_flanking_sequence_length) &&
        (distance_from_end($position, @arr_wildtype_sequence) >= $one_flanking_sequence_length)
    ) {
        @wildtype_arr = @arr_wildtype_sequence[ ($position - $one_flanking_sequence_length) ... ($position + $one_flanking_sequence_length) ];
        $mutation_position = ($peptide_sequence_length - 1) / 2;
    }
    else {
       die: Something went wrong during the retrieval of the wildtype sequence at $position $!;
  #          $protein_arr->[0],
   #         $protein_arr->[1],
    #        $protein_arr->[2],
     #       join('', @arr_wildtype_sequence)
      #  );
    }

    return $mutation_position, @wildtype_arr;
	}

#This subroutine is a bit funky but it was designed that way to mirror
#distance_from_end to increase code readability from the caller's perspective
sub distance_from_start {
    my $position = shift;
    my @array = @_;

    return $position;
}

sub distance_from_end {
    my $position = shift;
    my @array = @_;

    return $#array - $position;
}


sub print_to_output {
    my ($wildtype_arr_ref,$mutant_arr_ref,$protein_arr_ref,$position) = @_;

   if ($wildtype_arr_ref) {
        for my $designation ('WT', 'MT') {
            my $arr = $designation eq 'WT' ? $wildtype_arr_ref : $mutant_arr_ref;   
            # >WT.OR51I1.p.R124C
        print_fasta_entry(
                designation => $designation,
                arr         => $arr,
                protein_arr => $protein_arr_ref,
                );

        }
    }
    else {
        print  "NULL"."\t".$position."\n";
    }
}

sub print_fasta_entry {
    my %params = @_;
    open(OUTAPP,">>$output_file") || die "Can't open $output_file for writing.\n";
    
    my ($fasta_defline, $fasta_sequence);
    my $identifier = $params{'designation'} . "." . $params{protein_arr}->[5] . "." . $params{protein_arr}->[7];
    $fasta_defline = ">$identifier";
    $fasta_sequence = join "", @{$params{arr}};

    print OUTAPP  "$fasta_defline\n";
    print OUTAPP "$fasta_sequence\n";

    close(OUTAPP);
}


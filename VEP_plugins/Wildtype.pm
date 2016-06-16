=head1 NAME

 Wildtype

=head1 SYNOPSIS

 mv Wildtype.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin Wildtype

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 provides the wildtype protein sequence of a transcript.

=cut

package Wildtype;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub version {
    return '1.0';
}

sub feature_types {
    return ['Transcript'];
}

sub get_header_info {
    return {
        WildtypeProtein     => "The normal, non-mutated protein sequence",
    };
}

sub run {
    my ($self, $tva) = @_;

    my $tv = $tva->transcript_variation;
    my $tr = $tv->transcript;
    my $cds_seq = defined($tr->{_variation_effect_feature_cache}) ? $tr->{_variation_effect_feature_cache}->{translateable_seq} : $tr->translateable_seq;
    my $codon_seq = Bio::Seq->new(
        -seq      => $cds_seq,
        -moltype  => 'dna',
        -alphabet => 'dna'
    );

    #get codon table
    my $codon_table;
    if(defined($tr->{_variation_effect_feature_cache})) {
        $codon_table = $tr->{_variation_effect_feature_cache}->{codon_table} || 1;
    }
    else {
        my ($attrib) = @{$tr->slice->get_all_Attributes('codon_table')};
        $codon_table = $attrib ? $attrib->value || 1 : 1;
    }

    # translate
    my $new_pep = $codon_seq->translate(undef, undef, undef, $codon_table)->seq();
    $new_pep =~ s/\*.*//;

    return {
        WildtypeProtein => $new_pep,
    };
}

1;


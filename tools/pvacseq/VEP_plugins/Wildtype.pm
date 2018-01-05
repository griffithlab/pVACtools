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

    my $pep;

    # this handles the case where VEP hasn't precached the translation data
    if(!$tr->{_variation_effect_feature_cache}) {
        
        # initialize so we don't try again on failure
        $tr->{_variation_effect_feature_cache} = {};

        if(my $translation = $tr->translate) {

            # cache on transcript object so we're not refetching every call
            $tr->{_variation_effect_feature_cache}->{peptide} = $translation->seq;
        }
    }

    if(my $pep = $tr->{_variation_effect_feature_cache}->{peptide}) {
        return { WildtypeProtein => $pep };
    }
    # some transcripts won't (or shouldn't) have translations e.g. non-protein coding transcripts
    else {
        return {};
    }
}

1;


import vcfpy
import sys
import os
from pvactools.lib.csq_parser import CsqParser
from Bio.Seq import translate
import pvactools.lib.run_utils

class ProximalVariant:
    #flanking_bases is the number of bases (not amino acids!) to search on each side of a variant position
    def __init__(self, proximal_variants_vcf, pass_only, flanking_bases):
        if not os.path.exists(proximal_variants_vcf + '.tbi'):
            sys.exit('No .tbi file found for proximal variants VCF. Proximal variants VCF needs to be tabix indexed.')

        self.proximal_variants_vcf = vcfpy.Reader.from_path(proximal_variants_vcf)

        if len(self.proximal_variants_vcf.header.samples.names) > 1:
            raise Exception("Proximal variants vcf {} contains more than one sample.".format(proximal_variants_vcf))

        if 'CSQ' in self.proximal_variants_vcf.header.info_ids():
            csq_header = self.proximal_variants_vcf.header.get_info_field_info('CSQ')
        else:
            sys.exit('Input VCF does not contain a CSQ header. Please annotate the VCF with VEP before running it.')

        if csq_header is None:
            sys.exit('Failed to extract format string from info description for tag (CSQ)')
        else:
            self.csq_parser = CsqParser(csq_header.description)

        self.pass_only = pass_only
        self.flanking_bases = flanking_bases

    def extract(self, somatic_variant, alt, transcript):
        (phased_somatic_variant, potential_proximal_variants) = self.find_phased_somatic_variant_and_potential_proximal_variants(somatic_variant, alt, transcript)

        if phased_somatic_variant is None:
            print("Warning: Main somatic variant not found in phased variants file: {}, {}".format(somatic_variant, alt))
            return []

        if len(potential_proximal_variants) == 0:
            return []

        proximal_variants = []
        sample = self.proximal_variants_vcf.header.samples.names[0]
        phased_somatic_variant_genotype = phased_somatic_variant.call_for_sample[sample]
        if 'HP' in phased_somatic_variant.FORMAT:
            somatic_phasing = phased_somatic_variant_genotype.data['HP']
            for (entry, csq_entry, proximal_alt) in potential_proximal_variants:
                proximal_variant_genotype = entry.call_for_sample[sample]
                #identify variants that are in phase with the phased_somatic_variant
                if 'HP' in entry.FORMAT:
                    proximal_variant_phasing = proximal_variant_genotype.data['HP']
                    if proximal_variant_phasing == somatic_phasing:
                        proximal_variants.append([entry, csq_entry, proximal_alt])
                #proximal variant is hom var
                elif proximal_variant_genotype.is_variant and not proximal_variant_genotype.is_het:
                    #main somatic variant is het
                    if phased_somatic_variant_genotype.is_het:
                        proximal_variants.append([entry, csq_entry, proximal_alt])
                    #main somatic variant is hom var
                    if phased_somatic_variant_genotype.is_variant and not phased_somatic_variant_genotype.is_het:
                        proximal_variants.append([entry, csq_entry, proximal_alt])
        else:
            for (entry, csq_entry, proximal_alt) in potential_proximal_variants:
                proximal_variant_genotype = entry.call_for_sample[sample]
                #proximal variant is hom var
                if proximal_variant_genotype.is_variant and not proximal_variant_genotype.is_het:
                    #main somatic variant is het
                    if phased_somatic_variant_genotype.is_het:
                        proximal_variants.append([entry, csq_entry, proximal_alt])
                    #main somatic variant is hom var
                    if phased_somatic_variant_genotype.is_variant and not phased_somatic_variant_genotype.is_het:
                        proximal_variants.append([entry, csq_entry, proximal_alt])

        return proximal_variants

    def find_phased_somatic_variant_and_potential_proximal_variants(self, somatic_variant, alt, transcript):
        potential_proximal_variants = []
        phased_somatic_variant = None
        try:
            entries = self.proximal_variants_vcf.fetch(somatic_variant.CHROM, somatic_variant.begin - self.flanking_bases, somatic_variant.affected_end + self.flanking_bases)
        except ValueError as e:
            return (phased_somatic_variant, potential_proximal_variants)

        for entry in entries:
            if self.pass_only:
                filt = entry.FILTER
                if not (filt is None or len(filt) == 0 or filt == ['PASS']):
                    continue

            for proximal_alt in entry.ALT:
                proximal_alt = proximal_alt.value
                if entry.begin == somatic_variant.begin and entry.end == somatic_variant.end and proximal_alt == alt:
                    phased_somatic_variant = entry
                    continue

                if 'CSQ' not in entry.INFO:
                    print("Warning: Proximal variant is not VEP annotated and will be skipped: {} {}".format(entry.CHROM, entry.POS))
                    continue

                alleles_dict = self.csq_parser.resolve_alleles(entry)
                csq_entries = self.csq_parser.parse_csq_entries_for_allele(entry.INFO['CSQ'], proximal_alt)
                if len(csq_entries) == 0:
                    csq_allele = alleles_dict[str(proximal_alt)]
                    csq_entries = self.csq_parser.parse_csq_entries_for_allele(entry.INFO['CSQ'], csq_allele)

                    if len(csq_entries) == 0:
                        print("Warning: Proximal variant does not contain any VEP annotations for alternate allele and will be skipped: {} {} {}".format(entry.CHROM, entry.POS, proximal_alt))
                        continue

                picked_csq_entry = None
                for csq_entry in csq_entries:
                    if csq_entry['Feature'] == transcript:
                        picked_csq_entry = csq_entry
                if picked_csq_entry is None:
                    print("Warning: Proximal variant has no transcript annotation for somatic variant of interest transcript {} and will be skipped: {}".format(transcript, entry))
                    continue

                consequences = {consequence.lower() for consequence in picked_csq_entry['Consequence'].split('&')}
                if 'missense_variant' not in consequences:
                    print("Warning: Proximal variant is not a missense mutation and will be skipped: {} {}".format(entry.CHROM, entry.POS))
                    continue

                potential_proximal_variants.append([entry, picked_csq_entry, proximal_alt])

        return (phased_somatic_variant, potential_proximal_variants)

    @classmethod
    def combine_conflicting_variants(cls, codon_changes):
        codon = list(codon_changes[0].split('/')[0].lower())
        modified_positions = []
        for codon_change in codon_changes:
            (old_codon, new_codon) = codon_change.split('/')
            change_positions = [i for i in range(len(old_codon)) if old_codon[i] != new_codon[i]]
            for position in change_positions:
                if position in modified_positions:
                    print("Warning: position has already been modified")
                codon[position] = new_codon[position].lower()
                modified_positions.append(position)
        return translate("".join(codon))


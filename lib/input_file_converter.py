import vcf
import csv
import sys
import os
from abc import ABCMeta
from collections import OrderedDict
from lib.csq_parser import CsqParser
import lib.utils
from lib.proximal_variant import ProximalVariant
import lib.utils
import binascii
import re

class InputFileConverter(metaclass=ABCMeta):
    def __init__(self, **kwargs):
        self.input_file  = kwargs['input_file']
        self.output_file = kwargs['output_file']

    def output_headers(self):
        return[
            'chromosome_name',
            'start',
            'stop',
            'reference',
            'variant',
            'gene_name',
            'transcript_name',
            'transcript_support_level',
            'amino_acid_change',
            'codon_change',
            'ensembl_gene_id',
            'hgvsc',
            'hgvsp',
            'wildtype_amino_acid_sequence',
            'downstream_amino_acid_sequence',
            'fusion_amino_acid_sequence',
            'variant_type',
            'protein_position',
            'transcript_expression',
            'gene_expression',
            'normal_depth',
            'normal_vaf',
            'tdna_depth',
            'tdna_vaf',
            'trna_depth',
            'trna_vaf',
            'index',
            'protein_length_change',
        ]

class VcfConverter(InputFileConverter):
    def __init__(self, **kwargs):
        InputFileConverter.__init__(self, **kwargs)
        self.gene_expn_file              = kwargs.pop('gene_expn_file', None)
        self.transcript_expn_file        = kwargs.pop('transcript_expn_file', None)
        self.normal_snvs_coverage_file   = kwargs.pop('normal_snvs_coverage_file', None)
        self.normal_indels_coverage_file = kwargs.pop('normal_indels_coverage_file', None)
        self.tdna_snvs_coverage_file     = kwargs.pop('tdna_snvs_coverage_file', None)
        self.tdna_indels_coverage_file   = kwargs.pop('tdna_indels_coverage_file', None)
        self.trna_snvs_coverage_file     = kwargs.pop('trna_snvs_coverage_file', None)
        self.trna_indels_coverage_file   = kwargs.pop('trna_indels_coverage_file', None)
        self.pass_only                   = kwargs.pop('pass_only', False)
        self.sample_name        = kwargs.pop('sample_name', None)
        self.normal_sample_name = kwargs.pop('normal_sample_name', None)
        self.proximal_variants_vcf = kwargs.pop('proximal_variants_vcf', None)
        self.proximal_variants_tsv = kwargs.pop('proximal_variants_tsv', None)
        self.peptide_length = kwargs.pop('peptide_length', None)
        if self.proximal_variants_vcf and not (self.proximal_variants_tsv and self.peptide_length):
            sys.exit("A proximal variants TSV output path and peptide length need to be specified if a proximal variants input VCF is provided")
        if self.proximal_variants_vcf and not os.path.exists(self.proximal_variants_vcf + '.tbi'):
            sys.exit('No .tbi file found for proximal variants VCF {}. Proximal variants VCF needs to be tabix indexed'.format(self.proximal_variants_vcf))
        if self.proximal_variants_vcf and not os.path.exists(self.input_file + '.tbi'):
            sys.exit('No .tbi file found for input VCF {}. Input VCF needs to be tabix indexed if processing with proximal variants.'.format(self.input_file))
        if lib.utils.is_gz_file(self.input_file):
            mode = 'rb'
        else:
            mode = 'r'
        if self.proximal_variants_vcf:
            self.proximal_variants_tsv_fh = open(self.proximal_variants_tsv, 'w')
            self.proximal_variants_writer = csv.DictWriter(self.proximal_variants_tsv_fh, delimiter='\t', fieldnames=['chromosome_name', 'start', 'stop', 'reference', 'variant', 'amino_acid_change', 'codon_change', 'protein_position', 'type', 'main_somatic_variant'])
            self.proximal_variants_writer.writeheader()
            self.proximal_variant_parser = ProximalVariant(self.proximal_variants_vcf, self.pass_only)
            self.somatic_vcf_fh = open(self.input_file, mode)
            self.somatic_vcf_reader = vcf.Reader(self.somatic_vcf_fh)
        self.reader = open(self.input_file, mode)
        self.vcf_reader = vcf.Reader(self.reader)
        if len(self.vcf_reader.samples) > 1:
            if not self.sample_name:
                sys.exit("VCF contains more than one sample but sample_name is not set.")
            elif self.sample_name not in self.vcf_reader.samples:
                sys.exit("sample_name {} not in VCF {}".format(self.sample_name, self.input_file))
            if self.normal_sample_name is not None and self.normal_sample_name not in self.vcf_reader.samples:
                sys.exit("normal_sample_name {} not in VCF {}".format(self.normal_sample_name, self.input_file))
        elif len(self.vcf_reader.samples) ==  0:
            sys.exit("VCF doesn't contain any sample genotype information.")
        else:
            self.sample_name = self.vcf_reader.samples[0]
        self.writer = open(self.output_file, 'w')
        self.tsv_writer = csv.DictWriter(self.writer, delimiter='\t', fieldnames=self.output_headers(), restval='NA')
        self.tsv_writer.writeheader()
        self.csq_parser = self.create_csq_parser()

    def parse_bam_readcount_file(self, bam_readcount_file):
        with open(bam_readcount_file, 'r') as reader:
            coverage_tsv_reader = csv.reader(reader, delimiter='\t')
            coverage = {}
            for row in coverage_tsv_reader:
                chromosome     = row[0]
                position       = row[1]
                reference_base = row[2].upper()
                depth          = row[3]
                brct           = row[4:]
                if chromosome not in coverage:
                    coverage[chromosome] = {}
                if position not in coverage[chromosome]:
                    coverage[chromosome][position] = {}
                coverage[chromosome][position][reference_base] = self.parse_brct_field(brct)
                coverage[chromosome][position][reference_base]['depth'] = depth
        return coverage

    def parse_brct_field(self, brct_entry):
        parsed_brct = {}
        for brct in brct_entry:
            (base, count, rest) = brct.split(':', 2)
            parsed_brct[base.upper()] = count
        return parsed_brct

    def is_insertion(self, ref, alt):
        return len(alt) > len(ref)

    def is_deletion(self, ref, alt):
        return len(alt) < len(ref)

    def simplify_indel_allele(self, ref, alt):
        while len(ref)> 0 and len(alt) > 0 and ref[-1] == alt[-1]:
            ref = ref[0:-1]
            alt = alt[0:-1]
        while len(ref)> 0 and len(alt) > 0 and ref[0] == alt[0]:
            ref = ref[1:]
            alt = alt[1:]
        return ref, alt

    def create_csq_parser(self):
        info_fields = self.vcf_reader.infos

        if 'CSQ' not in info_fields:
            sys.exit('Input VCF does not contain a CSQ header. Please annotate the VCF with VEP before running it.')
        if info_fields['CSQ'] is None:
            sys.exit('Failed to extract format string from info description for tag (CSQ)')
        else:
            csq_header = info_fields['CSQ']
            return CsqParser(csq_header.desc)

    def resolve_consequence(self, consequence_string):
        if '&' in consequence_string:
            consequences = {consequence.lower() for consequence in consequence_string.split('&')}
        elif '.' in consequence_string:
            consequences = {consequence.lower() for consequence in consequence_string.split('.')}
        else:
            consequences = [consequence_string.lower()]

        if 'start_lost' in consequences:
            consequence = None
        elif 'frameshift_variant' in consequences:
            consequence = 'FS'
        elif 'missense_variant' in consequences:
            consequence = 'missense'
        elif 'inframe_insertion' in consequences:
            consequence = 'inframe_ins'
        elif 'inframe_deletion' in consequences:
            consequence = 'inframe_del'
        else:
            consequence = None
        return consequence

    def calculate_vaf(self, var_count, depth):
        if depth == 0:
            return 'NA'
        else:
            return (var_count / depth)

    def parse_gene_expns_file(self):
        gene_expns = {}
        if self.gene_expn_file is not None:
            with open(self.gene_expn_file, 'r') as reader:
                genes_tsv_reader = csv.DictReader(reader, delimiter='\t')
                for row in genes_tsv_reader:
                    if row['tracking_id'] not in gene_expns.keys():
                        gene_expns[row['tracking_id']] = {}
                    gene_expns[row['tracking_id']][row['locus']] = row
        return gene_expns

    def parse_transcript_expns_file(self):
        transcript_expns = {}
        if self.transcript_expn_file is not None:
            with open(self.transcript_expn_file, 'r') as reader:
                isoforms_tsv_reader = csv.DictReader(reader, delimiter='\t')
                for row in isoforms_tsv_reader:
                    transcript_expns[row['tracking_id']] = row
        return transcript_expns

    def parse_coverage_files(self):
        coverage = {}
        for variant_type in ['snvs', 'indels']:
            for data_type in ['normal', 'tdna', 'trna']:
                coverage_file_name = '_'.join([data_type, variant_type, 'coverage_file'])
                coverage_file = getattr(self, coverage_file_name)
                if coverage_file is not None:
                    if variant_type not in coverage:
                        coverage[variant_type] = {}
                    coverage[variant_type][data_type] = self.parse_bam_readcount_file(coverage_file)
        return coverage

    def determine_bam_readcount_bases(self, entry, reference, alt, start):
        if len(reference) == len(alt):
            bam_readcount_position = entry.POS
            variant_type = 'snvs'
            ref_base = reference
            var_base = alt
        else:
            if self.is_deletion(reference, alt):
                bam_readcount_position = start + 1
                (simplified_reference, simplified_alt) = self.simplify_indel_allele(reference, alt)
                ref_base = reference[1:2]
                var_base = '-' + simplified_reference
            elif self.is_insertion(reference, alt):
                bam_readcount_position = start
                (simplified_reference, simplified_alt) = self.simplify_indel_allele(reference, alt)
                ref_base = reference
                var_base = '+' + simplified_alt
            variant_type = 'indels'
        return (bam_readcount_position, ref_base, var_base, variant_type)

    def get_depth_from_vcf_genotype(self, genotype, tag):
        try:
            depth = genotype[tag]
            if depth is None or depth == "":
                depth = 'NA'
        except AttributeError:
            depth = 'NA'
        return depth

    def get_vaf_from_vcf_genotype(self, genotype, alts, alt, af_tag, ad_tag, dp_tag):
        try:
            allele_frequencies = genotype[af_tag]
            if isinstance(allele_frequencies, list):
                vaf = allele_frequencies[alts.index(alt)]
            else:
                vaf = allele_frequencies
            if vaf > 1:
                print("Warning: VAF is expected to be a fraction, but is larger than 1. If VAFs are encoded as percentages, please adjust the coverage cutoffs accordingly.")
        except (AttributeError, TypeError):
            try:
                allele_depths = genotype[ad_tag]
                if isinstance(allele_depths, list):
                    #sometimes AF is type R, sometimes it's A
                    if len(allele_depths) == len(alts):
                        var_count = allele_depths[alts.index(alt)]
                    elif len(allele_depths) == len(alts) + 1:
                        var_count = allele_depths[alts.index(alt) + 1]
                    else:
                        print("Warning: Mismatch between the number of alternate alleles and number of values in the AD field for genotype {}".format(genotype))
                        return 'NA'
                else:
                    var_count = allele_depths
                if var_count is None or var_count == "":
                    return 'NA'
                depth = genotype[dp_tag]
                if depth is None or depth == "":
                    return 'NA'
                vaf = self.calculate_vaf(int(var_count), int(depth))
            except AttributeError:
                vaf = 'NA'
        return vaf

    def calculate_coverage_for_entry(self, coverage, entry, reference, alt, start, chromosome, genotype):
        (bam_readcount_position, ref_base, var_base, variant_type) = self.determine_bam_readcount_bases(entry, reference, alt, start)
        coverage_for_entry = {}
        if variant_type in coverage:
            for coverage_type in coverage[variant_type]:
                if (
                    chromosome in coverage[variant_type][coverage_type]
                    and str(bam_readcount_position) in coverage[variant_type][coverage_type][chromosome]
                    and ref_base in coverage[variant_type][coverage_type][chromosome][str(bam_readcount_position)]
                ):
                    brct = coverage[variant_type][coverage_type][chromosome][str(bam_readcount_position)][ref_base]
                    if 'depth' in brct and var_base in brct:
                        coverage_for_entry[coverage_type + '_depth'] = int(brct['depth'])
                        coverage_for_entry[coverage_type + '_vaf']   = self.calculate_vaf(int(brct[var_base]), int(brct['depth']))
        else:
            coverage_for_entry['tdna_depth'] = self.get_depth_from_vcf_genotype(genotype, 'DP')
            coverage_for_entry['trna_depth'] = self.get_depth_from_vcf_genotype(genotype, 'RDP')
            coverage_for_entry['tdna_vaf'] = self.get_vaf_from_vcf_genotype(genotype, entry.ALT, alt, 'AF', 'AD', 'DP')
            coverage_for_entry['trna_vaf'] = self.get_vaf_from_vcf_genotype(genotype, entry.ALT, alt, 'RAF', 'RAD', 'RDP')
            if self.normal_sample_name is not None:
                normal_genotype = entry.genotype(self.normal_sample_name)
                coverage_for_entry['normal_depth'] = self.get_depth_from_vcf_genotype(normal_genotype, 'DP')
                coverage_for_entry['normal_vaf'] = self.get_vaf_from_vcf_genotype(normal_genotype, entry.ALT, alt, 'AF', 'AD', 'DP')
        return coverage_for_entry

    def write_proximal_variant_entries(self, entry, alt, transcript_name, index):
        proximal_variants = self.proximal_variant_parser.extract(entry, alt, transcript_name, self.peptide_length)
        for (proximal_variant, csq_entry) in proximal_variants:
            if len(list(self.somatic_vcf_reader.fetch(proximal_variant.CHROM, proximal_variant.POS - 1 , proximal_variant.POS))) > 0:
                proximal_variant_type = 'somatic'
            else:
                proximal_variant_type = 'germline'
            proximal_variant_entry = {
                'chromosome_name': proximal_variant.CHROM,
                'start': proximal_variant.affected_start,
                'stop': proximal_variant.affected_end,
                'reference': proximal_variant.REF,
                'variant': proximal_variant.ALT[0],
                'amino_acid_change': csq_entry['Amino_acids'],
                'codon_change': csq_entry['Codons'],
                'protein_position': csq_entry['Protein_position'],
                'type': proximal_variant_type,
                'main_somatic_variant': index,
            }
            self.proximal_variants_writer.writerow(proximal_variant_entry)

    def close_filehandles(self):
        self.writer.close()
        self.reader.close()
        if self.proximal_variants_vcf:
            self.proximal_variant_parser.fh.close()
            self.proximal_variants_tsv_fh.close()
            self.somatic_vcf_fh.close()

    def decode_hex(self, string):
        hex_string = string.group(0).replace('%', '')
        return binascii.unhexlify(hex_string).decode('utf-8')

    def execute(self):
        gene_expns = self.parse_gene_expns_file()
        transcript_expns = self.parse_transcript_expns_file()
        coverage = self.parse_coverage_files()

        indexes = []
        count = 1
        for entry in self.vcf_reader:
            chromosome = entry.CHROM
            start      = entry.affected_start
            stop       = entry.affected_end
            reference  = entry.REF
            alts       = entry.ALT

            genotype = entry.genotype(self.sample_name)
            if genotype.gt_type is None or genotype.gt_type == 0:
                #The genotype is uncalled or hom_ref
                continue

            filt = entry.FILTER
            if self.pass_only and not (filt is None or len(filt) == 0):
                continue

            if 'CSQ' not in entry.INFO:
                continue

            alleles_dict = self.csq_parser.resolve_alleles(entry)
            for alt in alts:
                alt = str(alt)
                if genotype.gt_bases and alt not in genotype.gt_bases.split('/'):
                    continue

                coverage_for_entry = self.calculate_coverage_for_entry(coverage, entry, reference, alt, start, chromosome, genotype)

                transcripts = self.csq_parser.parse_csq_entries_for_allele(entry.INFO['CSQ'], alt)
                if len(transcripts) == 0:
                    csq_allele = alleles_dict[alt]
                    transcripts = self.csq_parser.parse_csq_entries_for_allele(entry.INFO['CSQ'], csq_allele)
                if len(transcripts) == 0 and self.is_deletion(reference, alt):
                    transcripts = self.csq_parser.parse_csq_entries_for_allele(entry.INFO['CSQ'], 'deletion')

                for transcript in transcripts:
                    transcript_name = transcript['Feature']
                    consequence = self.resolve_consequence(transcript['Consequence'])
                    if consequence is None:
                        continue
                    elif consequence == 'FS':
                        if transcript['DownstreamProtein'] == '':
                            print("frameshift_variant transcript does not contain a DownstreamProtein sequence. Skipping.\n{} {} {} {} {}".format(entry.CHROM, entry.POS, entry.REF, alt, transcript['Feature']))
                            continue
                        else:
                            amino_acid_change_position = "%s%s/%s" % (transcript['Protein_position'], entry.REF, alt)
                    else:
                        if transcript['Amino_acids'] == '':
                            print("Transcript does not contain Amino_acids change information. Skipping.\n{} {} {} {} {}".format(entry.CHROM, entry.POS, entry.REF, alt, transcript['Feature']))
                            continue
                        else:
                            amino_acid_change_position = transcript['Protein_position'] + transcript['Amino_acids']
                    gene_name = transcript['SYMBOL']
                    index = '%s.%s.%s.%s.%s' % (count, gene_name, transcript_name, consequence, amino_acid_change_position)
                    if index in indexes:
                        sys.exit("Warning: TSV index already exists: {}".format(index))
                    else:
                        indexes.append(index)
                        count += 1

                    if self.proximal_variants_vcf:
                        self.write_proximal_variant_entries(entry, alt, transcript_name, index)

                    ensembl_gene_id = transcript['Gene']
                    hgvsc = re.sub(r'%[0-9|A-F][0-9|A-F]', self.decode_hex, transcript['HGVSc']) if 'HGVSc' in transcript else 'NA'
                    hgvsp = re.sub(r'%[0-9|A-F][0-9|A-F]', self.decode_hex, transcript['HGVSp']) if 'HGVSp' in transcript else 'NA'
                    if 'TSL' in transcript and transcript['TSL'] is not None and transcript['TSL'] != '':
                        tsl = transcript['TSL']
                    else:
                        tsl = 'NA'
                    output_row = {
                        'chromosome_name'                : entry.CHROM,
                        'start'                          : entry.affected_start,
                        'stop'                           : entry.affected_end,
                        'reference'                      : entry.REF,
                        'variant'                        : alt,
                        'gene_name'                      : gene_name,
                        'transcript_name'                : transcript_name,
                        'transcript_support_level'       : tsl,
                        'ensembl_gene_id'                : ensembl_gene_id,
                        'hgvsc'                          : hgvsc,
                        'hgvsp'                          : hgvsp,
                        'wildtype_amino_acid_sequence'   : transcript['WildtypeProtein'],
                        'downstream_amino_acid_sequence' : transcript['DownstreamProtein'],
                        'fusion_amino_acid_sequence'     : '',
                        'variant_type'                   : consequence,
                        'protein_position'               : transcript['Protein_position'],
                        'index'                          : index,
                        'protein_length_change'          : transcript['ProteinLengthChange'],
                    }
                    if transcript['Amino_acids']:
                        output_row['amino_acid_change'] = transcript['Amino_acids']

                    if transcript['Codons']:
                        output_row['codon_change'] =  transcript['Codons']
                    else:
                        output_row['codon_change'] = 'NA'

                    if transcript_name in transcript_expns.keys():
                        transcript_expn_entry = transcript_expns[transcript_name]
                        output_row['transcript_expression'] = transcript_expn_entry['FPKM']
                    elif 'TX' in self.vcf_reader.formats:
                        if 'TX' in genotype.data._asdict():
                            transcript_expressions = genotype['TX']
                            if isinstance(transcript_expressions, list):
                                for transcript_expression in transcript_expressions:
                                    (transcript, value) = transcript_expression.split('|')
                                    if transcript == transcript_name:
                                        output_row['transcript_expression'] = value
                            elif transcript_expressions is not None:
                                (transcript, value) = transcript_expressions.split('|')
                                if transcript == transcript_name:
                                    output_row['transcript_expression'] = value

                    if ensembl_gene_id in gene_expns.keys():
                        gene_expn_entries = gene_expns[ensembl_gene_id]
                        gene_fpkm = 0
                        for locus, gene_expn_entry in gene_expn_entries.items():
                            gene_fpkm += float(gene_expn_entry['FPKM'])
                        output_row['gene_expression'] = gene_fpkm
                    elif 'GX' in self.vcf_reader.formats:
                        if 'GX' in genotype.data._asdict():
                            gene_expressions = genotype['GX']
                            if isinstance(gene_expressions, list):
                                for gene_expression in gene_expressions:
                                    (gene, value) = gene_expression.split('|')
                                    if ensembl_gene_id == gene or gene_name == gene:
                                        output_row['gene_expression'] = value
                            elif gene_expressions is not None:
                                (gene, value) = gene_expressions.split('|')
                                if ensembl_gene_id == gene or gene_name == gene:
                                    output_row['gene_expression'] = value

                    output_row.update(coverage_for_entry)

                    self.tsv_writer.writerow(output_row)

        self.close_filehandles()

class IntegrateConverter(InputFileConverter):
    def input_fieldnames(self):
        return [
            'chr 5p',
            'start 5p',
            'end 5p',
            'chr 3p',
            'start 3p',
            'end 3p',
            'name of fusion',
            'tier of fusion',
            'strand 5p',
            'strand 3p',
            'quantitation',
            'is canonical boundary',
            'can be in-frame',
            'peptides',
            'fusion positions',
            'number of nucleotides in the break',
            'transcripts',
            'is canonical intron dinucleotide',
        ]

    def fusions_for_three_p_transcripts(self, five_p_transcript, three_p_transcripts):
        fusions = []
        if len(three_p_transcripts):
            for three_p_transcript in three_p_transcripts.split('|'):
                fusions.append("%s-%s"% (five_p_transcript, three_p_transcript))
        return fusions

    def execute(self):
        reader = open(self.input_file, 'r')
        csv_reader = csv.DictReader(reader, delimiter='\t', fieldnames=self.input_fieldnames())
        writer = open(self.output_file, 'w')
        tsv_writer = csv.DictWriter(writer, delimiter='\t', fieldnames=self.output_headers(), restval='NA')
        tsv_writer.writeheader()
        count = 1
        for entry in csv_reader:
            output_row = {
                'chromosome_name'            : "%s / %s" % (entry['chr 5p'], entry['chr 3p']),
                'start'                      : "%s / %s" % (entry['start 5p'], entry['start 3p']),
                'stop'                       : "%s / %s" % (entry['end 5p'], entry['end 3p']),
                'reference'                  : 'fusion',
                'variant'                    : 'fusion',
                'gene_name'                  : entry['name of fusion'],
                'wildtype_amino_acid_sequence'   : '',
                'downstream_amino_acid_sequence' : '',
                'protein_length_change'      : '',
                'amino_acid_change'          : 'NA',
                'codon_change'               : 'NA',
                'ensembl_gene_id'            : 'NA',
                'amino_acid_change'          : 'NA',
                'transcript_expression'      : 'NA',
                'gene_expression'            : 'NA',
                'normal_depth'               : 'NA',
                'normal_vaf'                 : 'NA',
                'tdna_depth'                 : 'NA',
                'tdna_vaf'                   : 'NA',
                'trna_depth'                 : 'NA',
                'trna_vaf'                   : 'NA',
            }

            if entry['fusion positions'] == 'NA' or entry['transcripts'] == 'NA' or entry['peptides'] == 'NA':
                continue
            for (fusion_position, transcript_set, fusion_amino_acid_sequence) in zip(entry['fusion positions'].split(','), entry['transcripts'].split(','), entry['peptides'].split(',')):
                (five_p_transcripts, three_p_inframe_transcripts, three_p_frameshift_transcripts) = transcript_set.split(';')
                fusions    = []
                for five_p_transcript in five_p_transcripts.split('|'):
                    fusions.extend(self.fusions_for_three_p_transcripts(five_p_transcript, three_p_inframe_transcripts))
                    fusions.extend(self.fusions_for_three_p_transcripts(five_p_transcript, three_p_frameshift_transcripts))

                if entry['can be in-frame']:
                    variant_type = 'inframe_fusion'
                else:
                    variant_type = 'frameshift_fusion'

                output_row['variant_type']               = variant_type
                output_row['protein_position']           = fusion_position
                output_row['fusion_amino_acid_sequence'] = fusion_amino_acid_sequence
                output_row['transcript_name']            = ';'.join(fusions)
                output_row['index']                      = '%s.%s.%s.%s' % (count, entry['name of fusion'], variant_type, fusion_position)
                tsv_writer.writerow(output_row)

                count += 1

        writer.close()
        reader.close()

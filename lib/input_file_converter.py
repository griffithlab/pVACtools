import vcf
import csv
import sys
import os
from abc import ABCMeta
from collections import OrderedDict
from lib.csq_parser import CsqParser
import lib.run_utils
from lib.proximal_variant import ProximalVariant
import binascii
import re
import glob
from Bio import SeqIO
import logging

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
            'frameshift_amino_acid_sequence',
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
        self.pass_only                   = kwargs.pop('pass_only', False)
        self.sample_name        = kwargs.pop('sample_name', None)
        self.normal_sample_name = kwargs.pop('normal_sample_name', None)
        self.proximal_variants_vcf = kwargs.pop('proximal_variants_vcf', None)
        self.proximal_variants_tsv = kwargs.pop('proximal_variants_tsv', None)
        self.flanking_bases = kwargs.pop('flanking_bases', None)
        if self.proximal_variants_vcf and not (self.proximal_variants_tsv and self.flanking_bases):
            sys.exit("A proximal variants TSV output path and number of flanking bases need to be specified if a proximal variants input VCF is provided.")
        if self.proximal_variants_vcf and not lib.run_utils.is_gz_file(self.input_file):
            sys.exit("Input VCF {} needs to be bgzipped when running with a proximal variants VCF.".format(self.input_file))
        if self.proximal_variants_vcf and not lib.run_utils.is_gz_file(self.proximal_variants_vcf):
            sys.exit("Proximal variants VCF {} needs to be bgzipped.".format(self.proximal_variants_vcf))
        if self.proximal_variants_vcf and not os.path.exists(self.proximal_variants_vcf + '.tbi'):
            sys.exit('No .tbi file found for proximal variants VCF {}. Proximal variants VCF needs to be tabix indexed.'.format(self.proximal_variants_vcf))
        if self.proximal_variants_vcf and not os.path.exists(self.input_file + '.tbi'):
            sys.exit('No .tbi file found for input VCF {}. Input VCF needs to be tabix indexed if processing with proximal variants.'.format(self.input_file))
        if lib.run_utils.is_gz_file(self.input_file):
            mode = 'rb'
        else:
            mode = 'r'
        if self.proximal_variants_vcf:
            self.proximal_variants_tsv_fh = open(self.proximal_variants_tsv, 'w')
            self.proximal_variants_writer = csv.DictWriter(self.proximal_variants_tsv_fh, delimiter='\t', fieldnames=['chromosome_name', 'start', 'stop', 'reference', 'variant', 'amino_acid_change', 'codon_change', 'protein_position', 'type', 'main_somatic_variant'])
            self.proximal_variants_writer.writeheader()
            self.proximal_variant_parser = ProximalVariant(self.proximal_variants_vcf, self.pass_only, self.flanking_bases)
            self.somatic_vcf_fh = open(self.input_file, mode)
            self.somatic_vcf_reader = vcf.Reader(self.somatic_vcf_fh)
        self.reader = open(self.input_file, mode)
        self.vcf_reader = vcf.Reader(self.reader)
        if len(self.vcf_reader.samples) > 1:
            if not self.sample_name:
                sys.exit("VCF contains more than one sample but sample_name is not set.")
            elif self.sample_name not in self.vcf_reader.samples:
                sys.exit("sample_name {} not a sample ID in the #CHROM header of VCF {}".format(self.sample_name, self.input_file))
            if self.normal_sample_name is not None and self.normal_sample_name not in self.vcf_reader.samples:
                sys.exit("normal_sample_name {} not a sample ID in the #CHROM header of VCF {}".format(self.normal_sample_name, self.input_file))
        elif len(self.vcf_reader.samples) ==  0:
            sys.exit("VCF doesn't contain any sample genotype information. Add a dummy sample using the vcf-genotype-annotator tool available as part of the vatools package.")
        else:
            if self.normal_sample_name is not None:
                sys.exit("normal_sample_name {} provided but the input file is a single-sample (tumor only) VCF".format(self.normal_sample_name))
            self.sample_name = self.vcf_reader.samples[0]
        if 'GT' not in self.vcf_reader.formats:
            sys.exit("VCF doesn't contain any sample genotype information. Add a dummy sample using the vcf-genotype-annotator tool available as part of the vatools package.")
        self.writer = open(self.output_file, 'w')
        self.tsv_writer = csv.DictWriter(self.writer, delimiter='\t', fieldnames=self.output_headers(), restval='NA')
        self.tsv_writer.writeheader()
        self.csq_parser = self.create_csq_parser()
        if 'FrameshiftSequence' not in self.csq_parser.csq_format:
            sys.exit("VCF doesn't contain VEP FrameshiftSequence annotations. Please re-annotate the VCF with VEP and the Wildtype and Frameshift plugins.")
        if 'WildtypeProtein' not in self.csq_parser.csq_format:
            sys.exit("VCF doesn't contain VEP WildtypeProtein annotations. Please re-annotate the VCF with VEP and the Wildtype and Frameshift plugins.")

    def is_insertion(self, ref, alt):
        return len(alt) > len(ref)

    def is_deletion(self, ref, alt):
        return len(alt) < len(ref)

    def create_csq_parser(self):
        info_fields = self.vcf_reader.infos

        if 'CSQ' not in info_fields:
            sys.exit('Input VCF does not contain a CSQ header. Please annotate the VCF with VEP before running it.')
        if info_fields['CSQ'] is None:
            sys.exit('Failed to extract format string from info description for tag (CSQ)')
        else:
            csq_header = info_fields['CSQ']
            return CsqParser(csq_header.desc)

    def resolve_consequence(self, consequence_string, ref, alt):
        if '&' in consequence_string:
            consequences = {consequence.lower() for consequence in consequence_string.split('&')}
        elif '.' in consequence_string:
            consequences = {consequence.lower() for consequence in consequence_string.split('.')}
        else:
            consequences = [consequence_string.lower()]

        if 'start_lost' in consequences:
            consequence = None
        elif 'stop_retained_variant' in consequences:
            consequence = None
        elif 'frameshift_variant' in consequences:
            consequence = 'FS'
        elif 'missense_variant' in consequences:
            consequence = 'missense'
        elif 'inframe_insertion' in consequences:
            consequence = 'inframe_ins'
        elif 'inframe_deletion' in consequences:
            consequence = 'inframe_del'
        elif 'protein_altering_variant' in consequences:
            if len(ref) > len(alt) and (len(ref) - len(alt)) % 3 == 0:
                consequence = 'inframe_del'
            elif len(alt) > len(ref) and (len(alt) - len(ref)) % 3 == 0:
                consequence = 'inframe_ins'
            else:
                consequence = None
        else:
            consequence = None
        return consequence

    def calculate_vaf(self, var_count, depth):
        if depth == 0:
            return 'NA'
        else:
            return (var_count / depth)

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

    def calculate_coverage_for_entry(self, entry, reference, alt, start, chromosome, genotype):
        coverage_for_entry = {}
        coverage_for_entry['tdna_depth'] = self.get_depth_from_vcf_genotype(genotype, 'DP')
        coverage_for_entry['trna_depth'] = self.get_depth_from_vcf_genotype(genotype, 'RDP')
        alts = list(map(lambda x: str(x) , entry.ALT))
        coverage_for_entry['tdna_vaf'] = self.get_vaf_from_vcf_genotype(genotype, alts, alt, 'AF', 'AD', 'DP')
        coverage_for_entry['trna_vaf'] = self.get_vaf_from_vcf_genotype(genotype, alts, alt, 'RAF', 'RAD', 'RDP')
        if self.normal_sample_name is not None:
            normal_genotype = entry.genotype(self.normal_sample_name)
            coverage_for_entry['normal_depth'] = self.get_depth_from_vcf_genotype(normal_genotype, 'DP')
            coverage_for_entry['normal_vaf'] = self.get_vaf_from_vcf_genotype(normal_genotype, alts, alt, 'AF', 'AD', 'DP')
        return coverage_for_entry

    def write_proximal_variant_entries(self, entry, alt, transcript_name, index):
        proximal_variants = self.proximal_variant_parser.extract(entry, alt, transcript_name)
        for (proximal_variant, csq_entry) in proximal_variants:
            if len(list(self.somatic_vcf_reader.fetch(proximal_variant.CHROM, proximal_variant.POS - 1 , proximal_variant.POS))) > 0:
                proximal_variant_type = 'somatic'
            else:
                proximal_variant_type = 'germline'
            if '/' in csq_entry['Protein_position']:
                protein_position = csq_entry['Protein_position'].split('/')[0]
                if protein_position == '-':
                    protein_position = csq_entry['Protein_position'].split('/')[1]
            else:
                protein_position = csq_entry['Protein_position']
            if protein_position == '-' or protein_position == '':
                print("Proximal variant doesn't have protein position information. Skipping.\n{} {} {} {} {}".format(entry.CHROM, entry.POS, entry.REF, alt, csq_entry['Feature']))
                continue
            proximal_variant_entry = {
                'chromosome_name': proximal_variant.CHROM,
                'start': proximal_variant.affected_start,
                'stop': proximal_variant.affected_end,
                'reference': proximal_variant.REF,
                'variant': proximal_variant.ALT[0],
                'amino_acid_change': csq_entry['Amino_acids'],
                'codon_change': csq_entry['Codons'],
                'protein_position': protein_position,
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
        indexes = []
        count = 1
        while True:
            try:
                entry = next(self.vcf_reader)
            except StopIteration:
                break
            except ValueError as e:
                raise Exception("VCF is truncated in the middle of an entry near string '{}'".format(str(e).split("'")[1]))
            except IndexError as e:
                raise Exception("VCF is truncated at the end of the file")
            except Exception as e:
                raise Exception("Error while reading VCF entry: {}".format(str(e)))
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

                coverage_for_entry = self.calculate_coverage_for_entry(entry, reference, alt, start, chromosome, genotype)

                transcripts = self.csq_parser.parse_csq_entries_for_allele(entry.INFO['CSQ'], alt)
                if len(transcripts) == 0:
                    csq_allele = alleles_dict[alt]
                    transcripts = self.csq_parser.parse_csq_entries_for_allele(entry.INFO['CSQ'], csq_allele)
                if len(transcripts) == 0 and self.is_deletion(reference, alt):
                    transcripts = self.csq_parser.parse_csq_entries_for_allele(entry.INFO['CSQ'], 'deletion')

                for transcript in transcripts:
                    if '/' in transcript['Protein_position']:
                        protein_position = transcript['Protein_position'].split('/')[0]
                        if protein_position == '-':
                            protein_position = transcript['Protein_position'].split('/')[1]
                    else:
                        protein_position = transcript['Protein_position']
                    if protein_position == '-' or protein_position == '':
                        print("Variant doesn't have protein position information. Skipping.\n{} {} {} {} {}".format(entry.CHROM, entry.POS, entry.REF, alt, transcript['Feature']))
                        continue
                    transcript_name = transcript['Feature']
                    consequence = self.resolve_consequence(transcript['Consequence'], reference, alt)
                    if consequence is None:
                        continue
                    elif consequence == 'FS':
                        if transcript['FrameshiftSequence'] == '':
                            print("frameshift_variant transcript does not contain a FrameshiftSequence. Skipping.\n{} {} {} {} {}".format(entry.CHROM, entry.POS, entry.REF, alt, transcript['Feature']))
                            continue
                        else:
                            amino_acid_change_position = "%s%s/%s" % (protein_position, entry.REF, alt)
                    else:
                        if transcript['Amino_acids'] == '':
                            print("Transcript does not contain Amino_acids change information. Skipping.\n{} {} {} {} {}".format(entry.CHROM, entry.POS, entry.REF, alt, transcript['Feature']))
                            continue
                        else:
                            amino_acid_change_position = protein_position + transcript['Amino_acids']
                    gene_name = transcript['SYMBOL']
                    index = lib.run_utils.construct_index(count, gene_name, transcript_name, consequence, amino_acid_change_position)
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

                    wildtype_amino_acid_sequence = transcript['WildtypeProtein']
                    if '*' in wildtype_amino_acid_sequence:
                        logging.warning("Transcript WildtypeProtein sequence contains internal stop codon. These can occur in Ensembl transcripts of the biotype polymorphic_pseudogene. Skipping.\n{} {} {} {} {}".format(entry.CHROM, entry.POS, entry.REF, alt, transcript['Feature']))
                        continue

                    if transcript['FrameshiftSequence'] != '':
                        wt_len = len(wildtype_amino_acid_sequence)
                        mt_len = len(transcript['FrameshiftSequence'])
                        if mt_len > wt_len:
                            protein_length_change = "+{}".format(mt_len - wt_len)
                        elif mt_len < wt_len:
                            protein_length_change = "-{}".format(wt_len - mt_len)
                        else:
                            protein_length_change = "0"
                    else:
                        protein_length_change = ""

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
                        'wildtype_amino_acid_sequence'   : wildtype_amino_acid_sequence,
                        'frameshift_amino_acid_sequence' : transcript['FrameshiftSequence'],
                        'fusion_amino_acid_sequence'     : '',
                        'variant_type'                   : consequence,
                        'protein_position'               : protein_position,
                        'index'                          : index,
                        'protein_length_change'          : protein_length_change,
                    }
                    if transcript['Amino_acids']:
                        output_row['amino_acid_change'] = transcript['Amino_acids']

                    if transcript['Codons']:
                        output_row['codon_change'] =  transcript['Codons']
                    else:
                        output_row['codon_change'] = 'NA'

                    for (tag, key, comparison_fields) in zip(['TX', 'GX'], ['transcript_expression', 'gene_expression'], [[transcript_name], [ensembl_gene_id, gene_name]]):
                        if tag in self.vcf_reader.formats:
                            if tag in genotype.data._asdict():
                                expressions = genotype[tag]
                                if isinstance(expressions, list):
                                    for expression in expressions:
                                        (item, value) = expression.split('|')
                                        for comparison_field in comparison_fields:
                                            if item == comparison_field:
                                                output_row[key] = value
                                elif expressions is not None:
                                    (item, value) = expressions.split('|')
                                    for comparison_field in comparison_fields:
                                        if item == comparison_field:
                                            output_row[key] = value

                    output_row.update(coverage_for_entry)

                    self.tsv_writer.writerow(output_row)

        self.close_filehandles()

class FusionInputConverter(InputFileConverter):
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

    def parse_integrate_neo_file(self):
        reader = open(self.input_file, 'r')
        csv_reader = csv.DictReader(reader, delimiter='\t', fieldnames=self.input_fieldnames())
        count = 1
        output_rows = []
        for entry in csv_reader:
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

                output_row = {
                    'chromosome_name'            : "%s / %s" % (entry['chr 5p'], entry['chr 3p']),
                    'start'                      : "%s / %s" % (entry['start 5p'], entry['start 3p']),
                    'stop'                       : "%s / %s" % (entry['end 5p'], entry['end 3p']),
                    'reference'                  : 'fusion',
                    'variant'                    : 'fusion',
                    'gene_name'                  : entry['name of fusion'],
                    'wildtype_amino_acid_sequence'   : '',
                    'frameshift_amino_acid_sequence' : '',
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
                output_row['variant_type']               = variant_type
                output_row['protein_position']           = fusion_position
                output_row['fusion_amino_acid_sequence'] = fusion_amino_acid_sequence
                transcripts                              = ';'.join(fusions)
                output_row['transcript_name']            = transcripts
                output_row['index']                      = lib.run_utils.construct_index(count, entry['name of fusion'], transcripts, variant_type, fusion_position)
                output_rows.append(output_row)
                count += 1

        reader.close()
        return output_rows

    def determine_fusion_sequence(self, full_sequence, variant_type):
        if '*' not in full_sequence:
            sys.exit("Fusion position marker '*' not found in fusion sequence. Please rerun AGFusion using the `--middlestar` option.")
        else:
            fusion_position = full_sequence.find('*')
            sequence = full_sequence.replace('*', '')
            return (fusion_position, sequence)

    def parse_exon_file(self, input_file):
        exon_file = input_file.replace('_protein.fa', '.exons.txt')
        if not os.path.exists(exon_file):
            exon_file = exon_file.replace('.txt', '.csv')
        five_prime_positions, three_prime_positions = [], []
        with open(exon_file, 'r') as fh:
            dialect = csv.Sniffer().sniff(fh.read())
            fh.seek(0)
            reader = csv.DictReader(fh, delimiter=dialect.delimiter)
            for record in reader:
                exon_start = int(record['exon_start'])
                exon_end = int(record['exon_end'])
                if record['exon_gene_source'] == "'5 gene":
                    five_prime_chr = record['exon_chr']
                    five_prime_positions.append(exon_start)
                    five_prime_positions.append(exon_end)
                else:
                    three_prime_chr = record['exon_chr']
                    three_prime_positions.append(exon_start)
                    three_prime_positions.append(exon_end)
        five_prime_start = min(five_prime_positions)
        five_prime_end = max(five_prime_positions)
        three_prime_start = min(three_prime_positions)
        three_prime_end = max(three_prime_positions)
        return (five_prime_chr, five_prime_start, five_prime_end, three_prime_chr, three_prime_start, three_prime_end)

    def parse_agfusion_files(self):
        output_rows = []
        count = 1
        for input_file in sorted(glob.glob(os.path.join(self.input_file, '*', '*_protein.fa'))):
            (five_prime_chr, five_prime_start, five_prime_end, three_prime_chr, three_prime_start, three_prime_end) = self.parse_exon_file(input_file)
            for record in SeqIO.parse(input_file, "fasta"):
                record_info = dict(map(lambda x: x.split(': '), record.description.split(', ')[1:]))
                if record_info['effect'] == 'in-frame' or record_info['effect'] == 'in-frame (with mutation)':
                    variant_type = 'inframe_fusion'
                elif record_info['effect'] == 'out-of-frame':
                    variant_type = 'frameshift_fusion'
                else:
                    sys.exit('Effect "{}" not supported'.format(record_info['effect']))

                (fusion_position, fusion_amino_acid_sequence) = self.determine_fusion_sequence(str(record.seq), variant_type)
                output_row = {
                    'chromosome_name'            : "{} / {}".format(five_prime_chr, three_prime_chr),
                    'start'                      : "{} / {}".format(five_prime_start, three_prime_start),
                    'stop'                       : "{} / {}".format(five_prime_end, three_prime_end),
                    'reference'                  : 'fusion',
                    'variant'                    : 'fusion',
                    'gene_name'                  : record_info['genes'],
                    'wildtype_amino_acid_sequence'   : '',
                    'frameshift_amino_acid_sequence' : '',
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
                    'variant_type'               : variant_type,
                    'protein_position'           : fusion_position,
                    'fusion_amino_acid_sequence' : fusion_amino_acid_sequence,
                    'transcript_name'            : record_info['transcripts'],
                    'index'                      : lib.run_utils.construct_index(count, record_info['genes'], record_info['transcripts'], variant_type, fusion_position)
                }
                output_rows.append(output_row)
                count += 1
        return output_rows

    def execute(self):
        writer = open(self.output_file, 'w')
        tsv_writer = csv.DictWriter(writer, delimiter='\t', fieldnames=self.output_headers(), restval='NA')
        tsv_writer.writeheader()
        if os.path.isfile(self.input_file):
            output_rows = self.parse_integrate_neo_file()
        elif os.path.isdir(self.input_file):
            output_rows = self.parse_agfusion_files()
        tsv_writer.writerows(output_rows)
        writer.close()

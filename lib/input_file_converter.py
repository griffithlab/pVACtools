import vcf
import csv
import sys
import re
from abc import ABCMeta
from collections import OrderedDict

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
            'amino_acid_change',
            'ensembl_gene_id',
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
                coverage[chromosome][position][reference_base] = brct
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

    def parse_csq_format(self, vcf_reader):
        info_fields = vcf_reader.infos

        if 'CSQ' not in info_fields:
            sys.exit('Input VCF does not contain a CSQ header. Please annotate the VCF with VEP before running it.')
        if info_fields['CSQ'] is None:
            sys.exit('Failed to extract format string from info description for tag (CSQ)')
        else:
            csq_header = info_fields['CSQ']
            format_pattern = re.compile('Format: (.*)')
            match = format_pattern.search(csq_header.desc)
            return match.group(1)

    def resolve_alleles(self, entry):
        alleles = {}
        if entry.is_indel:
            for alt in entry.ALT:
                alt = str(alt)
                if alt[0:1] != entry.REF[0:1]:
                    alleles[alt] = alt
                elif alt[1:] == "":
                    alleles[alt] = '-'
                else:
                    alleles[alt] = alt[1:]
        else:
            for alt in entry.ALT:
                alt = str(alt)
                alleles[alt] = alt

        return alleles

    def parse_csq_entries_for_allele(self, csq_entries, csq_format, csq_allele):
        csq_format_array = csq_format.split('|')

        transcripts = []
        for entry in csq_entries:
            values = entry.split('|')
            transcript = {}
            for key, value in zip(csq_format_array, values):
                transcript[key] = value
            if transcript['Allele'] == csq_allele:
                transcripts.append(transcript)

        return transcripts

    def resolve_consequence(self, consequence_string):
        consequences = {consequence.lower() for consequence in consequence_string.split('&')}
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

    def calculate_coverage(self, ref, var):
        return ref + var

    def calculate_vaf(self, ref, var):
        return (var / (self.calculate_coverage(ref, var)+0.00001)) * 100

    def execute(self):
        gene_expns = {}
        if self.gene_expn_file is not None:
            with open(self.gene_expn_file, 'r') as reader:
                genes_tsv_reader = csv.DictReader(reader, delimiter='\t')
                for row in genes_tsv_reader:
                    if row['tracking_id'] not in gene_expns.keys():
                        gene_expns[row['tracking_id']] = {}
                    gene_expns[row['tracking_id']][row['locus']] = row

        transcript_expns = {}
        if self.transcript_expn_file is not None:
            with open(self.transcript_expn_file, 'r') as reader:
                isoforms_tsv_reader = csv.DictReader(reader, delimiter='\t')
                for row in isoforms_tsv_reader:
                    transcript_expns[row['tracking_id']] = row

        coverage = {}
        for variant_type in ['snvs', 'indels']:
            for data_type in ['normal', 'tdna', 'trna']:
                coverage_file_name = '_'.join([data_type, variant_type, 'coverage_file'])
                coverage_file = getattr(self, coverage_file_name)
                if coverage_file is not None:
                    if variant_type not in coverage:
                        coverage[variant_type] = {}
                    coverage[variant_type][data_type] = self.parse_bam_readcount_file(coverage_file)

        reader = open(self.input_file, 'r')
        vcf_reader = vcf.Reader(reader)
        if len(vcf_reader.samples) > 1:
            sys.exit('ERROR: VCF file contains more than one sample')
        writer = open(self.output_file, 'w')
        tsv_writer = csv.DictWriter(writer, delimiter='\t', fieldnames=self.output_headers())
        tsv_writer.writeheader()

        csq_format = self.parse_csq_format(vcf_reader)
        for entry in vcf_reader:
            chromosome = entry.CHROM
            start      = entry.affected_start
            stop       = entry.affected_end
            reference  = entry.REF
            alts       = entry.ALT

            genotype = entry.genotype(vcf_reader.samples[0])
            if genotype.gt_type is None or genotype.gt_type == 0:
                #The genotype is uncalled or hom_ref
                continue

            if len(vcf_reader.samples) == 1:
                genotype = entry.genotype(vcf_reader.samples[0])
                if genotype.gt_type is None or genotype.gt_type == 0:
                    #The genotype is uncalled or hom_ref
                    continue

            alleles_dict = self.resolve_alleles(entry)
            indexes = []
            for alt in alts:
                alt = str(alt)
                if genotype.gt_bases and alt not in genotype.gt_bases.split('/'):
                    continue

                if entry.is_indel:
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
                else:
                    bam_readcount_position = entry.POS
                    variant_type = 'snvs'
                    ref_base = reference
                    var_base = alt
                coverage_for_entry = {}
                for coverage_type in ['normal', 'tdna', 'trna']:
                    coverage_for_entry[coverage_type + '_depth'] = 'NA'
                    coverage_for_entry[coverage_type + '_vaf'] = 'NA'
                if variant_type in coverage:
                    for coverage_type in coverage[variant_type]:
                        if (
                            chromosome in coverage[variant_type][coverage_type]
                            and str(bam_readcount_position) in coverage[variant_type][coverage_type][chromosome]
                            and ref_base in coverage[variant_type][coverage_type][chromosome][str(bam_readcount_position)]
                        ):
                            brct = self.parse_brct_field(coverage[variant_type][coverage_type][chromosome][str(bam_readcount_position)][ref_base])
                            if ref_base in brct and var_base in brct:
                                coverage_for_entry[coverage_type + '_depth'] = self.calculate_coverage(int(brct[ref_base]), int(brct[var_base]))
                                coverage_for_entry[coverage_type + '_vaf']   = self.calculate_vaf(int(brct[ref_base]), int(brct[var_base]))

                transcripts = self.parse_csq_entries_for_allele(entry.INFO['CSQ'], csq_format, alt)
                if len(transcripts) == 0:
                    csq_allele = alleles_dict[alt]
                    transcripts = self.parse_csq_entries_for_allele(entry.INFO['CSQ'], csq_format, csq_allele)

                for transcript in transcripts:
                    transcript_name = transcript['Feature']
                    consequence = self.resolve_consequence(transcript['Consequence'])
                    if consequence is None:
                        continue
                    elif consequence == 'FS':
                        if transcript['DownstreamProtein'] == '':
                            continue
                        else:
                            amino_acid_change_position = "%s%s/%s" % (transcript['Protein_position'], entry.REF, alt)
                    else:
                        amino_acid_change_position = transcript['Protein_position'] + transcript['Amino_acids']
                    gene_name = transcript['SYMBOL']
                    index = '%s.%s.%s.%s' % (gene_name, transcript_name, consequence, amino_acid_change_position)
                    if index in indexes:
                        sys.exit("Warning: TSV index already exists: {}".format(index))
                    else:
                        indexes.append(index)
                    ensembl_gene_id = transcript['Gene']
                    output_row = {
                        'chromosome_name'                : entry.CHROM,
                        'start'                          : entry.affected_start,
                        'stop'                           : entry.affected_end,
                        'reference'                      : entry.REF,
                        'variant'                        : alt,
                        'gene_name'                      : gene_name,
                        'transcript_name'                : transcript_name,
                        'ensembl_gene_id'                : ensembl_gene_id,
                        'wildtype_amino_acid_sequence'   : transcript['WildtypeProtein'],
                        'downstream_amino_acid_sequence' : transcript['DownstreamProtein'],
                        'variant_type'                   : consequence,
                        'protein_position'               : transcript['Protein_position'],
                        'index'                          : index,
                        'protein_length_change'          : transcript['ProteinLengthChange'],
                    }
                    if transcript['Amino_acids']:
                        output_row['amino_acid_change'] = transcript['Amino_acids']
                    else:
                        output_row['amino_acid_change'] = 'NA'

                    if transcript_name in transcript_expns.keys():
                        transcript_expn_entry = transcript_expns[transcript_name]
                        output_row['transcript_expression'] = transcript_expn_entry['FPKM']
                    else:
                        output_row['transcript_expression'] = 'NA'
                    if ensembl_gene_id in gene_expns.keys():
                        gene_expn_entries = gene_expns[ensembl_gene_id]
                        gene_fpkm = 0
                        for locus, gene_expn_entry in gene_expn_entries.items():
                            gene_fpkm += float(gene_expn_entry['FPKM'])
                        output_row['gene_expression'] = gene_fpkm
                    else:
                        output_row['gene_expression'] = 'NA'

                    output_row.update(coverage_for_entry)

                    tsv_writer.writerow(output_row)

        writer.close()
        reader.close()

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
        tsv_writer = csv.DictWriter(writer, delimiter='\t', fieldnames=self.output_headers())
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
                'amino_acid_change'          : 'NA',
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

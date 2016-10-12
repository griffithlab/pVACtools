import argparse
import vcf
import csv
import sys
import re

def parse_bam_readcount_file(bam_readcount_file):
    coverage_tsv_reader = csv.reader(bam_readcount_file, delimiter='\t')
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

def parse_brct_field(brct_entry):
    parsed_brct = {}
    for brct in brct_entry:
        (base, count, rest) = brct.split(':', 2)
        parsed_brct[base.upper()] = count
    return parsed_brct

def is_insertion(ref, alt):
    return len(alt) > len(ref)

def is_deletion(ref, alt):
    return len(alt) < len(ref)

def simplify_indel_allele(ref, alt):
    while len(ref)> 0 and len(alt) > 0 and ref[-1] == alt[-1]:
        ref = ref[0:-1]
        alt = alt[0:-1]
    while len(ref)> 0 and len(alt) > 0 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
    return ref, alt

def parse_csq_format(vcf_reader):
    info_fields = vcf_reader.infos

    if info_fields['CSQ'] is None:
        sys.exit('Failed to extract format string from info description for tag (CSQ)')
    else:
        csq_header = info_fields['CSQ']
        format_pattern = re.compile('Format: (.*)')
        match = format_pattern.search(csq_header.desc)
        return match.group(1)

def resolve_alleles(entry):
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

def parse_csq_entries_for_allele(csq_entries, csq_format, csq_allele):
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

def resolve_consequence(consequence_string):
    consequences = {consequence.lower() for consequence in consequence_string.split('&')}
    if 'frameshift_variant' in consequences:
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

def calculate_coverage(ref, var):
    return ref + var

def calculate_vaf(ref, var):
    return (var / (calculate_coverage(ref, var)+0.00001)) * 100

def output_headers():
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
        'index'
    ]

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser('pvacseq convert_vcf')
    parser.add_argument('input_file', type=argparse.FileType('r'), help='input VCF',)
    parser.add_argument('output_file', type=argparse.FileType('w'), help='output list of variants')
    parser.add_argument('-g', '--gene-expn-file', type=argparse.FileType('r'), help='genes.fpkm_tracking file from Cufflinks')
    parser.add_argument('-i', '--transcript-expn-file', type=argparse.FileType('r'), help='isoforms.fpkm_tracking file from Cufflinks')
    parser.add_argument('--normal-snvs-coverage-file', type=argparse.FileType('r'), help='bam-readcount output file for normal BAM and snvs'),
    parser.add_argument('--normal-indels-coverage-file', type=argparse.FileType('r'), help='bam-readcount output file for normal BAM and indels'),
    parser.add_argument('--tdna-snvs-coverage-file', type=argparse.FileType('r'), help='bam-readcount output file for tumor DNA BAM and snvs'),
    parser.add_argument('--tdna-indels-coverage-file', type=argparse.FileType('r'), help='bam-readcount output file for tumor DNA BAM and indels'),
    parser.add_argument('--trna-snvs-coverage-file', type=argparse.FileType('r'), help='bam-readcount output file for tumor RNA BAM and snvs'),
    parser.add_argument('--trna-indels-coverage-file', type=argparse.FileType('r'), help='bam-readcount output file for tumor RNA BAM and indels'),
    args = parser.parse_args(args_input)

    gene_expns = {}
    if args.gene_expn_file is not None:
        genes_tsv_reader = csv.DictReader(args.gene_expn_file, delimiter='\t')
        for row in genes_tsv_reader:
            if row['tracking_id'] not in gene_expns.keys():
                gene_expns[row['tracking_id']] = {}
            gene_expns[row['tracking_id']][row['locus']] = row
        args.gene_expn_file.close()

    transcript_expns = {}
    if args.transcript_expn_file is not None:
        isoforms_tsv_reader = csv.DictReader(args.transcript_expn_file, delimiter='\t')
        for row in isoforms_tsv_reader:
            transcript_expns[row['tracking_id']] = row
        args.transcript_expn_file.close()

    coverage = {}
    for variant_type in ['snvs', 'indels']:
        for data_type in ['normal', 'tdna', 'trna']:
            coverage_file_name = '_'.join([data_type, variant_type, 'coverage_file'])
            coverage_file = getattr(args, coverage_file_name)
            if coverage_file is not None:
                if variant_type not in coverage:
                    coverage[variant_type] = {}
                coverage[variant_type][data_type] = parse_bam_readcount_file(coverage_file)
                coverage_file.close()

    vcf_reader = vcf.Reader(args.input_file)
    if len(vcf_reader.samples) > 1:
        sys.exit('ERROR: VCF file contains more than one sample')
    tsv_writer = csv.DictWriter(args.output_file, delimiter='\t', fieldnames=output_headers())
    tsv_writer.writeheader()
    csq_format = parse_csq_format(vcf_reader)
    transcript_count = {}
    for entry in vcf_reader:
        chromosome = entry.CHROM
        start      = entry.affected_start
        stop       = entry.affected_end
        reference  = entry.REF
        alts       = entry.ALT

        alleles_dict = resolve_alleles(entry)
        for alt in alts:
            alt = str(alt)
            if entry.is_indel:
                if is_deletion(reference, alt):
                    bam_readcount_position = start + 1
                    (simplified_reference, simplified_alt) = simplify_indel_allele(reference, alt)
                    ref_base = reference[1:2]
                    var_base = '-' + simplified_reference
                elif is_insertion(reference, alt):
                    bam_readcount_position = start
                    (simplified_reference, simplified_alt) = simplify_indel_allele(reference, alt)
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
                        brct = parse_brct_field(coverage[variant_type][coverage_type][chromosome][str(bam_readcount_position)][ref_base])
                        if ref_base in brct and var_base in brct:
                            coverage_for_entry[coverage_type + '_depth'] = calculate_coverage(int(brct[ref_base]), int(brct[var_base]))
                            coverage_for_entry[coverage_type + '_vaf']   = calculate_vaf(int(brct[ref_base]), int(brct[var_base]))

            csq_allele = alleles_dict[alt]
            transcripts = parse_csq_entries_for_allele(entry.INFO['CSQ'], csq_format, csq_allele)
            for transcript in transcripts:
                transcript_name = transcript['Feature']
                if transcript_name in transcript_count:
                    transcript_count[transcript_name] += 1
                else:
                    transcript_count[transcript_name] = 1
                consequence = resolve_consequence(transcript['Consequence'])
                if consequence is None:
                    continue
                elif consequence == 'FS':
                    amino_acid_change_position = transcript['Protein_position']
                else:
                    amino_acid_change_position = transcript['Protein_position'] + transcript['Amino_acids']
                gene_name = transcript['SYMBOL']
                index = '%s_%s_%s.%s.%s' % (gene_name, transcript_name, transcript_count[transcript_name], consequence, amino_acid_change_position)
                ensembl_gene_id = transcript['Gene']
                output_row = {
                    'chromosome_name'                : entry.CHROM,
                    'start'                          : entry.affected_start,
                    'stop'                           : entry.affected_end,
                    'reference'                      : entry.REF,
                    'variant'                        : alt,
                    'gene_name'                      : gene_name,
                    'transcript_name'                : transcript_name,
                    'amino_acid_change'              : transcript['Amino_acids'],
                    'ensembl_gene_id'                : ensembl_gene_id,
                    'wildtype_amino_acid_sequence'   : transcript['WildtypeProtein'],
                    'downstream_amino_acid_sequence' : transcript['DownstreamProtein'],
                    'variant_type'                   : consequence,
                    'protein_position'               : transcript['Protein_position'],
                    'index'                          : index
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

    args.input_file.close()
    args.output_file.close()

if __name__ == '__main__':
    main()

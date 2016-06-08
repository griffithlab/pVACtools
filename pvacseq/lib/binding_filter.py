import argparse
import sys
import re
import os
import csv


def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser("Binding Filter")
    parser.add_argument('input', type=argparse.FileType('r'),
                        help="Input list of variants")
    parser.add_argument('fof', type=argparse.FileType('r'),
                        help="FOF containing list of parsed epitope files " +
                        "for different allele-length combinations (same sample)")
    parser.add_argument('output', type=argparse.FileType('w'),
                        help="Output .xls file containing list of filtered " +
                        "epitopes based on binding affinity for each " +
                        "allele-length combination per gene")
    parser.add_argument('-c', '--fold-change', type=int,
                        help="Minimum fold change between mutant binding " +
                        "score and wild-type score. The default is 0, which " +
                        "filters no results, but 1 is often a sensible " +
                        "default (requiring that binding is better to the MT " +
                        "than WT)",
                        default=0,
                        dest="minimum_fold_change")
    parser.add_argument('-b', '--binding-threshold', type=int,
                        help="Report only epitopes where the mutant allele " +
                        "has ic50 binding scores below this value; default 500",
                        default=500,
                        dest="binding_threshold")

    args = parser.parse_args(args_input)

    #precompile regex patterns used later
    chromosome_name = re.compile(r'^chromosome_name')
    netmhc_subber = re.compile(r"_netmhc")

    #open the variants file and parse into variants dictionary
    reader = csv.reader(args.input, delimiter='\t')
    variants = {}

    for data in reader:
        if chromosome_name.match(data[0]):
            variants['header'] = [word for word in data]
            continue

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

        gene = data[5]
        amino_acid_change = data[7]
        variants[ gene + '\t' + amino_acid_change] = [word for word in data]

    #dump header data to output file
    args.input.close()
    if 'header' not in variants:
        # usage("Header not defined in variant input file")
        parser.print_help()
        print("\nError: Header not defined in variant input file")
        exit(1)
    output = csv.writer(args.output, delimiter='\t', lineterminator='\n')

    output.writerow([
        *variants['header'],
        "GeneName",
        "HLAallele",
        "PeptideLength",
        "SubPeptidePosition",
        "MTScore",
        "WTScore",
        "MTEpitopeSeq",
        "WTEpitopeSeq",
        "FoldChange"
    ])

    #Read netmhc files from the fof, and parse into predictions
    prediction = {}

    intake = args.fof.readline().rstrip()


    while intake != "":

        #parse the filepath
        basename = os.path.basename(intake)
        data = basename.split(".")
        sample = netmhc_subber.sub("", data[0])
        allele = data[1]
        length = data[2]
        mode = 'filtered'

        if sample not in prediction:
            prediction[sample] = {}
        #open each file listed in the fof, and read gene data
        netmhc_file = open(intake, mode = 'r')
        netmhc_file.readline() #skip first line
        netmhc_reader = csv.reader(netmhc_file, delimiter='\t')
        for line_data in netmhc_reader:

            gene = {
                'gene_name' : line_data[0],
                'allele' : allele,
                'mode' : mode,
                'length' : length,
                'point_mutation' : line_data[1],
                'sub_peptide_position' : line_data[2],
                'mt_score' : line_data[3],
                'wt_score' : line_data[4],
                'mt_epitope_seq' : line_data[5],
                'wt_epitope_seq' : line_data[6],
                'fold_change' : line_data[7]
                }

            gene_name = gene['gene_name']
            gene_score = float(gene['mt_score'])

            #create the nested dictionary structure if necessary
            if (gene_name not in prediction[sample] or
                gene_score < prediction[sample][gene_name]['score']):
                prediction[sample][gene_name] = {
                    'genes' : [gene],
                    'score' : gene_score
                }
            elif gene_score == prediction[sample][gene_name]['score']:
                prediction[sample][gene_name]['genes'].append(gene)
        netmhc_file.close()

        intake = args.fof.readline().rstrip()
    args.fof.close()

    #Finish off by dumping the selections to the output file
    for entry in (
        #flatten the dictionary structure and iterate over each entry
        item for sample
        in sorted(prediction) for gene
        in sorted(prediction[sample]) for item
        in prediction[sample][gene]['genes']
        ):
        if (float(entry['mt_score']) < args.binding_threshold and
                float(entry['fold_change']) > args.minimum_fold_change):
            key = entry['gene_name'] + "\t" + entry['point_mutation']
            if key in variants:
                output.writerow([
                    *variants[key],
                    entry['gene_name'],
                    entry['allele'],
                    entry['length'],
                    entry['sub_peptide_position'],
                    entry['mt_score'],
                    entry['wt_score'],
                    entry['mt_epitope_seq'],
                    entry['wt_epitope_seq'],
                    entry['fold_change']
                ])
            else:
                print("Couldn't find variant for", gene,
                      entry['point_mutation'], "in variant file")
    args.output.close()


if __name__ == "__main__":
    main()

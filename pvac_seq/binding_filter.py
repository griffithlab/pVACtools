import argparse
import sys
import re
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=argparse.FileType('r'),
                        help="Input list of variants",
                        required=True)
    parser.add_argument('-f', '--fof', type=argparse.FileType('r'),
                        help="OF containing list of parsed epitope files " +
                        "for different allele-length combinations (same sample)",
                        required=True)
    parser.add_argument('-o', '--output', type=argparse.FileType('w'),
                        help="Output .xls file containing list of filtered " +
                        "epitopes based on binding affinity for each " +
                        "allele-length combination per gene",
                        required=True)
    parser.add_argument('-c', '--fold-change', type=int,
                        help="Minimum fold change between mutant binding " +
                        "score and wild-type score. The default is 0, which " +
                        "filters no results, but 1 is often a sensible " +
                        "default (requiring that binding is better to the MT " +
                        "than WT)",
                        default=0,
                        required=False,
                        dest="minimum_fold_change")
    parser.add_argument('-b', '--binding-threshold', type=int,
                        help="Report only epitopes where the mutant allele " +
                        "has ic50 binding scores below this value; default 500",
                        default=500,
                        required=False,
                        dest="binding_threshold")

    args = parser.parse_args()

    #precompile regex patterns used later
    chromosome_name = re.compile(r'^chromosome_name')
    netmhc_subber = re.compile(r"_netmhc")

    #open the variants file and parse into variants dictionary
    intake = args.input.readline().rstrip()
    variants = {}

    while intake != '':
        if chromosome_name.match(intake):
            variants['header'] = intake
            intake = args.input.readline().rstrip()
            continue

        data = intake.split('\t')

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
        variants[ gene + '\t' + amino_acid_change] = intake

        intake = args.input.readline().rstrip()

    #dump header data to output file
    args.input.close()
    if 'header' not in variants:
        # usage("Header not defined in variant input file")
        parser.print_help()
        print("\nError: Header not defined in variant input file")
        exit(1)
    args.output.write('\t'.join([
        variants['header'], "GeneName", "HLAallele",
        "PeptideLength", "SubPeptidePosition",
        "MTScore", "WTScore", "MTEpitopeSeq",
        "WTEpitopeSeq", "FoldChange"
        ]))
    args.output.write("\n")

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

        #open each file listed in the fof, and read gene data
        netmhc_reader = open(intake, mode = 'r')
        netmhc_reader.readline() #skip first line
        netmhc_intake = netmhc_reader.readline().rstrip()
        while netmhc_intake != "":

            line_data = netmhc_intake.split("\t")[:8]

            gene_name = {
                'gene_name' : line_data[0],
                'allele' : allele,
                'point_mutation' : line_data[1],
                'sub_peptide_position' : line_data[2],
                'mt_score' : line_data[3],
                'wt_score' : line_data[4],
                'mt_epitope_seq' : line_data[5],
                'wt_epitope_seq' : line_data[6],
                'fold_change' : line_data[7]
                }

            #create the nested dictionary structure if necessary
            if mode not in prediction:
                prediction[mode] = {}
            if sample not in prediction[mode]:
                prediction[mode][sample] = {}
            if length not in prediction[mode][sample]:
                prediction[mode][sample][length] = {}
            if 'genes' not in prediction[mode][sample][length]:
                prediction[mode][sample][length]['genes']=[]

            prediction[mode][sample][length]['genes'].append(gene_name)
            netmhc_intake = netmhc_reader.readline().rstrip()
        netmhc_reader.close()

        intake = args.fof.readline().rstrip()
    args.fof.close()

    best = {}
    #now select the best predictions
    for mode in sorted(prediction):
        for sample in sorted(prediction[mode]):
            for length in sorted(prediction[mode][sample]):
                for gene in prediction[mode][sample][length]['genes']:
                    if (sample in best and
                            gene['gene_name'] in best[sample] and
                            'SCORE' in best[sample][gene['gene_name']]):
                        if (float(gene['mt_score']) <
                                best[sample][gene['gene_name']]['SCORE']):
                            best[sample][gene['gene_name']]['SCORE'] = float(gene['mt_score'])
                            gene['sample'] = sample
                            gene['length'] = length
                            gene['mode'] = mode
                            best[sample][gene['gene_name']]['GENES'] = [gene]
                        elif (float(gene['mt_score']) ==
                                best[sample][gene['gene_name']]['SCORE']):
                            gene['sample'] = sample
                            gene['length'] = length
                            gene['mode'] = mode
                            best[sample][gene['gene_name']]['GENES'].append(gene)
                    else:
                        if sample not in best:
                            best[sample] = {}
                        if gene['gene_name'] not in best[sample]:
                            best[sample][gene['gene_name']] = {}
                        best[sample][gene['gene_name']]['SCORE'] = float(gene['mt_score'])
                        gene['sample'] = sample
                        gene['length'] = length
                        gene['mode'] = mode
                        best[sample][gene['gene_name']]['GENES'] = [gene]

    #Finish off by dumping the selections to the output file
    for sample in sorted(best):
        for gene in sorted(best[sample]):
            for entry in best[sample][gene]['GENES']:
                if (float(entry['mt_score']) < args.binding_threshold and
                        float(entry['fold_change']) > args.minimum_fold_change):
                    key = gene + "\t" + entry['point_mutation']
                    if key in variants:
                        args.output.write("\t".join([
                            variants[key],
                            gene,
                            entry['allele'],
                            entry['length'],
                            entry['sub_peptide_position'],
                            entry['mt_score'],
                            entry['wt_score'],
                            entry['mt_epitope_seq'],
                            entry['wt_epitope_seq'],
                            entry['fold_change']
                            ]))
                        args.output.write("\n")
                    else:
                        print("Couldn't find variant for", gene,
                              entry['point_mutation'], "in variant file")
    args.output.close()


if __name__ == "__main__":
    main()

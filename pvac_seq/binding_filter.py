import getopt
import sys
import re
import os

def usage(msg=""):
    print("Usage:", sys.argv[0])
    print("-i <input list of variants")
    print("-f <FOF containing list of parsed epitope files for different allele-length combinations (same sample)>")
    print("-o <Output .xls file containing list of filtered epitopes based on binding affinity for each allele-length combination per gene>")
    print("-c <Minimum fold change between mutant binding score and wild-type score. The default is 0, which filters no results, but 1 is often a sensible default (requiring that binding is better to the MT than WT)>")
    print("-b <report only epitopes where the mutant allele has ic50 binding scores below this value ; default 500>")
    if(msg != ""):
        print()
        print("Error:", msg)
    exit(1)

def main():
    #getopt argument parsing
    try:
        options, arguments = getopt.getopt(sys.argv[1:], "i:f:o:c:b:")
    except getopt.GetoptError as e:
        usage(e.msg)
    input_filename = ""
    fof_filename = ""
    output_filename = ""
    minimum_fold_change = 0
    binding_threshold = 500
    for opt, arg in options:
        if opt == '-i' :
            input_filename = arg
        elif opt == '-f' :
            fof_filename = arg
        elif opt == '-o':
            output_filename = arg
        elif opt == '-c' :
            try:
                minimum_fold_change = int(arg)
            except ValueError as e:
                usage("Minimum fold change must be an integer value")
        elif opt == '-b' :
            try:
                binding_threshold = int(arg)
            except ValueError as e:
                usage("Binding threshold must be an integer value")
        else:
            usage("unrecognized option " + opt)
    if len(arguments) > 0:
        usage("unrecognized trailing arguments: [" + ", ".join(arguments) + "]")
    if input_filename == "" or fof_filename == "" or output_filename == "":
        usage("Please provide all required parameters (-i -o -f)")

    #precompile regex patterns used later
    chromosome_name = re.compile(r'^chromosome_name')
    netmhc_subber = re.compile(r"_netmhc")

    #open the variants file and parse into variants dictionary
    input_file = open(input_filename, mode='r')
    intake = input_file.readline().rstrip()
    variants = {}

    while intake != '':
        if chromosome_name.match(intake):
            variants['header'] = intake
            intake = input_file.readline().rstrip()
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

        intake = input_file.readline().rstrip()

    #dump header data to output file
    input_file.close()
    if 'header' not in variants:
        usage("Header not defined in variant input file")
    output = open(output_filename, mode='w')
    output.write('\t'.join([
        variants['header'], "GeneName", "HLAallele",
        "PeptideLength", "SubPeptidePosition",
        "MTScore", "WTScore", "MTEpitopeSeq",
        "WTEpitopeSeq", "FoldChange"
        ]))
    output.write("\n")

    #open the fof file and parse into prediction dictionary
    fof_file = open(fof_filename, mode='r')

    prediction = {}

    intake = fof_file.readline().rstrip()


    while intake != "":

        #parse the filepath
        basename = os.path.basename(intake)
        data = basename.split(".")
        sample = netmhc_subber.sub("", data[0])
        allele = data[1]
        length = data[2]
        mode = 'filtered'
        # i = 0

        #open each file listed in the fof, and read gene data
        netmhc_reader = open(intake, mode = 'r')
        netmhc_reader.readline() #skip first line
        netmhc_intake = netmhc_reader.readline().rstrip()
        while netmhc_intake != "":
            # i+=1
            # if i==1:
            #     netmhc_intake = netmhc_reader.readline().rstrip()
            #     continue

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

        intake = fof_file.readline().rstrip()
    fof_file.close()

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
                if (float(entry['mt_score']) < binding_threshold and
                        float(entry['fold_change']) > minimum_fold_change):
                    key = gene + "\t" + entry['point_mutation']
                    if key in variants:
                        output.write("\t".join([
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
                        output.write("\n")
                    else:
                        print("Couldn't find variant for", gene,
                              entry['point_mutation'], "in variant file")
    output.close()


if __name__ == "__main__":
    main()

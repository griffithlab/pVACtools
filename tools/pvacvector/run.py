#  Input would be a protein fasta or both pvac-seq run output final.tsv and input annotated vcf
# python lib/vaccine_design.py test --generate-input-fasta -t tests/test_data/vaccine_design/input_parse_test_input.tsv -v tests/test_data/vaccine_design/input_parse_test_input.vcf ann H-2-Kb -o . -n 25 -l 8
# python lib/vaccine_design.py test -f tests/test_data/vaccine_design/Test.vaccine.results.input.fa ann H-2-Kb -o . -n 25 -l 8

import shutil
import sys
import argparse
import os
from pathlib import Path
root = str(Path(__file__).resolve().parents[1])
sys.path.append(root)

import pandas
import networkx as nx
import random
import lib
from lib.optimal_peptide import *
from lib.vector_visualization import *
from lib.prediction_class import *
from lib.run_argument_parser import *
from lib.pvacvector_input_fasta_generator import *
from lib.fasta_generator import *

def define_parser():
    return PvacvectorRunArgumentParser().parser

def main(args_input=sys.argv[1:]):

    parser = define_parser()
    args = parser.parse_args(args_input)

    if "." in args.sample_name:
        sys.exit("Run name cannot contain '.'")

    if args.iedb_retries > 100:
        sys.exit("The number of IEDB retries must be less than or equal to 100")

    if (os.path.splitext(args.input_file))[1] == '.fa':
        input_file = args.input_file
        generate_input_fasta = False
    elif (os.path.splitext(args.input_file))[1] == '.tsv':
        input_tsv = args.input_file
        input_vcf = args.input_vcf
        if input_vcf is None:
            sys.exit("Input VCF is required when using a pVACseq TSV as input file")
        generate_input_fasta = True
    else:
        sys.exit("Input file type not as expected. Needs to be a .fa or a .vcf file")
    input_n_mer = args.input_n_mer
    iedb_method = args.prediction_algorithms[0]
    ic50_cutoff = args.binding_threshold
    alleles = args.allele
    epl = args.epitope_length
    print("IC50 cutoff: " + str(ic50_cutoff))
    runname = args.sample_name
    outdir = args.output_dir

    base_output_dir = os.path.abspath(outdir)
    base_output_dir = os.path.join(base_output_dir, runname)
    tmp_dir = os.path.join(base_output_dir, runname + '_tmp')
    os.makedirs(tmp_dir, exist_ok=True)

    if os.environ.get('TEST_FLAG') or os.environ.get('TEST_FLAG') == '1':
        random.seed(0.5)
    if generate_input_fasta:
        generator = PvacvectorInputFastaGenerator(input_tsv, input_vcf, base_output_dir, input_n_mer)
        generator.execute()
        input_file = generator.output_file

    epitopes_file = os.path.join(tmp_dir, runname + "_epitopes.fa")
    params = {
        'input_file': input_file,
        'output_file': epitopes_file,
    }
    generator= VectorFastaGenerator(**params)
    generator.execute()
    epitopes = generator.epitopes
    seq_tuples = generator.seq_tuples
    seq_dict = generator.seq_dict
    seq_keys = generator.seq_keys

    outfile = os.path.join(tmp_dir, runname + '_iedb_out.csv')
    split_out = []

    for a in alleles:
        for l in epl:
            print ("Calling iedb for " + a + " of length " + str(l))
            prediction_class = globals()[iedb_method]
            prediction = prediction_class()
            translated_iedb_method = prediction.iedb_prediction_method
            lib.call_iedb.main([
                epitopes_file,
                outfile,
                translated_iedb_method,
                a,
                '-l', str(l),
                '-r', str(args.iedb_retries), 
                '-e', args.iedb_install_directory
            ])
            with open(outfile, 'rU') as sheet:
                split_out.append(pandas.read_csv(sheet, delimiter='\t'))

    print("IEDB calls complete. Merging data.")

    with open(outfile, 'rU') as sheet:
        split_out.append(pandas.read_csv(sheet, delimiter='\t'))
    epitope_binding = pandas.concat(split_out)
    problematic_neoepitopes = epitope_binding[epitope_binding.ic50 < ic50_cutoff]
    merged = pandas.DataFrame(pandas.merge(epitope_binding, problematic_neoepitopes, how='outer',
                                           indicator=True).query('_merge == "left_only"').drop(['_merge'], axis=1))
    merged = merged.sort_values('ic50', ascending=False)
    peptides = merged.set_index('peptide').T.to_dict('dict')

    keyErrorCount = 0
    successCount = 0
    iedb_results = dict()
    for seqID in epitopes:
        for l in epl:
            for i in range(0, len(epitopes[seqID]) - (l-1)):
                key = epitopes[seqID][i:i+l]
                try:
                    peptides[key]
                except KeyError:
                    keyErrorCount += 1
                    continue

                if seqID not in iedb_results:
                    iedb_results[seqID] = {}
                allele = peptides[key]['allele']
                if allele not in iedb_results[seqID]:
                    iedb_results[seqID][allele] = {}
                    if 'total_score' not in iedb_results[seqID][allele]:
                        iedb_results[seqID][allele]['total_score'] = list()
                        iedb_results[seqID][allele]['total_score'].append(peptides[key]['ic50'])
                    else:
                        iedb_results[seqID][allele]['total_score'].append(peptides[key]['ic50'])

                if 'min_score' in iedb_results[seqID][allele]:
                    iedb_results[seqID][allele]['min_score'] = min(iedb_results[seqID][allele]['min_score'], peptides[key]['ic50'])
                else:
                    iedb_results[seqID][allele]['min_score'] = peptides[key]['ic50']
                    successCount += 1

    print("Successful ic50 mappings: " + str(successCount) + " errors: " + str(keyErrorCount))

    Paths = nx.DiGraph()
    spacers = [None, "HH", "HHC", "HHH", "HHHD", "HHHC", "AAY", "HHHH", "HHAA", "HHL", "AAL"]
    for ep in seq_tuples:
        ID_1 = ep[0]
        ID_2 = ep[1]
        Paths.add_node(ID_1)
        Paths.add_node(ID_2)
        for space in spacers:
            if space is None:
                key = str(ID_1 + "|" + ID_2)
            else:
                key = str(ID_1 + "|" + space + "|" + ID_2)
            worst_case = sys.maxsize
            for allele in iedb_results[key]:
                if iedb_results[key][allele]['min_score'] < worst_case:
                    worst_case = iedb_results[key][allele]['min_score']
            if Paths.has_edge(ID_1, ID_2) and Paths[ID_1][ID_2]['weight'] < worst_case:
                Paths[ID_1][ID_2]['weight'] = worst_case
                if space is not None:
                    Paths[ID_1][ID_2]['spacer'] = space
                else:
                    Paths[ID_1][ID_2]['spacer'] = ''
            elif not Paths.has_edge(ID_1, ID_2):
                if space is not None:
                    Paths.add_edge(ID_1, ID_2, weight=worst_case, spacer=space)
                else:
                    Paths.add_edge(ID_1, ID_2, weight=worst_case, spacer='')

    print("Graph contains " + str(len(Paths)) + " nodes and " + str(Paths.size()) + " edges.")
    print("Finding path.")

    distance_matrix = {}
    for ID_1 in Paths:
        try:
            distance_matrix[ID_1]
        except KeyError:
            distance_matrix[ID_1] = {}
        for ID_2 in Paths[ID_1]:
            distance_matrix[ID_1][ID_2] = Paths[ID_1][ID_2]['weight']

    init_state = sorted(seq_dict)
    if not os.environ.get('TEST_FLAG') or os.environ.get('TEST_FLAG') == '0':
        random.shuffle(init_state)
    peptide = OptimalPeptide(init_state, distance_matrix)
    peptide.copy_strategy = "slice"
    peptide.save_state_on_exit = False
    state, e = peptide.anneal()
    while state[0] != seq_keys[0]:
        state = state[1:] + state[:1] 
    print("%i distance :" % e)

    for id in state:
        print("\t", id)

    results_file = os.path.join(base_output_dir, runname + '_results.fa')
    with open(results_file, 'w') as f:
        name = list()
        min_score = Paths[state[0]][state[1]]['weight']
        cumulative_weight = 0
        all_scores = list()

        for i in range(0, len(state)):
            name.append(state[i])
            try:
                min_score = min(min_score, Paths[state[i]][state[i + 1]]['weight'])
                cumulative_weight += Paths[state[i]][state[i + 1]]['weight']
                all_scores.append(str(Paths[state[i]][state[i + 1]]['weight']))
                spacer = Paths[state[i]][state[i + 1]]['spacer']
                if spacer is not '':
                    name.append(spacer)
            except IndexError:
                continue
        median_score = str(cumulative_weight/len(all_scores))
        peptide_id_list = ','.join(name)
        score_list = ','.join(all_scores)
        output = list()
        output.append(">")
        output.append(peptide_id_list)
        output.append("|Median_Junction_Score:")
        output.append(median_score)
        output.append("|Lowest_Junction_Score:")
        output.append(str(min_score))
        output.append("|All_Junction_Scores:")
        output.append(score_list)
        output.append("\n")
        for id in name:
            try:
                output.append(seq_dict[id])
            except KeyError:
                output.append(id)
            output.append("\n")
        f.write(''.join(output))

    if not args.keep_tmp_files:
        shutil.rmtree(tmp_dir)

    VectorVisualization(results_file, base_output_dir).draw()

if __name__ == "__main__":
    main()


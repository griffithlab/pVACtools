#  Input would be a protein fasta or both pvac-seq run output final.tsv and input annotated vcf
# python lib/vaccine_design.py test --generate-input-fasta -t tests/test_data/vaccine_design/input_parse_test_input.tsv -v tests/test_data/vaccine_design/input_parse_test_input.vcf ann H-2-Kb -o . -n 25 -l 8
# python lib/vaccine_design.py test -f tests/test_data/vaccine_design/Test.vaccine.results.input.fa ann H-2-Kb -o . -n 25 -l 8

import shutil
import sys
import argparse
import os
import pandas
import networkx as nx
import random
from Bio import SeqIO

from lib.optimal_peptide import *
from lib.vector_visualization import *
from lib.run_argument_parser import *
from lib.pvacvector_input_fasta_generator import *
from lib.pipeline import *
import lib.call_iedb

def define_parser():
    return PvacvectorRunArgumentParser().parser

def run_pipelines(input_file, base_output_dir, args):
    class_i_prediction_algorithms = []
    class_ii_prediction_algorithms = []
    for prediction_algorithm in sorted(args.prediction_algorithms):
        prediction_class = globals()[prediction_algorithm]
        prediction_class_object = prediction_class()
        if isinstance(prediction_class_object, MHCI):
            class_i_prediction_algorithms.append(prediction_algorithm)
        elif isinstance(prediction_class_object, MHCII):
            class_ii_prediction_algorithms.append(prediction_algorithm)

    class_i_alleles = []
    class_ii_alleles = []
    for allele in sorted(set(args.allele)):
        valid = 0
        if allele in MHCI.all_valid_allele_names():
            class_i_alleles.append(allele)
            valid = 1
        if allele in MHCII.all_valid_allele_names():
            class_ii_alleles.append(allele)
            valid = 1
        if not valid:
            print("Allele %s not valid. Skipping." % allele)

    shared_arguments = {
        'input_file'      : input_file,
        'input_file_type' : 'pvacvector_input_fasta',
        'sample_name'     : args.sample_name,
    }

    parsed_output_files = []
    if len(class_i_prediction_algorithms) > 0 and len(class_i_alleles) > 0:
        if args.epitope_length is None:
            sys.exit("Epitope length is required for class I binding predictions")

        if args.iedb_install_directory:
            iedb_mhc_i_executable = os.path.join(args.iedb_install_directory, 'mhc_i', 'src', 'predict_binding.py')
            if not os.path.exists(iedb_mhc_i_executable):
                sys.exit("IEDB MHC I executable path doesn't exist %s" % iedb_mhc_i_executable)
        else:
            iedb_mhc_i_executable = None

        print("Executing MHC Class I predictions")

        output_dir = os.path.join(base_output_dir, 'MHC_Class_I')
        os.makedirs(output_dir, exist_ok=True)

        class_i_arguments = shared_arguments.copy()
        class_i_arguments['alleles']                 = class_i_alleles
        class_i_arguments['peptide_sequence_length'] = ""
        class_i_arguments['iedb_executable']         = iedb_mhc_i_executable
        class_i_arguments['epitope_lengths']         = args.epitope_length
        class_i_arguments['prediction_algorithms']   = class_i_prediction_algorithms
        class_i_arguments['output_dir']              = output_dir
        pipeline_i = MHCIPipeline(**class_i_arguments)
        pipeline_i.generate_fasta([[1, 1]])
        parsed_output_files.extend(pipeline_i.call_iedb_and_parse_outputs([[1, 1]]))

    if len(class_ii_prediction_algorithms) > 0 and len(class_ii_alleles) > 0:
        if args.iedb_install_directory:
            iedb_mhc_ii_executable = os.path.join(args.iedb_install_directory, 'mhc_ii', 'mhc_II_binding.py')
            if not os.path.exists(iedb_mhc_ii_executable):
                sys.exit("IEDB MHC II executable path doesn't exist %s" % iedb_mhc_ii_executable)
        else:
            iedb_mhc_ii_executable = None

        print("Executing MHC Class II predictions")

        output_dir = os.path.join(base_output_dir, 'MHC_Class_II')
        os.makedirs(output_dir, exist_ok=True)

        class_ii_arguments = shared_arguments.copy()
        class_ii_arguments['alleles']               = class_ii_alleles
        class_ii_arguments['prediction_algorithms'] = class_ii_prediction_algorithms
        class_ii_arguments['iedb_executable']       = iedb_mhc_ii_executable
        class_ii_arguments['output_dir']            = output_dir
        class_ii_arguments['netmhc_stab']           = False
        pipeline_ii = MHCIIPipeline(**class_ii_arguments)
        pipeline_ii.generate_fasta([[1, 1]])
        parsed_output_files.extend(pipeline_ii.call_iedb_and_parse_outputs([[1, 1]]))

    return parsed_output_files

def find_min_scores(parsed_output_files, args):
    iedb_results = {}
    epitopes = []
    indexes = []
    for parsed_output_file in parsed_output_files:
        with open(parsed_output_file, 'r') as parsed:
            reader = csv.DictReader(parsed, delimiter="\t")
            for row in reader:
                index = row['Index']
                allele = row['HLA Allele']

                score = float(row['Best MT Score'])
                if args.allele_specific_binding_thresholds:
                    threshold = PredictionClass.cutoff_for_allele(entry[allele])
                    threshold = float(args.binding_threshold) if threshold is None else float(threshold)
                else:
                    threshold = float(args.binding_threshold)
                if score < threshold:
                    continue

                if index not in iedb_results:
                    iedb_results[index] = {}
                if allele not in iedb_results[index]:
                    iedb_results[index][allele] = {}
                if 'min_score' in iedb_results[index][allele]:
                    iedb_results[index][allele]['min_score'] = min(iedb_results[index][allele]['min_score'], score)
                else:
                    iedb_results[index][allele]['min_score'] = score
                epitopes.append(row['MT Epitope Seq'])
    return (iedb_results, epitopes)

def create_graph(iedb_results, seq_tuples):
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
            if key in iedb_results:
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
    return Paths

def create_distance_matrix(Paths):
    print("Finding path.")
    distance_matrix = {}
    for ID_1 in Paths:
        try:
            distance_matrix[ID_1]
        except KeyError:
            distance_matrix[ID_1] = {}
        for ID_2 in Paths[ID_1]:
            distance_matrix[ID_1][ID_2] = Paths[ID_1][ID_2]['weight']
    return distance_matrix

def find_optimal_path(Paths, distance_matrix, seq_dict, seq_keys, base_output_dir, args):
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

    results_file = os.path.join(base_output_dir, args.sample_name + '_results.fa')
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
    return results_file

def main(args_input=sys.argv[1:]):

    parser = define_parser()
    args = parser.parse_args(args_input)

    if "." in args.sample_name:
        sys.exit("Run name cannot contain '.'")

    if args.iedb_retries > 100:
        sys.exit("The number of IEDB retries must be less than or equal to 100")

    if args.iedb_install_directory:
        lib.call_iedb.setup_iedb_conda_env()

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
        sys.exit("Input file type not as expected. Needs to be a .fa or a .tsv file")

    base_output_dir = os.path.abspath(args.output_dir)
    os.makedirs(base_output_dir, exist_ok=True)

    if os.environ.get('TEST_FLAG') or os.environ.get('TEST_FLAG') == '1':
        random.seed(0.5)
    if generate_input_fasta:
        generator = PvacvectorInputFastaGenerator(input_tsv, input_vcf, base_output_dir, args.input_n_mer, args.sample_name)
        generator.execute()
        input_file = generator.output_file

    seq_dict = dict()
    for record in SeqIO.parse(input_file, "fasta"):
        seq_dict[record.id] = str(record.seq)
    seq_keys = sorted(seq_dict)
    seq_tuples = list(itertools.permutations(seq_keys, 2))

    parsed_output_files = run_pipelines(input_file, base_output_dir, args)
    (iedb_scores, epitopes) = find_min_scores(parsed_output_files, args)
    Paths = create_graph(iedb_scores, seq_tuples)
    distance_matrix = create_distance_matrix(Paths)
    results_file = find_optimal_path(Paths, distance_matrix, seq_dict, seq_keys, base_output_dir, args)
    if 'DISPLAY' in os.environ.keys():
        VectorVisualization(results_file, base_output_dir).draw()

    if not args.keep_tmp_files:
        shutil.rmtree(os.path.join(base_output_dir, 'MHC_Class_I'), ignore_errors=True)
        shutil.rmtree(os.path.join(base_output_dir, 'MHC_Class_II'), ignore_errors=True)

if __name__ == "__main__":
    main()


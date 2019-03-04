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
        'n_threads'       : args.n_threads,
        'spacers'         : args.spacers,
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
        pipeline_i = Pipeline(**class_i_arguments)
        pipeline_i.generate_fasta([[1, 1]])
        pipeline_i.call_iedb([[1, 1]])
        parsed_output_files.extend(pipeline_i.parse_outputs([[1, 1]]))

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
        class_ii_arguments['alleles']                 = class_ii_alleles
        class_ii_arguments['prediction_algorithms']   = class_ii_prediction_algorithms
        class_ii_arguments['peptide_sequence_length'] = ""
        class_ii_arguments['iedb_executable']         = iedb_mhc_ii_executable
        class_ii_arguments['epitope_lengths']         = [15]
        class_ii_arguments['output_dir']              = output_dir
        class_ii_arguments['netmhc_stab']             = False
        pipeline_ii = Pipeline(**class_ii_arguments)
        pipeline_ii.generate_fasta([[1, 1]])
        pipeline_ii.call_iedb([[1, 1]])
        parsed_output_files.extend(pipeline_ii.parse_outputs([[1, 1]]))

    return parsed_output_files

def find_min_scores(parsed_output_files, args):
    min_scores = {}
    indexes_with_good_binders = []
    #find indexes that contain a good binder so that they can be excluded from further processing
    #we don't want any peptide-spacer-peptide combination (aka index) that contains a good binder
    #Find min score of all the epitopes of each of the remaining peptide-spacer-peptide combinations 
    for parsed_output_file in parsed_output_files:
        with open(parsed_output_file, 'r') as parsed:
            reader = csv.DictReader(parsed, delimiter="\t")
            for row in reader:
                index = row['Index']
                if index in indexes_with_good_binders:
                    continue

                if args.top_score_metric == 'lowest':
                    score = float(row['Best MT Score'])
                elif args.top_score_metric == 'median':
                    score = float(row['Median MT Score'])
                if args.allele_specific_binding_thresholds:
                    allele = row['HLA Allele']
                    threshold = PredictionClass.cutoff_for_allele(allele)
                    threshold = float(args.binding_threshold) if threshold is None else float(threshold)
                else:
                    threshold = float(args.binding_threshold)
                if score < threshold:
                    indexes_with_good_binders.append(index)
                    continue

                if index in min_scores:
                    min_scores[index] = min(min_scores[index], score)
                else:
                    min_scores[index] = score

    for index in indexes_with_good_binders:
        if index in min_scores:
            del min_scores[index]

    return min_scores

def create_graph(iedb_results, seq_tuples, spacers):
    Paths = nx.DiGraph()
    for ep in seq_tuples:
        ID_1 = ep[0]
        ID_2 = ep[1]
        for spacer in spacers:
            if spacer == 'None':
                key = str(ID_1 + "|" + ID_2)
            else:
                key = str(ID_1 + "|" + spacer + "|" + ID_2)
            if key in iedb_results:
                worst_case = iedb_results[key]
            else:
                continue

            if not Paths.has_node(ID_1):
                Paths.add_node(ID_1)

            if not Paths.has_node(ID_2):
                Paths.add_node(ID_2)

            if Paths.has_edge(ID_1, ID_2) and Paths[ID_1][ID_2]['weight'] < worst_case:
                Paths[ID_1][ID_2]['weight'] = worst_case
                Paths[ID_1][ID_2]['spacer'] = spacer
            elif not Paths.has_edge(ID_1, ID_2):
                Paths.add_edge(ID_1, ID_2, weight=worst_case, spacer=spacer)

    print("Graph contains " + str(len(Paths)) + " nodes and " + str(Paths.size()) + " edges.")
    return Paths

def check_graph_valid(Paths):
    error_text = ('A vaccine design using the parameters specified could not be found.  Some options that you may want to consider:\n' +
                 '1) increasing the acceptable junction binding score to allow more possible connections (-b parameter)\n' +
                 '2) using the "median" binding score instead of the "best" binding score for each junction, (best may be too conservative, -m parameter)')

    n_nodes_without_outgoing_edges = 0
    for node in Paths.nodes():
        if len(Paths.out_edges(node)) == 0:
            n_nodes_without_outgoing_edges += 1
    if n_nodes_without_outgoing_edges > 1:
        raise Exception("Unable to create valid graph. No outgoing edges for more than one node.\n {}".format(error_text))

    n_nodes_without_incoming_edges = 0
    for node in Paths.nodes():
        if len(Paths.in_edges(node)) == 0:
            n_nodes_without_incoming_edges += 1
    if n_nodes_without_incoming_edges > 1:
        raise Exception("Unable to create valid graph. No incoming edges for more than one node.\n {}".format(error_text))

    n_nodes_without_any_edges = 0
    for node in Paths.nodes():
        if len(Paths.in_edges(node)) == 0 and len(Paths.out_edges(node)) == 0:
            n_nodes_without_any_edges += 1
    if n_nodes_without_any_edges > 0:
        raise Exception("Unable to create valid graph. No edges for at least one node.\n {}".format(error_text))

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
    init_state = sorted(Paths.nodes())
    if not os.environ.get('TEST_FLAG') or os.environ.get('TEST_FLAG') == '0':
        random.shuffle(init_state)
    peptide = OptimalPeptide(init_state, distance_matrix)
    peptide.copy_strategy = "slice"
    peptide.save_state_on_exit = False
    state, e = peptide.anneal()
    print("%i distance :" % e)

    for id in state:
        print("\t", id)

    results_file = os.path.join(base_output_dir, args.sample_name + '_results.fa')
    with open(results_file, 'w') as f:
        name = list()
        min_score = Paths[state[0]][state[1]]['weight']
        cumulative_weight = 0
        all_scores = list()

        for i in range(0, (len(state) - 1)):
            name.append(state[i])
            if Paths.has_edge(state[i], state[i + 1]):
                edge = Paths[state[i]][state[i + 1]]
                min_score = min(min_score, edge['weight'])
                cumulative_weight += edge['weight']
                all_scores.append(str(edge['weight']))
                spacer = edge['spacer']
                if spacer != 'None':
                    name.append(spacer)
            else:
                sys.exit("Unable to find path. All possible peptides for edge '{} - spacer - {}' contain at least one epitope that is a good binder.".format(state[i], state[i + 1]))
        name.append(state[-1])
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
    min_scores = find_min_scores(parsed_output_files, args)
    Paths = create_graph(min_scores, seq_tuples, args.spacers)
    check_graph_valid(Paths)
    distance_matrix = create_distance_matrix(Paths)
    results_file = find_optimal_path(Paths, distance_matrix, seq_dict, seq_keys, base_output_dir, args)
    if 'DISPLAY' in os.environ.keys():
        VectorVisualization(results_file, base_output_dir).draw()

    if not args.keep_tmp_files:
        shutil.rmtree(os.path.join(base_output_dir, 'MHC_Class_I'), ignore_errors=True)
        shutil.rmtree(os.path.join(base_output_dir, 'MHC_Class_II'), ignore_errors=True)

if __name__ == "__main__":
    main()


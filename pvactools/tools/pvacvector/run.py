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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import itertools
import json
import platform

from pvactools.lib.optimal_peptide import OptimalPeptide
from pvactools.lib.vector_visualization import VectorVisualization
from pvactools.lib.run_argument_parser import PvacvectorRunArgumentParser
from pvactools.lib.pvacvector_input_fasta_generator import PvacvectorInputFastaGenerator
from pvactools.lib.pipeline import *
from pvactools.lib.run_utils import *
from pvactools.lib.prediction_class_utils import *

def define_parser():
    return PvacvectorRunArgumentParser().parser

def run_pipelines(input_file, base_output_dir, args, junctions_to_test, spacer, tries, class_i_prediction_algorithms, class_ii_prediction_algorithms, class_i_alleles, class_ii_alleles):
    shared_arguments = {
        'input_file'      : input_file,
        'input_file_type' : 'pvacvector_input_fasta',
        'sample_name'     : args.sample_name,
        'n_threads'       : args.n_threads,
        'spacer'          : spacer,
        'clip_length'     : tries,
        'downstream_sequence_length': 200,
        'iedb_retries'    : args.iedb_retries,
        'additional_report_columns' : None,
        'junctions_to_test': junctions_to_test,
    }

    parsed_output_files = []
    if len(class_i_prediction_algorithms) > 0 and len(class_i_alleles) > 0:
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
        class_i_arguments['iedb_executable']         = iedb_mhc_i_executable
        class_i_arguments['epitope_lengths']         = args.class_i_epitope_length
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
        class_ii_arguments['iedb_executable']         = iedb_mhc_ii_executable
        class_ii_arguments['epitope_lengths']         = args.class_ii_epitope_length
        class_ii_arguments['output_dir']              = output_dir
        class_ii_arguments['netmhc_stab']             = False
        pipeline_ii = Pipeline(**class_ii_arguments)
        pipeline_ii.generate_fasta([[1, 1]])
        pipeline_ii.call_iedb([[1, 1]])
        parsed_output_files.extend(pipeline_ii.parse_outputs([[1, 1]]))

    return parsed_output_files

def write_junctions_file(graph, current_output_dir):
    junctions_file = os.path.join(current_output_dir, 'junctions.tsv')
    with open(junctions_file, 'w') as fh:
        fieldnames = ['left_peptide', 'left_partner_clip', 'spacer', 'right_partner_clip', 'right_peptide', 'junction_score', 'percentile']
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for (left_peptide, right_peptide, edge_data) in graph.edges.data():
            row = {
                'left_peptide': left_peptide,
                'left_partner_clip': edge_data['left_partner_trim'],
                'spacer': edge_data['spacer'],
                'right_partner_clip': edge_data['right_partner_trim'],
                'right_peptide': right_peptide,
                'junction_score': edge_data['weight'],
                'percentile': edge_data['percentile'],
            }
            writer.writerow(row)

def find_min_scores(parsed_output_files, current_output_dir, args, min_scores, min_percentiles):
    #min_scores_rows = {}
    junctions_with_good_binders = set()
    #find junctions that contain a good binder so that they can be excluded from further processing
    #we don't want any peptide|left_clip|spacer|right_clip|peptide combination (aka junctions) that contains a good binder
    #Find min score of all the epitopes of each of the remaining peptide-spacer-peptide combinations 
    processed_junctions = set()
    junction_min_scores = {}
    junction_min_percentiles = {}
    for parsed_output_file in parsed_output_files:
        with open(parsed_output_file, 'r') as parsed:
            reader = csv.DictReader(parsed, delimiter="\t")
            for row in reader:
                index = row['Mutation']
                processed_junctions.add(index)

                if args.top_score_metric == 'lowest':
                    score = float(row['Best IC50 Score'])
                    percentile = float(row['Best Percentile'])
                elif args.top_score_metric == 'median':
                    score = float(row['Median IC50 Score'])
                    percentile = float(row['Median Percentile'])
                if args.allele_specific_binding_thresholds:
                    allele = row['HLA Allele']
                    threshold = PredictionClass.cutoff_for_allele(allele)
                    threshold = float(args.binding_threshold) if threshold is None else float(threshold)
                else:
                    threshold = float(args.binding_threshold)
                if score < threshold:
                    junctions_with_good_binders.add(index)
                if args.percentile_threshold is not None and percentile < args.percentile_threshold:
                    junctions_with_good_binders.add(index)

                if index not in junction_min_scores:
                    #"initialize" with the first score encountered
                    junction_min_scores[index] = score
                elif score < junction_min_scores[index]:
                    #if the current score is lower than the saved one, update the saved one
                    junction_min_scores[index] = score
                else:
                    continue
                if index not in junction_min_percentiles:
                    #"initialize" with the first percentile encountered
                    junction_min_percentiles[index] = percentile
                elif percentile < junction_min_percentiles[index]:
                    #if the current percentile is lower than the saved one, update the saved one
                    junction_min_percentiles[index] = percentile
                else:
                    continue
    good_junctions = processed_junctions - junctions_with_good_binders
    for good_junction in good_junctions:
        min_scores[good_junction] = junction_min_scores[good_junction]
        min_percentiles[good_junction] = junction_min_percentiles[good_junction]

    return min_scores, min_percentiles

def initialize_graph(seq_keys):
    graph = nx.DiGraph()
    for key in seq_keys:
        graph.add_node(key)
    return graph

def add_valid_junctions_to_graph(graph, min_scores, min_percentiles):
    for key, worst_case in min_scores.items():
        percentile = min_percentiles[key]
        if (key.count("|") == 3):
            (id_1, left_partner_trim, right_partner_trim, id_2) = key.split("|")
            spacer = "None"
        elif (key.count("|") == 4):
            (id_1, left_partner_trim, spacer, right_partner_trim, id_2) = key.split("|")
        if graph.has_edge(id_1, id_2) and graph[id_1][id_2]['weight'] < worst_case:
            graph[id_1][id_2]['weight'] = worst_case
            graph[id_1][id_2]['percentile'] = percentile
            graph[id_1][id_2]['spacer'] = spacer
            graph[id_1][id_2]['left_partner_trim'] = int(left_partner_trim)
            graph[id_1][id_2]['right_partner_trim'] = int(right_partner_trim)
        elif not graph.has_edge(id_1, id_2):
            graph.add_edge(id_1, id_2, weight=worst_case, percentile=percentile, spacer=spacer, left_partner_trim=int(left_partner_trim), right_partner_trim=int(right_partner_trim))
    return graph

def check_graph_valid(Paths, seq_dict):
    graph_valid = True
    errors = []
    if len(Paths.nodes()) < len(seq_dict.keys()):
        graph_valid = False
        errors.append("No valid junctions found for peptides: {}".format(set(seq_dict.keys()) - set(Paths.nodes())))

    nodes_without_outgoing_edges = []
    nodes_without_incoming_edges = []
    nodes_without_any_edges = []
    for node in Paths.nodes():
        if len(Paths.out_edges(node)) == 0:
            nodes_without_outgoing_edges.append(node)
        if len(Paths.in_edges(node)) == 0:
            nodes_without_incoming_edges.append(node)
        if len(Paths.in_edges(node)) == 0 and len(Paths.out_edges(node)) == 0:
            nodes_without_any_edges.append(node)

    if len(nodes_without_outgoing_edges) > 1:
        graph_valid = False
        errors.append("More than one peptide without valid outgoing junction: {}".format(nodes_without_outgoing_edges))
    if len(nodes_without_incoming_edges) > 1:
        graph_valid = False
        errors.append("More than one peptide without valid incoming junction: {}".format(nodes_without_incoming_edges))
    if len(nodes_without_any_edges) > 0:
        graph_valid = False
        errors.append("No valid junctions found for peptides: {}".format(nodes_without_any_edges))

    if graph_valid:
        return(graph_valid, None)
    else:
        return(graph_valid, "\n".join(errors))

def identify_problematic_junctions(graph, seq_tuples):
    tuples_to_process = []
    for (seq1, seq2) in seq_tuples:
        if not graph.has_edge(seq1, seq2):
            tuples_to_process.append([seq1, seq2])
    return tuples_to_process

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

def find_optimal_path(graph, distance_matrix, seq_dict, base_output_dir, args):
    init_state = sorted(graph.nodes())
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
        names = list()
        cumulative_weight = 0
        all_scores = list()

        problematic_junctions = []
        for i in range(0, (len(state) - 1)):
            names.append(state[i])
            if graph.has_edge(state[i], state[i + 1]):
                edge = graph[state[i]][state[i + 1]]
                try:
                    min_score = min(min_score, edge['weight'])
                except:
                    min_score = edge['weight']
                cumulative_weight += edge['weight']
                all_scores.append(str(edge['weight']))
                spacer = edge['spacer']
                if spacer != 'None':
                    names.append(spacer)
            else:
                problematic_junctions.append("{} - {}".format(state[i], state[i+1]))
        if len(problematic_junctions) > 0:
            return (None, "No valid junction between peptides: {}".format(", ".join(problematic_junctions)))
        names.append(state[-1])
        median_score = str(cumulative_weight/len(all_scores))
        peptide_id_list = ','.join(names)
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
        for (i, name) in enumerate(names):
            if name in args.spacers:
                output.append(name)
            else:
                full_sequence = seq_dict[name]
                if i > 0:
                    left_name = names[i-1]
                    if left_name in args.spacers:
                        left_name = names[i-2]
                    left_clip = graph.edges[(left_name, name)]['right_partner_trim']
                else:
                    left_clip = 0
                if i < (len(names) - 1):
                    right_name = names[i+1]
                    if right_name in args.spacers:
                        right_name = names[i+2]
                    right_clip = graph.edges[(name, right_name)]['left_partner_trim']
                else:
                    right_clip = 0
                clipped_sequence = full_sequence[left_clip:len(full_sequence) - right_clip]
                output.append(clipped_sequence)

            output.append("\n")
        f.write(''.join(output))
    return (results_file, None)

def get_codon_for_amino_acid(amino_acid):
    amino_acid_to_codon = {
       'A': 'GCC',
       'C': 'TGC',
       'D': 'GAC',
       'E': 'GAG',
       'F': 'TTC',
       'G': 'GGC',
       'H': 'CAC',
       'I': 'ATC',
       'K': 'AAG',
       'L': 'CTG',
       'M': 'ATG',
       'N': 'AAC',
       'P': 'CCC',
       'Q': 'CAG',
       'R': 'AGA',
       'S': 'AGC',
       'T': 'ACC',
       'V': 'GTG',
       'W': 'TGG',
       'Y': 'TAC',
    }
    return amino_acid_to_codon[amino_acid]

def create_dna_backtranslation(results_file, dna_results_file):
    record = SeqIO.read(results_file, 'fasta')
    seq_num = record.id
    peptide = str(record.seq)
    dna_sequence = ""
    for amino_acid in peptide:
        dna_sequence += get_codon_for_amino_acid(amino_acid)
    output_record = SeqRecord(Seq(dna_sequence), id=str(seq_num), description=str(seq_num))
    SeqIO.write([output_record], dna_results_file, 'fasta')

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
        sys.exit("Input file type not as expected. Needs to be a .fa or a .tsv file")

    if args.n_threads > 1 and platform.system() == "Darwin":
        raise Exception("Multithreading is not supported on MacOS")

    (class_i_prediction_algorithms, class_ii_prediction_algorithms) = split_algorithms(args.prediction_algorithms)
    if len(class_i_prediction_algorithms) == 0:
        print("No MHC class I prediction algorithms chosen. Skipping MHC class I predictions.")
    elif len(class_ii_prediction_algorithms) == 0:
        print("No MHC class II prediction algorithms chosen. Skipping MHC class II predictions.")

    alleles = combine_class_ii_alleles(args.allele)
    (class_i_alleles, class_ii_alleles, species) = split_alleles(alleles)
    if len(class_i_alleles) == 0:
        print("No MHC class I alleles chosen. Skipping MHC class I predictions.")
    elif len(class_ii_alleles) == 0:
        print("No MHC class II alleles chosen. Skipping MHC class II predictions.")

    if len(class_i_prediction_algorithms) == 0 and len(class_i_alleles) == 0 and len(class_ii_prediction_algorithms) == 0 and len(class_ii_alleles) == 0:
        return

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
    graph = initialize_graph(seq_keys)

    seq_tuples = list(itertools.permutations(seq_keys, 2))
    junctions_to_process = seq_tuples

    max_tries = args.max_clip_length + 1
    tries = 0
    results_file = None
    min_scores = {}
    min_percentiles = {}
    while results_file is None and tries < max_tries:
        print("Processing clip length {}".format(tries))
        for spacer in args.spacers:
            print("Processing spacer {}".format(spacer))
            current_output_dir = os.path.join(base_output_dir, str(tries), spacer)
            parsed_output_files = run_pipelines(
                input_file,
                current_output_dir,
                args,
                junctions_to_process,
                spacer,
                tries,
                class_i_prediction_algorithms,
                class_ii_prediction_algorithms,
                class_i_alleles,
                class_ii_alleles
            )
            min_scores, min_percentiles = find_min_scores(parsed_output_files, current_output_dir, args, min_scores, min_percentiles)
            add_valid_junctions_to_graph(graph, min_scores, min_percentiles)
            write_junctions_file(graph, current_output_dir)
            (valid, error) = check_graph_valid(graph, seq_dict)
            if not valid:
                junctions_to_process = identify_problematic_junctions(graph, seq_tuples)
                print("No valid path found. {}".format(error))
                continue
            distance_matrix = create_distance_matrix(graph)
            (results_file, error) = find_optimal_path(graph, distance_matrix, seq_dict, base_output_dir, args)
            if results_file is None:
                print("No valid path found. {}".format(error))
                junctions_to_process = identify_problematic_junctions(graph, seq_tuples)
            else:
                junctions_file = os.path.join(current_output_dir, 'junctions.tsv')
                import shutil
                shutil.copy(junctions_file, base_output_dir)
                break
        tries += 1

    if results_file is None:
        print(
            'Unable to find path. ' +
            'A vaccine design using the parameters specified could not be found.  Some options that you may want to consider:\n' +
            '1) decreasing the acceptable junction binding score to allow more possible connections (-b parameter)\n' +
            '2) using the "median" binding score instead of the "best" binding score for each junction, (best may be too conservative, -m parameter)\n' +
            '3) if running with a percentile threshold set, either remove this parameter or reduce the acceptable threshold to allow more possible connections (--percentile-threshold parameter)'
        )
    else:
        if 'DISPLAY' in os.environ.keys():
            VectorVisualization(results_file, base_output_dir, args.spacers).draw()

        dna_results_file = os.path.join(base_output_dir, args.sample_name + '_results.dna.fa')
        create_dna_backtranslation(results_file, dna_results_file)

        if not args.keep_tmp_files:
            for subdirectory in range(tries):
                for spacer in processed_spacers:
                    shutil.rmtree(os.path.join(base_output_dir, str(subdirectory), spacer, 'MHC_Class_I'), ignore_errors=True)
                    shutil.rmtree(os.path.join(base_output_dir, str(subdirectory), spacer, 'MHC_Class_II'), ignore_errors=True)

    change_permissions_recursive(base_output_dir, 0o755, 0o644)

if __name__ == "__main__":
    main()


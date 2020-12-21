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
from Bio.Alphabet import IUPAC
from lib.optimal_peptide import *
from lib.vector_visualization import *
from lib.run_argument_parser import *
from lib.pvacvector_input_fasta_generator import *
from lib.pipeline import *
from lib.run_utils import *
import lib.call_iedb

def define_parser():
    return PvacvectorRunArgumentParser().parser

def run_pipelines(input_file, base_output_dir, args, spacer):
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
        'spacers'         : [spacer],
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

def write_min_scores(min_scores_rows, directory, args):
    #This will write the junction scores for all tested spacers
    min_scores_file = os.path.join(directory, 'junction_scores.tsv')
    rows = []
    with open(min_scores_file, 'w') as fh:
        fieldnames = ['left_peptide', 'spacer', 'right_peptide', 'junction_score', 'epitope', 'allele', 'method']
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in min_scores_rows:
            index_parts = row['Mutation'].split('|')
            left_peptide = index_parts[0]
            if len(index_parts) == 2:
                spacer = 'None'
                right_peptide = index_parts[1]
            else:
                spacer = index_parts[1]
                right_peptide = index_parts[2]
            new_row = {
                'left_peptide': left_peptide,
                'spacer': spacer,
                'right_peptide': right_peptide,
                'epitope': row['Epitope Seq'],
                'allele': row['HLA Allele'],
            }
            if args.top_score_metric == 'lowest':
                new_row['junction_score'] = float(row['Best Score'])
                new_row['method'] = row['Best Score Method']
            elif args.top_score_metric == 'median':
                new_row['junction_score'] = float(row['Median Score'])
                new_row['method'] = 'median'
            rows.append(new_row)
        sorted_rows = sorted(rows, key=lambda k: k['junction_score'])
        writer.writerows(sorted_rows)

    #This will filter `rows` to only the ones with the best spacer for each junction
    best_spacers_min_scores = {}
    best_spacers_min_scores_rows = {}
    for row in rows:
        index = (row['left_peptide'], row['right_peptide'])
        score = row['junction_score']
        if index in best_spacers_min_scores and score >= best_spacers_min_scores[index]:
            continue
        else:
            best_spacers_min_scores[index] = score
            best_spacers_min_scores_rows[index] = row
    sorted_best_spacers_min_scores_rows = sorted(best_spacers_min_scores_rows.values(), key=lambda k: k['junction_score'])
    best_spacers_min_scores_file = os.path.join(directory, 'junction_scores.best_spacers.tsv')
    with open(best_spacers_min_scores_file, 'w') as fh:
        fieldnames = ['left_peptide', 'spacer', 'right_peptide', 'junction_score', 'epitope', 'allele', 'method']
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(sorted_best_spacers_min_scores_rows)

def find_min_scores(parsed_output_files, current_output_dir, args):
    min_scores = {}
    min_scores_rows = {}
    indexes_with_good_binders = []
    #find indexes that contain a good binder so that they can be excluded from further processing
    #we don't want any peptide-spacer-peptide combination (aka index) that contains a good binder
    #Find min score of all the epitopes of each of the remaining peptide-spacer-peptide combinations 
    for parsed_output_file in parsed_output_files:
        with open(parsed_output_file, 'r') as parsed:
            reader = csv.DictReader(parsed, delimiter="\t")
            for row in reader:
                index = row['Mutation']

                if args.top_score_metric == 'lowest':
                    score = float(row['Best Score'])
                elif args.top_score_metric == 'median':
                    score = float(row['Median Score'])
                if args.allele_specific_binding_thresholds:
                    allele = row['HLA Allele']
                    threshold = PredictionClass.cutoff_for_allele(allele)
                    threshold = float(args.binding_threshold) if threshold is None else float(threshold)
                else:
                    threshold = float(args.binding_threshold)
                if score < threshold:
                    indexes_with_good_binders.append(index)

                if index in min_scores and score >= min_scores[index]:
                    continue
                else:
                    min_scores[index] = score
                    min_scores_rows[index] = row

    write_min_scores(min_scores_rows.values(), current_output_dir, args)

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

    print("Graph contains " + str(len(Paths)) + " nodes (peptides) and " + str(Paths.size()) + " edges (junctions).")
    return Paths

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
        problematic_ends = []
        problematic_starts = []
        problematic_junctions = []
        cumulative_weight = 0
        all_scores = list()

        for i in range(0, (len(state) - 1)):
            name.append(state[i])
            if Paths.has_edge(state[i], state[i + 1]):
                edge = Paths[state[i]][state[i + 1]]
                try:
                    min_score = min(min_score, edge['weight'])
                except:
                    min_score = edge['weight']
                cumulative_weight += edge['weight']
                all_scores.append(str(edge['weight']))
                spacer = edge['spacer']
                if spacer != 'None':
                    name.append(spacer)
            else:
                problematic_ends.append(state[i])
                problematic_starts.append(state[i+1])
                problematic_junctions.append("{} - {}".format(state[i], state[i+1]))
        if len(problematic_junctions) > 0:
            return (None, "No valid junction between peptides: {}".format(", ".join(problematic_junctions)), problematic_starts, problematic_ends)
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
    return (results_file, None, None, None)

def shorten_problematic_peptides(input_file, problematic_start, problematic_end, output_dir):
    print("Shortening problematic peptides")
    records = []
    for record in SeqIO.parse(input_file, "fasta"):
        if record.id in problematic_start and record.id in problematic_end:
            record_new = SeqRecord(Seq(str(record.seq)[1:-1], IUPAC.protein), id=record.id, description=record.description)
        elif record.id in problematic_start:
            record_new = SeqRecord(Seq(str(record.seq)[1:], IUPAC.protein), id=record.id, description=record.description)
        elif record.id in problematic_end:
            record_new = SeqRecord(Seq(str(record.seq)[:-1], IUPAC.protein), id=record.id, description=record.description)
        else:
            record_new = record
        records.append(record_new)
    os.makedirs(output_dir, exist_ok=True)
    new_input_file = os.path.join(output_dir, "vector_input.fa")
    SeqIO.write(records, new_input_file, "fasta")
    return new_input_file

def identify_problematic_peptides(Paths, seq_dict):
    problematic_start = set(seq_dict.keys()) - set(Paths.nodes())
    problematic_end = set(seq_dict.keys()) - set(Paths.nodes())
    for node in Paths.nodes():
        if len(Paths.out_edges(node)) == 0 and len(Paths.in_edges(node)) == 0:
            problematic_start.add(node)
            problematic_end.add(node)
        elif len(Paths.out_edges(node)) == 0:
            problematic_end.add(node)
        elif len(Paths.in_edges(node)) == 0:
            problematic_start.add(node)
    return (problematic_start, problematic_end)

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

    results_file = None
    max_tries = args.max_clip_length + 1
    tries = 0
    while results_file is None and tries < max_tries:
        if tries > 0:
            input_file = shorten_problematic_peptides(input_file, problematic_start, problematic_end, os.path.join(base_output_dir, str(tries)))
        seq_dict = dict()
        for record in SeqIO.parse(input_file, "fasta"):
            seq_dict[record.id] = str(record.seq)
        seq_keys = sorted(seq_dict)
        seq_tuples = list(itertools.permutations(seq_keys, 2))

        all_parsed_output_files = []
        processed_spacers = []
        results_file = None
        for spacer in args.spacers:
            print("Processing spacer {}".format(spacer))
            processed_spacers.append(spacer)
            current_output_dir = os.path.join(base_output_dir, str(tries), spacer)
            parsed_output_files = run_pipelines(input_file, current_output_dir, args, spacer)
            all_parsed_output_files.extend(parsed_output_files)
            min_scores = find_min_scores(all_parsed_output_files, current_output_dir, args)
            Paths = create_graph(min_scores, seq_tuples, processed_spacers)
            (valid, error) = check_graph_valid(Paths, seq_dict)
            if not valid:
                (problematic_start, problematic_end) = identify_problematic_peptides(Paths, seq_dict)
                print("No valid path found. {}".format(error))
                continue
            distance_matrix = create_distance_matrix(Paths)
            (results_file, error, problematic_start, problematic_end) = find_optimal_path(Paths, distance_matrix, seq_dict, seq_keys, base_output_dir, args)
            if results_file is not None:
                break
            else:
                print("No valid path found. {}".format(error))
        tries += 1

    if results_file is None:
        raise Exception(
            'Unable to find path. ' +
            'A vaccine design using the parameters specified could not be found.  Some options that you may want to consider:\n' +
            '1) increasing the acceptable junction binding score to allow more possible connections (-b parameter)\n' +
            '2) using the "median" binding score instead of the "best" binding score for each junction, (best may be too conservative, -m parameter)'
        )

    if 'DISPLAY' in os.environ.keys():
        VectorVisualization(results_file, base_output_dir, args.spacers).draw()

    if not args.keep_tmp_files:
        shutil.rmtree(os.path.join(base_output_dir, 'MHC_Class_I'), ignore_errors=True)
        shutil.rmtree(os.path.join(base_output_dir, 'MHC_Class_II'), ignore_errors=True)

    change_permissions_recursive(base_output_dir, 0o755, 0o644)

if __name__ == "__main__":
    main()


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
import itertools
from Bio import SeqIO
import lib

from lib.optimal_peptide import *
import random

from lib.prediction_class import *

import turtle

def define_parser():

    parser = argparse.ArgumentParser('pvacvector run')
    parser.add_argument(
        "run_name",
        help="The name of the run being processed." +
             " This will be used as a prefix for output files."
    )
    parser.add_argument('-g', "--generate-input-fasta", action="store_true")
    parser.add_argument('-f', "--input-fa", type=argparse.FileType('r'),
                       help="Path to input fasta file. " + "Required if not generating input fasta")
    parser.add_argument('-t', "--input_tsv",
                        help="Path to input tsv file with the epitopes selected for vector design. " + "Required if generating input fasta. ")
    parser.add_argument('-v', "--input_vcf",
                        help="Path to original pVAC-Seq input vcf file" + "Requiired if generating input fasta. ")
    parser.add_argument('method',
                        choices=PredictionClass.iedb_prediction_methods(),
                        help="The iedb analysis method to use")
    parser.add_argument('allele',
                        help="Allele for which to make prediction")
    parser.add_argument('-o', "--outdir", help="Output directory")
    parser.add_argument('-k', "--keep-tmp",
                        help="Option to store tmp files. ", action="store_true")
    parser.add_argument('-n', "--input-n-mer", default='25',
                        help="Length of peptide sequence to be generated for use in input to main vector design algorithm. " + "Default: 25")
    parser.add_argument(
        "-l", "--epitope-length",
        type=lambda s: [int(epl) for epl in s.split(',')],
        help="Length of subpeptides (neoepitopes) to predict. " +
             "Multiple epitope lengths can be specified " +
             "using a comma-separated list. " +
             "Typical epitope lengths vary between 8-11. " +
             "Required for Class I prediction algorithms",
    )
    parser.add_argument("-c", "--cutoff", type=int,
                        default=500,
                        help="Optional ic50 cutoff value." +
                             " Junctional neoepitopes with IC50 values below this " +
                             " value will be excluded. Default: 500")
    parser.add_argument(
        "-r", "--iedb-retries", type=int,
        default=5,
        help="Number of retries when making requests to the IEDB RESTful web interface. " +
             " Must be less than or equal to 100." +
             "Default: 5"
    )
    parser.add_argument(
        "-e", "--iedb-executable-path",
        help="The executable path of the local IEDB install"
    )
    parser.add_argument(
        "-s", "--seed-rng", action="store_true",
        help="Seed random number generator with default value 0.5 for unit test." +
        " Default: False")
    return parser

def tsvToFasta(n_mer, input_tsv, input_vcf, output_dir):

    def parse_choosen_epitopes(input_tsv):
        with open(input_tsv, 'r') as input_f:
            next(input_f)
            mut_IDs, mutations, mut_types, mt_epitope_seqs, wt_epitope_seqs, transcript_IDs = [], [], [], [], [], []
            for line in input_f:
                fields = line.split("\t")
                mut_type, mutation, pos, gene_name = fields[7], fields[8], fields[9], fields[10]
                mt_epitope_seq, wt_epitope_seq = fields[15], fields[16]
                mutations.append(mutation)
                mut_types.append(mut_type)
                mutation = mutation.split("/")
                #if position presented as a range, use higher end of range
                old_AA, new_AA = mutation[0], mutation[1]
                if "-" in pos:
                    pos = pos.split("-")
                    pos = pos[1]
                    mut_ID = ("MT." + gene_name + "." +  pos + "fs")
                elif mut_type == "FS":
                    mut_ID = "MT." + gene_name + "." + old_AA + pos + "fs"
                elif mut_type == "missense": 
                    mut_ID = "MT." + gene_name + "."  + old_AA + pos + new_AA
                mut_IDs.append(mut_ID)
                mt_epitope_seqs.append(mt_epitope_seq)
                wt_epitope_seqs.append(wt_epitope_seq)
                transcript_IDs.append(fields[5])
        input_f.close()
        return mut_IDs, mutations, mut_types, mt_epitope_seqs, wt_epitope_seqs, transcript_IDs

#get necessary data from initial pvacseq input vcf
    def parse_original_vcf(input_vcf):
        with open(input_vcf, 'r') as input_f:
            transcripts_dict = {}
            for line in input_f:
                attributes = []
                if line[0] != "#":
                    fields = line.split("\t")
                    info = fields[7]
                    info = info.split("|")
                    transcript_ID, downstr_seq, len_change, full_seq = info[6], info[23], info[24], info[25]
                    attributes.append(full_seq)
                    attributes.append(downstr_seq)
                    attributes.append(len_change)
                    transcripts_dict[transcript_ID] = attributes
        input_f.close()
        return(transcripts_dict)

    def edit_full_seq(i, mut_types, mutations, wt_epitope_seqs, mt_epitope_seqs, sub_seq, full_seq, transcripts_dict):
        if mut_types[i] == "FS":
            downstr_seq, len_change = transcripts_dict[transcript_IDs[i]][1], int(transcripts_dict[transcript_IDs[i]][2])
            parts = mutations[i].split("/")
            initial = parts[0]
            final = parts[1]

            #handle -/X mutations by appending downstr_seq to next position, instead of overwriting last position
            if initial == "-":
                new_end_of_full_seq = len(full_seq) + len_change - len(downstr_seq)
            #overwrites last position of seq with first position of
            #predicted downstr seq
            else:
                new_end_of_full_seq = len(full_seq) + len_change - len(downstr_seq) - 1
            full_seq = full_seq[:new_end_of_full_seq]
        #handles ex: L/LX mutations by adding sequence that is preserved
        #before the downstr predicted sequence
            if len(final) > 1:
                final = final.replace("X", "")
                full_seq = full_seq + final
            full_seq = full_seq + downstr_seq
        elif mut_types[i] == "missense":
            full_seq = full_seq.replace(wt_epitope_seqs[i], mt_epitope_seqs[i])
        else:
            sys.exit("Mutation not yet handled by this parser")
        return(full_seq)

    #get flanking peptides for the epitope chosen
    def get_sub_seq(full_seq, mt_seq, n_mer):
        beginning = full_seq.find(mt_seq)
        if beginning == -1:
            sys.exit("Error: could not find mutant epitope sequence in mutant full sequence")
        length = len(mt_seq)
        end = beginning + length
        #if eptitope sequence is too close to the beginning or end to get the
        #right amount of flanking peptides, get appropriate length from solely
        #ahead or behind
        len_needed = n_mer - length
        if len_needed % 2 != 0:
            front = int(beginning - len_needed / 2)
            back = int(end + len_needed / 2)
        else:
            front = int(beginning - len_needed / 2)
            back = int(end + len_needed / 2)
        if front < 0:
            sub_seq = full_seq[beginning:(beginning + n_mer)]
        elif back > len(full_seq):
            sub_seq = full_seq[(end - n_mer):end]
        else:
            sub_seq = full_seq[front:back]
        return(sub_seq)

    def write_output_fasta(output_f, n_mer, mut_IDs, mutations, mut_types, mt_epitope_seqs, wt_epitope_seqs, transcript_IDs, transcripts_dict):
        with open(output_f, 'w') as out_f:
            sub_seq = ""
            full_seq = ""
            n_mer = int(n_mer)
            for i in range(len(transcript_IDs)):
                full_seq = (transcripts_dict[transcript_IDs[i]])[0] 
            
                full_seq = edit_full_seq(i, mut_types, mutations, wt_epitope_seqs, mt_epitope_seqs, sub_seq, full_seq, transcripts_dict)

                sub_seq = get_sub_seq(full_seq, mt_epitope_seqs[i], n_mer)
                out_f.write(">" + mut_IDs[i] + "\n")
                out_f.write(sub_seq + "\n")
                print("ID: " + mut_IDs[i] + ", sequence: " + sub_seq)
        out_f.close()
        print("FASTA file written")
        return()

    output_f = os.path.join(output_dir, "vector_input.fa")

    mut_IDs, mutations, mut_types, mt_epitope_seqs, wt_epitope_seqs, transcript_IDs = parse_choosen_epitopes(input_tsv)

    transcripts_dict = parse_original_vcf(input_vcf)

    write_output_fasta(output_f, n_mer, mut_IDs, mutations, mut_types, mt_epitope_seqs, wt_epitope_seqs, transcript_IDs, transcripts_dict)
    return(output_f)

def parse_input(input_file):
    pep_seqs = []
    with open(input_file, 'r') as input_f:
        header = input_f.readline().strip()
        for line in input_f:
            pep_seqs.append(line.strip())
    #remove >, get peptide names and junction scores from FASTA input
    edited_header = header[1:]
    fields = edited_header.split("|")
    pep_ids_joined = fields[0].split(",")
    pep_ids = []
    for pep_id in pep_ids_joined:
        if "." in pep_id:
            mt, gene, var = pep_id.split(".")
            pep_id = "-".join((gene,var))
        pep_ids.append(pep_id)
    junct_scores = fields[3].split(":")
    junct_scores = junct_scores[1].split(",")
    return(pep_seqs, pep_ids, junct_scores)

#determine number of peptides (not including junction additions)
def get_peptide_num(pep_seqs, min_pep_length, max_pep_length):    
    num_peptides = 0
    for pep in pep_seqs:
        length = len(pep)
        if length > min_pep_length and length < max_pep_length:
            num_peptides += 1
    return(num_peptides)

#determine what proportion of circle each peptide should take up
def get_conversion_factor(pep_seqs):
    total_len = 0
    for pep in pep_seqs:
        total_len += len(pep)   
    #30 degrees reserved for white space
    conversion_factor = 330 / total_len
    return(conversion_factor)

def draw_header(t, header_pos):
    t.pu()
    t.setpos(header_pos)
    t.write("Vector Design", align="center", font=("Arial", 18, "bold"))
    t.pd()
    return()

def draw_wht_space(t, circle_radius, wht_space_angle):
    t.pencolor("white")
    t.circle(circle_radius,wht_space_angle)
    return()

#select color from scheme
def get_color(count):
    #option 1: 3 blue/green color scheme
    #scheme = [(161,218,180),(65,182,196),(44,127,184),(37,52,148)]
    #option 2: 6 color scheme
    scheme = [(31,120,180),(51,160,44),(227,26,28),(255,127,0),(106,61,154),(177,89,40)]
    count =  count % len(scheme)
    return scheme[count]\

def write_junct_score(t, junct_score, size):
    t.back(size)
    t.write(junct_score + 'nM', align="center")
    t.forward(size)

#draw perpindicular line to arc to mark junction
def draw_junction_w_label(junct_score, t, pen_thin, angle):
    reset = t.heading()
    t.rt(90)
    t.pencolor("black")
    t.pensize(pen_thin)
    t.forward(10)
    t.back(20)
    t.pu()
    if (angle >= 0 and angle < 70):
        write_junct_score(t, junct_score, 15)
    elif (angle >= 65 and angle < 115):
        write_junct_score(t, junct_score, 10)
    elif (angle >= 115 and angle < 165):
        write_junct_score(t, junct_score, 20)
    elif (angle >= 165 and angle < 195):
        write_junct_score(t, junct_score, 25)
    elif (angle >= 195 and angle < 245):
        write_junct_score(t, junct_score, 35)
    elif (angle >= 245 and angle < 295):
        write_junct_score(t, junct_score, 15)
    else:
        write_junct_score(t, junct_score, 20)
    t.pd()
    t.forward(10)
    t.setheading(reset)
    return()

#draw second line of junctions with amino acid additions
def draw_junction(t, pen_thin):
    reset = t.heading()
    t.rt(90)
    t.pencolor("black")
    t.pensize(pen_thin)
    t.forward(10)
    t.back(20)
    t.forward(10)
    t.setheading(reset)
    return()

#draw arc for peptide
def draw_arc_peptide(peptide, length, count, angle, t, circle_radius, conversion_factor, pep_id_space):
    t.pencolor(get_color(count))
    t.circle(circle_radius, (conversion_factor * length) / 2)
    t.pu()
    reset = t.heading()
    t.left(90)
    t.forward(pep_id_space)
    if (angle > 80 and angle < 100) or (angle > 260 and angle < 280):
        t.write(peptide, align="center", font=("Arial", 10, "bold"))
    elif (angle > 0 and angle < 90) or (angle > 270 and angle < 360):
        t.write(peptide, align="right", font=("Arial", 10, "bold"))
    else:
        t.write(peptide, align="left", font=("Arial", 10, "bold"))
    t.back(pep_id_space)
    t.setheading(reset)
    t.pd()
    t.circle(circle_radius, (conversion_factor * length) / 2)

#draw arc for amino acid inserts to junctions
def draw_arc_junct(peptide, length, t, conversion_factor, junct_seq_space, circle_radius):
    #t.pencolor("black")
    t.circle(circle_radius, (conversion_factor * length) / 2)
    t.pu()
    reset = t.heading()
    t.left(90)
    t.back(junct_seq_space)
    t.write(peptide, align="center")
    t.forward(junct_seq_space)
    t.setheading(reset)
    t.pd()
    t.circle(circle_radius, (conversion_factor * length) / 2)

def draw_peptide(t, pep, pep_ids, peptides_parsed, pen_thick, angle_parsed, conversion_factor, min_pep_length, max_pep_length, junctions_parsed, junct_scores, circle_radius, pep_id_space, junct_seq_space, pen_thin):
    junction_parsed = 0
    pep_length = len(pep)
    peptide = pep_ids[peptides_parsed]
    t.pensize(pen_thick)
    angle_parsed += conversion_factor * pep_length
    #if length within reasonable range, draw and label arc for peptide
    if pep_length > min_pep_length and pep_length < max_pep_length:
        draw_arc_peptide(peptide, pep_length, junctions_parsed, angle_parsed, t, circle_radius, conversion_factor, pep_id_space)
        if junctions_parsed < len(junct_scores):
            draw_junction_w_label(junct_scores[junctions_parsed], t, pen_thin, angle_parsed)
            junction_parsed += 1
    #if length is less than minimum peptide length, assume amino acid addition to junction
    elif pep_length < min_pep_length:
        draw_arc_junct(peptide, pep_length, t, conversion_factor, junct_seq_space, circle_radius)
        draw_junction(t, pen_thin)
    else:
        print("Error: Peptide sequence over 100 amino acids inputted")
        sys.exit()
    return(junction_parsed, angle_parsed)

#print turtle screen to a postscript file, convert to pdf
def output_screen(t, out_f):
    ps_file = "/".join((out_f, "vector.ps"))
    out_file = "/".join((out_f, "vector.jpg"))
    ts = t.getscreen()
    ts.getcanvas().postscript(file=ps_file)
    os.system('convert -density 300 -quality 100 ' + ps_file + " " + out_file)
    os.system('rm ' + ps_file)
    return()

def output_vaccine_png(input_file, out_f):
    min_pep_length = 8
    max_pep_length = 100

    pep_seqs, pep_ids, junct_scores = parse_input(input_file)

    #Error if not a peptide sequence for every peptide ID
    if len(pep_ids) != len(pep_seqs):
        print("Error: Not an equal number of peptide sequences and peptide IDs")
        sys.exit()

    conversion_factor = get_conversion_factor(pep_seqs)
    num_peptides = get_peptide_num(pep_seqs, min_pep_length, max_pep_length)

    #draw vaccine
    #800,600
    turtle.setup(800,600)
    t = turtle.Turtle()
    myWin = turtle.Screen()
    turtle.colormode(255)
    turtle.mode("logo")
    t.speed(0)
    t.hideturtle()

    #negative radius draws circle clockwise
    circle_radius = -200
    circle_pos = (-200,0)
    header_pos = (0,0)
    wht_space_angle = 15
    pen_thick = 5
    pen_thin = 2
    pep_id_space = 45 + num_peptides
    junct_seq_space = 25

    draw_header(t, header_pos)
    t.pu()
    t.setpos(circle_pos)
    t.pd()
    t.pensize(pen_thick)
    angle_parsed = 0
    #add white space in circle before genes
    draw_wht_space(t, circle_radius, wht_space_angle)
    angle_parsed += wht_space_angle
    draw_junction(t, pen_thin)

    #draw main part of circle
    junctions_parsed = 0
    peptides_parsed = 0
    for pep in pep_seqs:
        junction_parsed, angle_parsed = draw_peptide(t, pep, pep_ids, peptides_parsed, pen_thick, angle_parsed, conversion_factor, min_pep_length, max_pep_length, junctions_parsed, junct_scores, circle_radius, pep_id_space, junct_seq_space, pen_thin)
        junctions_parsed += junction_parsed 
        peptides_parsed += 1

    #add white space in circle after genes    
    draw_junction(t, pen_thin)
    draw_wht_space(t, circle_radius, wht_space_angle)
    output_screen(t, out_f)
    
    #keeps turtle screen open until closed by user
    #turtle.mainloop()

def main(args_input=sys.argv[1:]):

    parser = define_parser()
    args = parser.parse_args(args_input)

    if "." in args.run_name:
        sys.exit("Run name cannot contain '.'")

    if args.iedb_retries > 100:
        sys.exit("The number of IEDB retries must be less than or equal to 100")

    input_tsv = args.input_tsv
    input_vcf = args.input_vcf
    input_file = args.input_fa
    input_n_mer = args.input_n_mer
    iedb_method = args.method
    ic50_cutoff = args.cutoff
    alleles = args.allele.split(',')
    epl = args.epitope_length
    print("IC50 cutoff: " + str(ic50_cutoff))
    runname = args.run_name
    outdir = args.outdir
    generate_input_fasta = args.generate_input_fasta

    base_output_dir = os.path.abspath(outdir)
    base_output_dir = os.path.join(base_output_dir, runname)
    tmp_dir = os.path.join(base_output_dir, runname + '_tmp')
    os.makedirs(tmp_dir, exist_ok=True)

    if args.seed_rng:
        random.seed(0.5)
    if generate_input_fasta:
        input_file = tsvToFasta(input_n_mer, input_tsv, input_vcf, base_output_dir)

    peptides = SeqIO.parse(input_file, "fasta")
   
    seq_dict = dict()
    for record in peptides:
        seq_dict[record.id] = str(record.seq)

    seq_keys = sorted(seq_dict)
    seq_tuples = list(itertools.combinations_with_replacement(seq_keys, 2))
    combinations = list()

    for key in seq_tuples:
        if key[0] != key[1]:
            combinations.append((key[0], key[1]))
            combinations.append((key[1], key[0]))
    
    seq_tuples = combinations
    epitopes = dict()
    rev_lookup = dict()

    for comb in seq_tuples:
        seq1 = comb[0]
        seq2 = comb[1]
        for length in range(8, 11):
            seq_ID = seq1 + "|" + seq2
            trunc_seq1 = seq_dict[seq1][(len(seq_dict[seq1]) - length):len(seq_dict[seq1])]
            trunc_seq2 = seq_dict[seq2][0:(length - 1)]
            epitopes[seq_ID] = trunc_seq1 + trunc_seq2
            rev_lookup[(trunc_seq1 + trunc_seq2)] = seq_ID

            spacers = ["HH", "HHC", "HHH", "HHHD", "HHHC", "AAY", "HHHH", "HHAA", "HHL", "AAL"]
            for this_spacer in spacers:
                seq_ID = seq1 + "|" + this_spacer + "|" + seq2
                epitopes[seq_ID] = (trunc_seq1 + this_spacer + trunc_seq2)
                rev_lookup[(trunc_seq1 + this_spacer + trunc_seq2)] = seq_ID

    epitopes_file = os.path.join(tmp_dir, runname + "_epitopes.fa")
    with open(epitopes_file, "w") as tmp:
        for each in epitopes:
            tmp.write(">" + each + "\n" + epitopes[each] + "\n")

    outfile = os.path.join(tmp_dir, runname + '_iedb_out.csv')
    split_out = []

    for a in alleles:
        for l in epl:
            print ("Calling iedb for " + a + " of length " + str(l))
            lib.call_iedb.main([
                epitopes_file,
                outfile,
                iedb_method,
                a,
                '-l', str(l),
                '-r', str(args.iedb_retries), 
                '-e', args.iedb_executable_path
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
    for ep in combinations:
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

    init_state = seq_keys
    if not args.seed_rng:
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

    if not args.keep_tmp:
        shutil.rmtree(tmp_dir) 

    #keep this change when merging
    output_vaccine_png(results_file, base_output_dir)

if __name__ == "__main__":
    main()


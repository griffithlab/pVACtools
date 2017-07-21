#  Input would be a protein fasta
# python vaccine_design.py test peptides.fa ann H-2-Kb -o /Users/user/Desktop -l 8
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

import time
from simanneal import Annealer
import random
import math

from lib.prediction_class import *

import turtle

#gives graphics layer when running on a headless system
os.environ["DISPLAY"] = ':0.0'

def define_parser():

    parser = argparse.ArgumentParser('pvacseq vaccine_design')
    parser.add_argument(
        "run_name",
        help="The name of the run being processed." +
             " This will be used as a prefix for output files."
    )
    parser.add_argument('input_file', type=argparse.FileType('r'),
                        help="Path to input FASTA file")
    parser.add_argument('method',
                        choices=PredictionClass.iedb_prediction_methods(),
                        help="The iedb analysis method to use")
    parser.add_argument('allele',
                        help="Allele for which to make prediction")
    parser.add_argument('-o', "--outdir", help="Output directory")
    parser.add_argument('-k', "--keep-tmp",
                        help="Option to store tmp files. ", action="store_true")
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


#https://github.com/perrygeo/simanneal/blob/master/simanneal/anneal.py
class OptimalPeptide(Annealer):

    def __init__(self, state, distance_matrix):
        self.distance_matrix = distance_matrix
        super(OptimalPeptide, self).__init__(state)  # important!

    def move(self):
        """Swaps two peptides in the path."""
        a = random.randint(0, len(self.state) - 1)
        b = random.randint(0, len(self.state) - 1)
        self.state[a], self.state[b] = self.state[b], self.state[a]

    def energy(self):
        """Calculates the length of the route."""
        e = 0
        for i in range(len(self.state)):
            e += self.distance_matrix[self.state[i - 1]][self.state[i]]
        return e

    def anneal(self):
        """Minimizes the energy of a system by simulated annealing.
        Parameters
        state : an initial arrangement of the system
        Returns
        (state, energy): the best state and energy found.
        """
        step = 0
        self.start = time.time()

        # Precompute factor for exponential cooling from Tmax to Tmin
        if self.Tmin <= 0.0:
            raise Exception('Exponential cooling requires a minimum "\
                "temperature greater than zero.')
        Tfactor = math.log(self.Tmax / self.Tmin)

        # Note initial state
        T = self.Tmax
        E = self.energy()
        prevState = self.copy_state(self.state)
        prevEnergy = E
        self.best_state = self.copy_state(self.state)
        self.best_energy = E
        trials, accepts, improves = 0, 0, 0
        if self.updates > 0:
            updateWavelength = self.steps / self.updates
            self.update(step, T, E, None, None)

        # Attempt moves to new states
        while step < self.steps and not self.user_exit:
            step += 1
            T = self.Tmax * math.exp(Tfactor * step / self.steps)
            self.move()
            E = self.energy()
            dE = E - prevEnergy
            trials += 1
            if dE < 0.0 and math.exp(dE / T) < random.random():
                # Restore previous state
                self.state = self.copy_state(prevState)
                E = prevEnergy
            else:
                # Accept new state and compare to best state
                accepts += 1
                if dE > 0.0:
                    improves += 1
                prevState = self.copy_state(self.state)
                prevEnergy = E
                if E > self.best_energy:
                    self.best_state = self.copy_state(self.state)
                    self.best_energy = E
            if self.updates > 1:
                if (step // updateWavelength) > ((step - 1) // updateWavelength):
                    self.update(
                        step, T, E, accepts / trials, improves / trials)
                    trials, accepts, improves = 0, 0, 0

        # line break after progress output
        print('')

        self.state = self.copy_state(self.best_state)
        if self.save_state_on_exit:
            self.save_state()
        # Return best state and energy
        return self.best_state, self.best_energy

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
    t.write("Vaccine Design", align="center", font=("Arial", 18, "bold"))
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

#print turtle screen to a postscript file, convert to pdf
def output_screen(t, out_f):
    ps_file = "/".join((out_f, "vaccine.ps"))
    out_file = "/".join((out_f, "vaccine.jpg"))
    ts = t.getscreen()
    ts.getcanvas().postscript(file=ps_file)
    os.system('convert -density 300 -quality 100 ' + ps_file + " " + out_file)
    os.system('rm ' + ps_file)
    return()

def output_vaccine_png(input_file, out_f):
    min_pep_length = 8
    max_pep_length = 100

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

    #Error if not a peptide sequence for every peptide ID
    if len(pep_ids) != len(pep_seqs):
        print("Error: Not an equal number of peptide sequences and peptide IDs")
        sys.exit()

    conversion_factor = get_conversion_factor(pep_seqs)

    #determine number of peptides (not including junction additions)
    num_peptides = 0
    for pep in pep_seqs:
        length = len(pep)
        if length > 8 and length < 25:
            num_peptides += 1
            

    #draw vaccine
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
        pep_length = len(pep)
        peptide = pep_ids[peptides_parsed]
        t.pensize(pen_thick)
        angle_parsed += conversion_factor * pep_length
    #if length within reasonable range, draw and label arc for peptide
        if pep_length > min_pep_length and pep_length < max_pep_length:
            draw_arc_peptide(peptide, pep_length, junctions_parsed, angle_parsed, t, circle_radius, conversion_factor, pep_id_space)
            if junctions_parsed < len(junct_scores):
                draw_junction_w_label(junct_scores[junctions_parsed], t, pen_thin, angle_parsed)
                junctions_parsed += 1
    #if length is less than minimum peptide length, assume amino acid addition to junction
        elif pep_length < min_pep_length:
                draw_arc_junct(peptide, pep_length, t, conversion_factor, junct_seq_space, circle_radius)
                draw_junction(t, pen_thin)
        else:
                print("Error: Peptide sequence over 100 amino acids inputted")
                sys.exit()
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

    input_file = args.input_file
    iedb_method = args.method
    ic50_cutoff = args.cutoff
    alleles = args.allele.split(',')
    epl = args.epitope_length
    print("IC50 cutoff: " + str(ic50_cutoff))
    runname = args.run_name
    outdir = args.outdir

    base_output_dir = os.path.abspath(outdir)
    tmp_dir = os.path.join(base_output_dir, runname, runname + '_tmp')
    os.makedirs(tmp_dir, exist_ok=True)

    if args.seed_rng:
        random.seed(0.5)

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

    results_file = os.path.join(base_output_dir, runname, runname + '_results.fa')
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

    out_f = "/".join((base_output_dir, runname)) 

    output_vaccine_png(results_file, out_f)

if __name__ == "__main__":
    main()


import turtle
import os
import sys

class VectorVisualization:
    def __init__(self, input_fasta, output_directory):
        self.input_fasta = input_fasta
        self.output_directory = output_directory

        self.min_pep_length = 8
        self.max_pep_length = 100

        self.pep_seqs, self.pep_ids, self.junct_scores = self.parse_input()

        #Error if not a peptide sequence for every peptide ID
        if len(self.pep_ids) != len(self.pep_seqs):
            sys.exit("Error: Not an equal number of peptide sequences and peptide IDs")

        self.conversion_factor = self.get_conversion_factor()
        self.num_peptides = self.get_peptide_num()

        self.turtle = turtle.Turtle()

        self.circle_radius = -200
        self.circle_pos = (-200,0)
        self.header_pos = (0,0)
        self.wht_space_angle = 15
        self.pen_thick = 5
        self.pen_thin = 2
        self.pep_id_space = 45 + self.num_peptides
        self.junct_seq_space = 25

    def parse_input(self):
        pep_seqs = []
        with open(self.input_fasta, 'r') as input_f:
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
                desc = pep_id.split(".")
                if len(desc) == 3: #if in expected "MT.GENE.AACHANGE" format, simplify
                    mt, gene, var = pep_id.split(".")
                    pep_id = "-".join((gene,var))
            pep_ids.append(pep_id)
        junct_scores = fields[3].split(":")
        junct_scores = junct_scores[1].split(",")
        return(pep_seqs, pep_ids, junct_scores)

    #determine what proportion of circle each peptide should take up
    def get_conversion_factor(self):
        total_len = 0
        for pep in self.pep_seqs:
            total_len += len(pep)
        #30 degrees reserved for white space
        conversion_factor = 330 / total_len
        return conversion_factor

    #determine number of peptides (not including junction additions)
    def get_peptide_num(self):
        num_peptides = 0
        for pep in self.pep_seqs:
            length = len(pep)
            if length > self.min_pep_length and length < self.max_pep_length:
                num_peptides += 1
        return(num_peptides)

    def draw(self):
        self.setup_turtle()
        self.draw_header()
        self.reset_turtle()

        angle_parsed = 0
        #add white space in circle before genes
        self.draw_wht_space()
        angle_parsed += self.wht_space_angle
        self.draw_junction()

        #draw main part of circle
        junctions_parsed = 0
        peptides_parsed = 0
        for pep in self.pep_seqs:
            junction_parsed, angle_parsed = self.draw_peptide(pep, peptides_parsed, junctions_parsed, angle_parsed)
            junctions_parsed += junction_parsed
            peptides_parsed += 1

        #add white space in circle after genes
        self.draw_junction()
        self.draw_wht_space()
        self.output_screen()

        #keeps turtle screen open until closed by user
        #turtle.mainloop()

    def setup_turtle(self):
        turtle.setup(800,600)
        #myWin = turtle.Screen()
        turtle.colormode(255)
        turtle.mode("logo")
        self.turtle.speed(0)
        self.turtle.hideturtle()

    def draw_header(self):
        self.turtle.pu()
        self.turtle.setpos(self.header_pos)
        self.turtle.write("Vector Design", align="center", font=("Arial", 18, "bold"))
        self.turtle.pd()

    def reset_turtle(self):
        self.turtle.pu()
        self.turtle.setpos(self.circle_pos)
        self.turtle.pd()
        self.turtle.pensize(self.pen_thick)

    def draw_wht_space(self):
        self.turtle.pencolor("white")
        self.turtle.circle(self.circle_radius, self.wht_space_angle)

    #draw second line of junctions with amino acid additions
    def draw_junction(self):
        reset = self.turtle.heading()
        self.turtle.rt(90)
        self.turtle.pencolor("black")
        self.turtle.pensize(self.pen_thin)
        self.turtle.forward(10)
        self.turtle.back(20)
        self.turtle.forward(10)
        self.turtle.setheading(reset)

    def draw_peptide(self, pep, peptides_parsed, junctions_parsed, angle_parsed):
        junction_parsed = 0
        pep_length = len(pep)
        peptide = self.pep_ids[peptides_parsed]
        self.turtle.pensize(self.pen_thick)
        angle_parsed += self.conversion_factor * pep_length
        #if length within reasonable range, draw and label arc for peptide
        if pep_length >= self.min_pep_length and pep_length <= self.max_pep_length:
            self.draw_arc_peptide(peptide, pep_length, junctions_parsed, angle_parsed)
            if junctions_parsed < len(self.junct_scores):
                self.draw_junction_w_label(self.junct_scores[junctions_parsed], angle_parsed)
                junction_parsed += 1
        #if length is less than minimum peptide length, assume amino acid addition to junction
        elif pep_length < self.min_pep_length:
            self.draw_arc_junct(peptide, pep_length)
            self.draw_junction()
        else:
            sys.exit("Error: Peptide sequence over 100 amino acids inputted")
        return(junction_parsed, angle_parsed)

    #draw arc for peptide
    def draw_arc_peptide(self, peptide, length, count, angle):
        self.turtle.pencolor(self.get_color(count))
        self.turtle.circle(self.circle_radius, (self.conversion_factor * length) / 2)
        self.turtle.pu()
        reset = self.turtle.heading()
        self.turtle.left(90)
        self.turtle.forward(self.pep_id_space)
        if (angle > 80 and angle < 100) or (angle > 260 and angle < 280):
            self.turtle.write(peptide, align="center", font=("Arial", 10, "bold"))
        elif (angle > 0 and angle < 90) or (angle > 270 and angle < 360):
            self.turtle.write(peptide, align="right", font=("Arial", 10, "bold"))
        else:
            self.turtle.write(peptide, align="left", font=("Arial", 10, "bold"))
        self.turtle.back(self.pep_id_space)
        self.turtle.setheading(reset)
        self.turtle.pd()
        self.turtle.circle(self.circle_radius, (self.conversion_factor * length) / 2)

    #draw perpindicular line to arc to mark junction
    def draw_junction_w_label(self, junct_score, angle):
        reset = self.turtle.heading()
        self.turtle.rt(90)
        self.turtle.pencolor("black")
        self.turtle.pensize(self.pen_thin)
        self.turtle.forward(10)
        self.turtle.back(20)
        self.turtle.pu()
        if (angle >= 0 and angle < 70):
            self.write_junct_score(junct_score, 15)
        elif (angle >= 65 and angle < 115):
            self.write_junct_score(junct_score, 10)
        elif (angle >= 115 and angle < 165):
            self.write_junct_score(junct_score, 20)
        elif (angle >= 165 and angle < 195):
            self.write_junct_score(junct_score, 25)
        elif (angle >= 195 and angle < 245):
            self.write_junct_score(junct_score, 35)
        elif (angle >= 245 and angle < 295):
            self.write_junct_score(junct_score, 15)
        else:
            self.write_junct_score(junct_score, 20)
        self.turtle.pd()
        self.turtle.forward(10)
        self.turtle.setheading(reset)

    def write_junct_score(self, junct_score, size):
        self.turtle.back(size)
        self.turtle.write("{}nM".format(round(float(junct_score), 5)) , align="center")
        self.turtle.forward(size)

    #draw arc for amino acid inserts to junctions
    def draw_arc_junct(self, peptide, length):
        #t.pencolor("black")
        self.turtle.circle(self.circle_radius, (self.conversion_factor * length) / 2)
        self.turtle.pu()
        reset = self.turtle.heading()
        self.turtle.left(90)
        self.turtle.back(self.junct_seq_space)
        self.turtle.write(peptide, align="center")
        self.turtle.forward(self.junct_seq_space)
        self.turtle.setheading(reset)
        self.turtle.pd()
        self.turtle.circle(self.circle_radius, (self.conversion_factor * length) / 2)

    #print turtle screen to a postscript file, convert to pdf
    def output_screen(self):
        ps_file = os.path.join(self.output_directory, "vector.ps")
        out_file = os.path.join(self.output_directory, "vector.jpg")
        ts = self.turtle.getscreen()
        ts.getcanvas().postscript(file=ps_file)
        os.system('convert -density 300 -quality 100 ' + ps_file + " " + out_file)
        os.remove(ps_file)

    #select color from scheme
    def get_color(self, count):
        #option 1: 3 blue/green color scheme
        #scheme = [(161,218,180),(65,182,196),(44,127,184),(37,52,148)]
        #option 2: 6 color scheme
        scheme = [(31,120,180),(51,160,44),(227,26,28),(255,127,0),(106,61,154),(177,89,40)]
        count =  count % len(scheme)
        return scheme[count]\

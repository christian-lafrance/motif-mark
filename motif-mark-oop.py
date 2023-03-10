#!/usr/bin/env python

import re
import cairo
import argparse
from bioinfo import onelineFasta

class Gene:

    """
    The Gene class is used to create gene objects from FASTA sequences. One Gene object
    is created for each FASTA sequence. Two methods, getExonPos and findMotifs, are used
    to retrieve relevant sequence info. These methods return None, but set relevant 
    attributes for each Gene instance. 
    """
    
    def __init__(self, header:str, seq:str):
        self.header = header # the fasta header
        self.seq = seq # the gene sequence
        self.seq_len = len(seq) # the total length of the sequence
        self.name = header.split(" ")[0][1:] # the gene name from the header


    def getExonPos(self) -> None:
        '''
        Parses a gene sequence, searches for the exons which are in all caps,
        and returns the 1-based start position of each exon. Returns None but sets
        the exonPos attribute of the Gene object. This is a dictionary
        with exon number as the key and a tuple (start pos, end pos) as the value.
        Positions are in 1-based numbering. 
        '''
        exons = {}
        i = 0 # counter for current position in sequence
        e = 1 # keep track of number of exons encountered
        found_exon = False
        while i < len(self.seq):
            if self.seq[i].isupper() == False:
                if found_exon:
                    e += 1
                    found_exon = False
                next
            else:
                found_exon = True
                if f"Exon_{e}" in exons:
                    exons[f"Exon_{e}"][1] += 1
                else: 
                    # add 1 to i for 1-based position numbering
                    exons[f"Exon_{e}"] = [i+1, i+1]
            i += 1
        self.exonPos = exons
        return


    def findMotifs(self, motifList: list) -> None:
        """
        Parses the gene sequence for instances of a motif. A list of motif sequences
        are passed into the function. Returns None but sets the motifs attribute of 
        the Gene object. This attribute is a dictionary where the motif sequence is
        the key and the value is a list of tuples representing the (start, end) 
        positions of each occurance of that motif. Positions are in 1-based numbering. 
        """
        self.motifs = {}
        for motif in motifList:
            motif_len = len(motif)
            motif = motif.upper()
            motif1 = motif.replace("U", "[UT]")
            motif2 = motif1.replace("Y", "[CTU]")
            motif3 = motif2.replace("R", "[AG]")
            pattern = f'(?=({motif3}))'
            motif_locations = [(m.start()+1, m.start()+motif_len) for m in re.finditer(pattern.upper(), self.seq.upper())]
            if len(motif_locations) > 0:
                self.motifs[motif] = motif_locations
        return


class GeneDraw:

    """
    Class used to create the GeneDraw object which is used to draw out the pycairo images.
    Has no attributes but uses two methods to generate the png gene image based on the data
    in the Gene objects. 
    """
    
    def findMaxSeqLen(self, genes: list) -> None:
        """
        Iterates over a list of Gene objects and finds the length of the longest 
        sequence and returns None but sets the max_seq_len attribute for the object. 
        """
        max_len = 0
        for gene in genes:
            if gene.seq_len > max_len:
                max_len = gene.seq_len
        
        self.max_seq_len = max_len

        return


    def drawGene(self, genes: list, filename:str, motifs:list): 
        #positions: list[tuple], seq_length:int):
        """
        This method will take a list of tuples where each tuple contains 
        (start, end) position of each occurence of a given motif. Seq_length 
        is the length of the gene. 

        Need to draw exons and motifs in a different color. 
        """

        # get file name prefix
        file_prefix = re.findall(r'.+(?=\.)', filename)[0]

        n_genes = len(genes) # number of genes to plot. Used for png dims. 
        n_motifs = len(motifs) # total number of motifs. Used to adjust png dims for legend. 
        png_len = n_genes * 100 + n_motifs*40 # adjusts the length of the png for number of genes. 

        # This list of R, G, B values will be used to color each motif in the drawing. 
        colors = [(0,255,0), (0,0,255), (255,255,0), (0,255,255), (255,0,255), 
                (192,192,192), (128,0,0), (128,128,0), (0,128,0), (128,0,128)]
        
        legend = {} # stores motifs and their colors for the figure legend. 

        # initialize the drawing
        with cairo.SVGSurface(f"{file_prefix}.svg", self.max_seq_len, png_len) as surface:

            # fill background
            context = cairo.Context(surface)
            context.set_source_rgb(1,1,1) # white background
            context.rectangle(0, 0, self.max_seq_len, png_len)
            context.fill()

            n = 50 # used for spacing genes. 
            for gene in genes:
                # add gene name as text
                context.set_source_rgb(0, 0, 0)
                context.set_font_size(15)
                context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
                context.move_to(0, n-20)
                context.show_text(gene.name)

                # darw 1 line per gene, representing introns. 
                context.set_source_rgb(0,0,0)
                context.set_line_width(1)
                context.move_to(0,n)        #(x,y)
                context.line_to(gene.seq_len,n)
                context.stroke()

                # draw 1 rectangle per exon for each gene. 
                for position in gene.exonPos.values():
                    context.set_source_rgb(0,0,0)
                    start = position[0]
                    end = position[1]
                    length = end-start
                    context.rectangle(start,n-20,length,40)       
                    context.stroke()

                # draw 1 rectangle per motif. 
                color = 0 # used to cycle through the colors list. 
                for motif_name in gene.motifs.keys():

                    # keep track of what color was used for each motif in legend dict.
                    if motif_name not in legend:
                        r, g, b = colors[color]
                        context.set_source_rgb(r, g, b)
                        legend[motif_name] = (r, g, b)
                    else:
                        r, g, b = legend[motif_name]
                        context.set_source_rgb(r, g, b)
                    for motif_occurance in gene.motifs[motif_name]:
                        start = motif_occurance[0]
                        end = motif_occurance[1]
                        length = end-start
                        context.rectangle(start,n-10,length,20) 
                        context.stroke()
                    color += 1

                n += 100

            # add legend
            # add legend text
            context.move_to(0, n)
            context.set_source_rgb(0, 0, 0)
            context.show_text("Motif sequence | Color")

            # add motif sequences and corresponding colors to legend. 
            for motif in legend.keys():
                n += 20
                context.set_source_rgb(0, 0, 0)
                context.set_font_size(10)
                context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
                context.move_to(0, n)
                context.show_text(motif)

                r, g, b = legend[motif]
                context.set_source_rgb(r, g, b)
                x_start_pos = len(max(legend.keys(), key=len))*13 # adjust for seq length
                context.rectangle(x_start_pos, n-10,10,10)        #(start x pos, top y pos, x length, height)
                context.stroke()

            surface.write_to_png(f"{file_prefix}.png")

        return 


def get_args():
    """
    Parse command line for arguments and return. 
    """
    parser = argparse.ArgumentParser(description="""
    Parses a FASTA file and motif sequence file and draws a line and box representation of
    introns and exons with the location of each motif. Motifs are color coded. The size of
    each intro, exon, and motif is drawn to scale. Assumes that introns are in lower case
    and exons are in upper case in the FASTA file.
    """)
    parser.add_argument("-f", help="fasta file", required=True, type=str)
    parser.add_argument("-m", help="motif file", required=True, type=str)
    return parser.parse_args()
args = get_args()


def readFiles(fasta: str, motif_file:str) -> tuple[dict,list]:
    """
    Reads in the fasta file, makes the sequence one line, and creates an 
    object based on the Gene class. 

    reads is a dictionary where the keyes are fasta headers and the 
    values are the oneline sequences. 

    returns a tuple storing reads (dictionary of headers and sequences),
    and motifs (list of motif sequences).
    """
    reads = onelineFasta(fasta)
    with open(motif_file, "r") as mf:
        motifs = [line.rstrip('\n') for line in mf]

    return reads, motifs




################# Function calls ###########################

# Parse the FASTA and motifs files. 
seq_dict, motifs = readFiles(args.f, args.m)


# instantiate a gene object for each read in the FASTA file using the Gene class. 
# store each gene object in a list. 
gene_objects = []
for read in seq_dict.keys():
    gene = Gene(read, seq_dict[read])
    gene_objects.append(gene)
    gene.getExonPos()
    gene.findMotifs(motifs)


# instantiate the plotting object. 
plot_object = GeneDraw()
plot_object.findMaxSeqLen(gene_objects)


# use the drawGene method to draw each gene and the motifs. 
plot_object.drawGene(gene_objects, args.f, motifs)


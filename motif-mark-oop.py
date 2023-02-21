#!/usr/bin/env python

import re
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", help="fasta file", required=True, type=str)
    parser.add_argument("-m", help="motif file", required=True, type=str)
    return parser.parse_args()
args = get_args()


class Gene:
    
    def __init__(self, seq:str):
        self.seq = seq

    def getExonStartPos(self):
        '''
        Parses a gene sequence, searches for the exons which are in all caps,
        and returns the 1-based start position of each exon. Returns a dictionary
        with exon number as the key and a tuple (start pos, end pos) as the value.
        '''
        exons = {}
        i = 1 # counter for current 1-based position
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
                    exons[f"Exon_{e}"] = [i, i]
            i += 1
        self.exonPos = exons
        return


    def findMotifs(self, motifList: list):
        self.motifs = {}
        for motif in motifList:
            self.motifs[motif] = [(m.start()+1, m.end()) for m in re.finditer(motif.upper(), self.seq.upper())]
        return





class Motif:
    pass



gene1 = Gene("aaaaaaaaaaGGGGGGATGCaaaaaaatgcGGGGGGGGGGaaaCCC")
gene1.getExonStartPos()

motifs = ["atgc", "cc"]
gene1.findMotifs(motifs)

print(gene1.exonPos)
print(gene1.motifs)
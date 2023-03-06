#!/usr/bin/env python

'''
Author: Christian La France clafranc@uoregon.edu
This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
'''
__version__ = "0.9"         
                            
DNA_bases = "ATGCNatgcn"
RNA_bases = "AUGCNatgcn"

def convert_phred(letter:str) -> int:
    '''
    Takes as an argument a single ASCII character and 
    returns the corresponding Phred 33 score. 
    '''
    score = ord(letter) - 33
    return score


def qual_score(phred_score: list) -> float:
    """
    Takes as an argument a list of ASCII characters, converts 
    each into the corresponding numerical score and returns 
    their average.
    """
    total = 0
    for score in phred_score:
        total += convert_phred(score)
        
    return total/len(phred_score)


def validate_base_seq(seq:str, RNA=False) -> bool:
    '''
    Takes as an argument a sequence and returns a bool, True if
    the sequence contains valid bases and False if not. An optional
    second argument specifies if the input sequence is RNA (bool).
    Default is False.
    '''

    seq = seq.upper()

    return len(seq) == seq.count("A") + seq.count("U" if RNA else "T") + seq.count("G") + seq.count("C") 



def gc_content(seq:str) -> float:
    '''
    Takes as an argument a string of DNA bases and calculates
    the GC content. 
    '''
    seq = seq.upper()

    gc_content = (seq.count("G") + seq.count("C"))/len(seq)

    return gc_content


def onelineFasta(fasta_file: str):
    '''
    Accepts as an argument a fasta file name as a string and
    converts it to a fasta file without newline characters in
    the sequences. Returns a dictionary where the keyes are 
    fasta headers and the values are the oneline sequences. 
    '''

    reads = {}
    fasta_header = None
    with open(fasta_file, "r") as fa:
        for line in fa:
            if line[0] == ">":
                fasta_header = line.strip("\n")
                reads[fasta_header] = "" # add the header
            else:
                reads[fasta_header] += line.strip("\n")

    return reads



if __name__ == "__main__":
    assert validate_base_seq("AATAGAT", False) == True
    assert validate_base_seq("AAUAGAU", True) == True
    assert validate_base_seq("TATUC",False) == False
    assert validate_base_seq("UCUGCU", False) == False

    assert gc_content("GTCA") == 0.5
    assert gc_content("ATTC") == 0.25

    assert convert_phred("E") == 36

    assert qual_score(["E", "A"]) == 34.0
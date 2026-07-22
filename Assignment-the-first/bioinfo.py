#!/usr/bin/env python

# Author: hanana2000 hannahkapoor00@gmail.com

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "6.0"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNAbases = set('ATGCNatcgn')
RNAbases = set('AUGCNaucgn')
DNACOMPDICT = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}

print(f"right now the name is {__name__}")

def validate_base_seq(seq,RNAflag=False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    if type(letter) is int: return letter
    else: return ord(letter)-33

def calc_median(lst: list) -> float:
    '''
    if odd, then return the middle value. if even, return the average of the two middle values.
    expects an already sorted list 
    '''
    listlen = len(lst)
    return lst[(listlen//2)] if listlen%2 else (lst[listlen//2 - 1]+lst[listlen//2])/2

def qual_score(phred_score: str) -> float:
    '''gets the average score of a line of phred scores'''
    sum = 0
    for pos in phred_score: sum+=convert_phred(pos)
    return sum/len(phred_score)

def gc_content(seq: str) -> float:
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    GCcount = 0
    if validate_base_seq(seq) or validate_base_seq(seq, True): 
        for base in seq: 
            if base in ('G','C','g','c'): GCcount+=1 
    else: assert False, "Not a DNA or RNA seq!"
    return GCcount/len(seq)

def oneline_fasta(filename: str, outfile: str) -> bool:
    '''
    converts a fasta file so that each sequence is on one line
    with a one line header (2 lines per entry)
    takes two arguments, the fasta file and the name for an output file 
    '''
    with open(filename, "r") as infile, open(outfile, "w") as outf: 
        record = ""
        for line in infile: 
            if line.startswith(">"):
                if record: outf.write(f"{record}\n")
                outf.write(line)
                record = ""
            else: record+=line.strip()
        outf.write(record)

def reverse_compliment(seq: str) -> str: 
    """
    return the reverse compliment of a string 
    """
    comp = ""
    if not validate_base_seq(seq): return ""
    for base in seq.upper(): 
        comp+=DNACOMPDICT[base]
    return comp[::-1]

def corrected_revcomp(seq1: str, seq2: str) -> str: 
    """
    takes two sequences with 'N's and checks if they are reverse complimented 
    returns the corrected seq1 (no N's) if they are 
    returns an empty string if they are not 
    """ 
    out, seq1, seq2 = "", seq1.upper(), seq2.upper()[::-1]
    if len(seq1) != len(seq2) or not validate_base_seq(seq1) or not validate_base_seq(seq2): return ""
    for x in range(len(seq1)): 
        if seq1[x] == 'N' and seq2[x] == 'N': return ""
        if seq1[x] == 'N' or seq2[x] == 'N': 
            out += seq1[x] if seq1[x] != 'N' else DNACOMPDICT[seq2[x]]
            continue 
        if seq1[x] != DNACOMPDICT[seq2[x]]: return ""
        out+=seq1[x]
    return out

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    print("Passed DNA and RNA tests")
    assert qual_score("@B=>1") == 27.4, "wrong phred score for '@B=>1'"
    assert qual_score("FG-#'\"#") == 14, "wrong phred score for 'FG-#'""'"
    # 31 + 33 + 28 + 29 + 16 == 137, /5 = 27.4
    # 37 + 38 + 12 + 2 + 6 + 1 + 2 = 98, /7 = 14
    print("Passed qual score tests!")
    assert gc_content("atgcgcgc") == .75, "wrong GC content for 'atgcgcgc'"
    assert gc_content("ATGcGcATat") == .4, "wrong GC content for 'ATGcGcATat'"
    print("Passed GC content tests!")
    assert calc_median([1,3,10,10,10]) == 10, "wrong median for '[1,3,10,10,10]'"
    assert calc_median([25,60,65,75]) == 62.5, "wrong median for '[25,60,65,75]'"
    print("Passed median tests!")
    assert reverse_compliment("AGCTA") == "TAGCT", "wrong reverse compliment for 'AGCTA'"
    assert reverse_compliment("AGBTF") == "", "wrong reverse compliment for 'AGBTF'"
    assert reverse_compliment("atcGaG") == "CTCGAT", "wrong reverse compliment for 'atcGaG'"
    print("passed reverse compliment tests!")
    assert corrected_revcomp("NGTANC", "GATACN") == "", "wrong corrected reverse compliment for 'NGTANC', 'GATACN'"
    assert corrected_revcomp("NGTANC", "GATNCA") == "TGTATC", "wrong corrected reverse compliment for 'NGTANC', 'GATNCA'"
    assert corrected_revcomp("AGTG", "GATCGAT") == "", "wrong corrected reverse compliment for seqs of diff length"
    print("passed corrected reverse compliment tests!")

# NGTANC
# GATNCA

# NGTANC
# ACNTAG

# TGTATC <- expected 





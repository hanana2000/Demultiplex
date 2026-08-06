#!/usr/bin/env python 

# python3 -m pip install numpy
# python3 -m pip install matplotlib
# srun --account=bgmp --partition=bgmp --cpus-per-task=8 --time=2:00:00 --pty bash

import bioinfo
import argparse
import math
import gzip
import numpy as np
import matplotlib.pyplot as plt
import itertools
import os
import re 

# ./visualizing.py -i stats_demultiplex_nocut.txt -o . -k known_barcodes.tsv -f test.png

# commands run: 
"""
./barchart.py -i stats_demultiplex_nocut.txt -o . -k known_barcodes.tsv -f nocut_chart.png
./barchart.py -i stats_demultiplex.txt -o . -k known_barcodes.tsv -f cut25_chart.png
./barchart.py -i stats_demultiplex_30.txt -o . -k known_barcodes.tsv -f cut30_chart.png

./barchart.py -i stats_demultiplex_Ncorrect_ -o . -k known_barcodes.tsv -f Ncorrect_nocut_chart.png
./barchart.py -i stats_demultiplex_Ncorrect.txt -o . -k known_barcodes.tsv -f Ncorrect_cut25_chart.png
./barchart.py -i stats_demultiplex_Ncorrect_30.txt -o . -k known_barcodes.tsv -f Ncorrect_cut30_chart.png
"""

def get_args():
    parser = argparse.ArgumentParser(description="A program to visualize the demultiplexing stats")
    parser.add_argument("-i", "--infile", help="stats file path", required=True, type=str)
    parser.add_argument("-o", "--outpath", help="path to output dir where subfolder will be created, with NO '/' at the end", required=True, type=str, default=None)
    parser.add_argument("-k", "--knownbarcs", help="path to known barcodes tsv with format of 'A12	TCGACAAG' on each line", required=True, type=str, default=None)
    parser.add_argument("-f", "--outfilename", help="name of out file", required=True, type=str, default=None)
    return parser.parse_args()
	
args = get_args()

infile, knownbarcsfile, outpath, outfilename = args.infile, args.knownbarcs, args.outpath, args.outfilename

def barcodepop(knownbarcsfile: str) -> (dict, list): 
    """
    This function takes a file path to the barcodes list file
    it will then populate a dict with all the barcodes (with the names in the files as keys, i.e. A1, C3, etc.)
    and a list with just the barcodes 
    """
    known_barcs, knownbarcs_list = {}, []
    with open(knownbarcsfile, 'r') as fh: 
        while True: 
            line = fh.readline()
            if not line: break
            key,val = line.strip().split('\t')
            known_barcs[key] = val
            knownbarcs_list.append(val)
    return known_barcs, knownbarcs_list

def pop_comparedict(barcs_list: list, infile: str) -> dict: 
    """
    This function populates the comparedict dictionary from the infile
    the keys are tuples of barcodes 
    the values are instances of the barcode pair
    """
    comparedict = {}
    with open(infile, 'r') as stats: 
        for line in stats: 
            matches = re.match(r"\('([ATGC]+)',\s'([ATGC]+)'\):\s([0-9]+)", line)
            barc1, barc2, num = "", "", ""
            if matches: 
                barc1, barc2, num = matches.group(1), matches.group(2), matches.group(3)
                comparedict[(barc1, barc2)] = num
    return comparedict

known_barcs, barcs_list = barcodepop(knownbarcsfile)
barcs_list.sort()
# for item in barcs_list: print(item)

def get_percentages(infile: str) -> dict: 
    percdict = {}
    with open(infile, 'r') as stats: 
        for line in stats: 
            matches = re.match(r"percentage of ([ATCG]+) reads: ([0-9.]+)%", line)
            barc, perc = "", ""
            if matches: 
                barc, perc = matches.group(1), float(matches.group(2))
                percdict[barc] = perc
    return percdict

percdict = get_percentages(infile)
# for k,v in percdict.items(): 
#     print(k,v) 
xvals, yvals = [x for x in percdict.keys()], [y for y in percdict.values()]

fig, ax = plt.subplots()
im = ax.bar(xvals, yvals)
ax.set_xticks(range(len(barcs_list)), labels=barcs_list,
              rotation=55, rotation_mode="xtick")
# ax.set_yticks(range(len(barcs_list)), labels=barcs_list)
# plt.bar(xvals, yvals)
plt.tight_layout()
plt.title(f"Barchart of known barc percentages for {outfilename.replace(".png","").replace("_chart","")}")
ax.set_xlabel("Barcodes")
ax.set_ylabel("Percentages")
plt.subplots_adjust(top=0.90, left=.12, bottom=.25)

plt.savefig(outfilename)


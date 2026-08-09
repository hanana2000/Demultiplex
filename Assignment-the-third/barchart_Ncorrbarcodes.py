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
./barchart_Ncorrbarcodes.py -i1 stats_txt_files/stats_demultiplex_nocut.txt -i2 stats_txt_files/stats_demultiplex_Ncorrect_nocut.txt -o . -k known_barcodes.tsv -f noN_N_nocutcompare_chart.png
"""

def get_args():
    parser = argparse.ArgumentParser(description="A program to visualize the demultiplexing stats between two runs")
    parser.add_argument("-i1", "--infile1", help="stats1 file path", required=True, type=str)
    parser.add_argument("-i2", "--infile2", help="stats2 file path", required=True, type=str)
    parser.add_argument("-o", "--outpath", help="path to output dir where subfolder will be created, with NO '/' at the end", required=True, type=str, default=None)
    parser.add_argument("-k", "--knownbarcs", help="path to known barcodes tsv with format of 'A12	TCGACAAG' on each line", required=True, type=str, default=None)
    parser.add_argument("-f", "--outfilename", help="name of out file", required=True, type=str, default=None)
    return parser.parse_args()
	
args = get_args()

infile1, infile2, knownbarcsfile, outpath, outfilename = args.infile1, args.infile2, args.knownbarcs, args.outpath, args.outfilename

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

def get_matched(sorted_comparedict: dict) -> dict: 
    matched_dict = {}
    for k,v in sorted_comparedict.items(): 
        if k[0] == k[1]: 
            matched_dict[k] = int(v)
    return matched_dict

known_barcs, barcs_list = barcodepop(knownbarcsfile)
barcs_list.sort()
# for item in barcs_list: print(item)

comparedict1, comparedict2 = pop_comparedict(barcs_list, infile1), pop_comparedict(barcs_list, infile2)
sorted_comparedict1, sorted_comparedict2 = dict(sorted(comparedict1.items())), dict(sorted(comparedict2.items()))

matched_dict1, matched_dict2 = get_matched(sorted_comparedict1), get_matched(sorted_comparedict2)

xvals = barcs_list
y_noN, y_N = [y for y in matched_dict1.values()],[y for y in matched_dict2.values()]

w, x = 0.45, np.arange(len(xvals))
fig, ax = plt.subplots()
ax.bar(x - w/2, y_noN, width=w, label='no N correction')
ax.bar(x + w/2, y_N, width=w, label='N correction')

ax.set_ylabel("Number Observed")
ax.set_xticks(x, labels=barcs_list, rotation=90, rotation_mode="xtick")
ax.set_xticklabels(xvals)
ax.set_title(f"Barchart of matched barcs for {outfilename.replace(".png","").replace("_chart","").replace("_", "-").split('/')[-1]}")
ax.legend()
plt.subplots_adjust(top=0.90, left=.12, bottom=.25)

plt.savefig(f"{outpath}/{outfilename}")


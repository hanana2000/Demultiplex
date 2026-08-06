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
./visualizing.py -i stats_demultiplex_30.txt -o . -k known_barcodes.tsv -f cut30.png
./visualizing.py -i stats_demultiplex.txt -o . -k known_barcodes.tsv -f cut25.png
./visualizing.py -i stats_demultiplex_nocut.txt -o . -k known_barcodes.tsv -f nocut.png

./visualizing.py -i stats_demultiplex_Ncorrect_30.txt -o . -k known_barcodes.tsv -f Ncorr_cut30.png
./visualizing.py -i stats_demultiplex_Ncorrect.txt -o . -k known_barcodes.tsv -f Ncorr_cut25.png
./visualizing.py -i stats_demultiplex_Ncorrect_nocut.txt -o . -k known_barcodes.tsv -f Ncorr_nocut.png
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

def pop_nump(sorted_comparedict: dict, barcs_list: list): 
    """
    this function will populate the numpy array with the sorted compare dict
    it populates with the barcodes sorted alphabetically 
    on both the x and y axis
    """
    totalarray = np.zeros((24,24),dtype=int) # create 24 lists of 24 positions 
    for k,v in sorted_comparedict.items():
        totalarray[barcs_list.index(k[0]), barcs_list.index(k[1])] = v
    return totalarray

# # Create a 2D numpy array
# data = np.random.rand(10, 10)

# # Plot the heatmap
# plt.imshow(data, cmap='viridis')
# plt.colorbar()  # Add a colorbar to show value mapping
# plt.savefig(f"{outpath}/{outfilename}")

comparedict = pop_comparedict(barcs_list, infile)
sorted_comparedict = dict(sorted(comparedict.items()))
totalarray = pop_nump(sorted_comparedict, barcs_list) # create 24 lists of 24 positions 

# print(totalarray)

# Plot the heatmap
# cmap='hot'
# cmap='viridis'
# plt.imshow(totalarray, cmap='viridis')
# plt.colorbar()  # Add a colorbar to show value mapping
# plt.title(f'2-D Heat Map for {outfilename.strip('.png')}')

fig, ax = plt.subplots()
im = ax.imshow(totalarray)
ax.set_xticks(range(len(barcs_list)), labels=barcs_list,
              rotation=55, rotation_mode="xtick")
ax.set_yticks(range(len(barcs_list)), labels=barcs_list)
# # Loop over data dimensions and create text annotations.
# for i in range(len(barcs_list)):
#     for j in range(len(barcs_list)):
#         print(round(totalarray[i, j], 1))
#         text = ax.text(j, i, round(totalarray[i, j], 1),
#                        ha="center", va="center", color="w")
ax.set_xlabel('Second Barcode')
ax.set_ylabel('First Barcode')
plt.tight_layout()
plt.colorbar(im) # Add a colorbar to show value mapping
plt.title(f'2-D Heat Map for {outfilename.replace('.png','')}')
plt.subplots_adjust(top=0.90)
plt.savefig(f"{outpath}/{outfilename}")




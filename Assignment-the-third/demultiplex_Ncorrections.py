#!/usr/bin/env python 


# python3 -m pip install numpy
# python3 -m pip install matplotlib
# srun --account=bgmp --partition=bgmp --cpus-per-task=8 --time=2:00:00 --pty bash

# command to run test files: 
# ./demultiplex.py -R1 ../TEST-input_FASTQ/R1_test.fq, -R2 ../TEST-input_FASTQ/R2_test.fq, -R3 ../TEST-input_FASTQ/R3_test.fq, -R4 ../TEST-input_FASTQ/R4_test.fq -o /projects/bgmp/hankap/bioinfo/Bi622/Demultiplex/Assignment-the-third -k ./known_barcodes.tsv

# command to run actual files: 
# ./demultiplex.py -R1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz, -R2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz, -R3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz, -R4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -o /projects/bgmp/hankap/bioinfo/Bi622/Demultiplex/Assignment-the-third -k ./known_barcodes.tsv


import bioinfo
import argparse
import math
import gzip
import numpy
import matplotlib
import itertools
import os

def get_args():
    parser = argparse.ArgumentParser(description="A program to demultiplex 4 illumina sequencing files")
    parser.add_argument("-R1", "--R1", help="R1 file path", required=True, type=str)
    parser.add_argument("-R2", "--R2", help="R2 file path", required=True, type=str)
    parser.add_argument("-R3", "--R3", help="R3 file path", required=True, type=str)
    parser.add_argument("-R4", "--R4", help="R4 file path", required=True, type=str)
    parser.add_argument("-o", "--outpath", help="path to output dir where subfolder will be created, with NO '/' at the end", required=True, type=str, default=None)
    parser.add_argument("-k", "--knownbarcs", help="path to known barcodes tsv with format of 'A12	TCGACAAG' on each line", required=True, type=str, default=None)
    parser.add_argument("-c", "--phredcutoff", help="int indicated minimum average phred score for barcodes", required=False, type=str, default=25)
    parser.add_argument("-f", "--outfilename", help="name of stats file", required=True, type=str, default=None)
    return parser.parse_args()
	
args = get_args()

R1, R2, R3, R4, knownbarcsfile, outpath, phredcutoff, outfilename = args.R1, args.R2, args.R3, args.R4, args.knownbarcs, args.outpath, int(args.phredcutoff), args.outfilename

# R1_test, R2_test, R3_test, R4_test = "../TEST-input_FASTQ/R1_test.fq","../TEST-input_FASTQ/R2_test.fq", "../TEST-input_FASTQ/R3_test.fq", "../TEST-input_FASTQ/R4_test.fq"
# knownbarcsfile = "./known_barcodes.tsv"

# R1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
# R2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
# R3 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
# R4 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"

# outputpath = "/projects/bgmp/hankap/bioinfo/Bi622/Demultiplex/Assignment-the-third/output_with_Ncorrection"
outputpath = f"{outpath}/out_test_Ncorrect"
os.makedirs(outputpath, exist_ok=True)

# R2 and R3 are barcode files 
# R1 and R4 are sequence files 
# No file names should have Ns in them 

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

# we can have all possible combos of barcodes 
# if we have ATGC and GCAA
# can have a ATGC-GCAA, GCAA-ATGC, ATGC-ATGC, GCAA-GCAA
# will never add the reverse compliment of the first barcode to this dict
# always reverse compliment the second barcode to add to this dict 

def permutepop(knownbarcs_list: list) -> dict: 
    """
    this function will output a dict of all the possible permutations of barcodes 
    permutations of barcodes will be the keys 
    and all values will be set to 0 (to be updated as we go)
    """
    permutedict = {}
    perms = itertools.product(knownbarcs_list, repeat=2)
    for item in perms: permutedict[item] = 0  
    return permutedict

def createfiledict(knownbarcs): 
    """
    This function will create a dictionary of the file objects
    with each barcode having 2 files (R1 and R2)
    as well as unknown and hopped files (R1 and R2)
    keys are the barcodes and values are the file objects 
    """
    file_obj_dict = {}
    for barc in knownbarcs: 
        R1 = open(f"{outputpath}/{barc}_R1.fq", 'w')
        R2 = open(f"{outputpath}/{barc}_R2.fq", 'w')
        file_obj_dict[barc] = (R1, R2)
    unk1, unk2 = open(f"{outputpath}/unk_R1.fq", 'w'), open(f"{outputpath}/unk_R2.fq", 'w')
    hopped1, hopped2 = open(f"{outputpath}/hopped_R1.fq", 'w'), open(f"{outputpath}/hopped_R2.fq", 'w')
    statsfile = open(f"./{outfilename}.txt", 'w')
    file_obj_dict["unk"] = (unk1, unk2)
    file_obj_dict["hopped"] = (hopped1, hopped2)
    file_obj_dict["stats"] = statsfile
    return file_obj_dict 

def writeto(file_obj_dict: dict, fileprefix: str, record1: list, record2: list) -> bool: 
    """
    this function will write a record of four lines out to the specified file 
    by looking up the fileprefix in the file_obj_dict 
    it will write to both R1 and R2
    """
    try: 
        for line in record1: 
            line = line+'\n'
            file_obj_dict[fileprefix][0].write(line)
        for line in record2: 
            line = line+'\n'
            file_obj_dict[fileprefix][1].write(line)
        return True 
    except: 
        return False 

def writestats(permutedict: dict, file_obj_dict: dict, hopped: int, known: int, unk: int, total: int) -> None: 
    file_obj_dict["stats"].write(f"total num records: {total}\nknown barcode records: {known}\nhopped barcode records: {hopped}\nunknown barcode records: {unk}\n\n")
    knowndict = {}
    for k,v in permutedict.items(): 
        file_obj_dict["stats"].write(f"{k}: {v}\n")
        if v != 0: 
            if k[0] == k[1]: 
                knowndict[k[0]] = v
    file_obj_dict["stats"].write('\n')
    for k,v in knowndict.items(): 
        file_obj_dict["stats"].write(f"percentage of {k} reads: {(v/total)*100}%\n")

def corr_barcode(barc: str, barclist: list, revcomp: bool = False) -> str:
    """
    this function will take a barcode with N in it and check if it is in the list of barcodes 
    if it is passed True for revcomp it will check if it is the rev comp of any barcode in the list 
    returns the barcode that was matched 
    """
    if not revcomp: 
        for currbar in barclist: 
            if bioinfo.corrected_comp_nlessseq2(barc, currbar): return currbar
    elif revcomp: 
        for currbar in barclist: 
            if bioinfo.corrected_revcomp_nlessseq2(barc, currbar): return currbar
    return ""
# the forward barcode (R2) should be reverse compliment of reverse (R3)
# should not see any known barcodes in the R3 file
# the 24 known barcodes can be the file names (48 files total for R1, R2)
# then 4 unk + hopped files (2 each)
# 52 total 
def demultiplex(R1, R2, R3, R4): 
    """
    This is the main function that will execute the demultiplexing by calling all other functions  
    
    """
    known_barcs, knownbarcs_list = barcodepop(knownbarcsfile)
    permutedict = permutepop(knownbarcs_list)
    file_obj_dict = createfiledict(knownbarcs_list)
    
    # file_obj_dict["unk"][1].write("hello world!")

    with gzip.open(R2, 'rt') as bar1, gzip.open(R3, 'rt') as bar2, gzip.open(R1, 'rt') as seq1, gzip.open(R4, 'rt') as seq2: # opening all the read files and fixed name write files
    # with open(R2, 'r') as bar1, open(R3, 'r') as bar2, open(R1, 'r') as seq1, open(R4, 'r') as seq2: # opening all the read files
        hopped, known, unk, total = 0,0,0,0
        currbar1, currbar2, currseq1, currseq2 = [],[],[],[] #four lists for current R2, R3, R1, and R4, will have max 4 items at any time


        while True: 
            # readline for R2, R3, R1, and R4 and add to current record lists 
            # add barcodes to ends of header lines, remember to strip!
            currbar1.append(bar1.readline().strip())
            currbar2.append(bar2.readline().strip())
            currseq1.append(seq1.readline().strip())
            currseq2.append(seq2.readline().strip())

            if currbar1[0] == "": 
                print("end of file!")
                break 
                        
            if len(currbar1) == 4: 
                total += 1
                forbar, revbar = currbar1[1], currbar2[1]
                # add barcodes to end of header line 
                currseq1[0] += f" {forbar}-{revbar}"
                currseq2[0] += f" {forbar}-{revbar}"

                # check quality cuttoff first for barcode (set by histogram results)
                if bioinfo.qual_score(currbar1[3]) < phredcutoff or bioinfo.qual_score(currbar2[3]) < phredcutoff: 
                    # write seq to R1_unk.fastq and R2_unk.fastq
                    writeto(file_obj_dict, "unk", currseq1, currseq2)
                    unk+=1

                # check if has N
                elif 'N' in forbar or 'N' in revbar: 
                    # check if it has more than 2 Ns (2 Ns or less is ok)
                    # if so, is unk 
                    if forbar.count('N') > 2 or revbar.count('N') > 2: 
                        writeto(file_obj_dict, "unk", currseq1, currseq2)
                        unk+=1
                    else: 
                        # correct the forward and reverse barcodes if contain N 
                        # the reverse barcodes will now be unreversed and complimented 
                        corrbar1, corrbar2 = "" if 'N' in forbar else forbar,"" if 'N' in revbar else bioinfo.reverse_compliment(revbar)
                        # check if the barcode had no N's but does not match a known barcode 
                        if corrbar1 and corrbar1 not in knownbarcs_list: 
                            writeto(file_obj_dict, "unk", currseq1, currseq2)
                            unk+=1
                        elif corrbar2 and corrbar2 not in knownbarcs_list: 
                            writeto(file_obj_dict, "unk", currseq1, currseq2)
                            unk+=1
                        else: # barcodes without N's that match known barcodes, and barcodes with N's get to this point 
                            # now correct the barcodes that had N's 
                            if not corrbar1: corrbar1 = corr_barcode(forbar, knownbarcs_list)
                            if not corrbar2: corrbar2 = corr_barcode(revbar, knownbarcs_list, True)
                            
                            # check if a barcode is still blank, meaning it did not match 
                            if not corrbar1 or not corrbar2: 
                                writeto(file_obj_dict, "unk", currseq1, currseq2)
                                unk+=1
                            # check if barcodes match
                            elif corrbar1 == corrbar2: 
                                writeto(file_obj_dict, corrbar1, currseq1, currseq2)
                                known+=1
                                permutedict[(corrbar1,corrbar2)]+=1
                            # check if hopped 
                            else: 
                                writeto(file_obj_dict, "hopped", currseq1, currseq2)
                                hopped+=1
                        

                # check if barcodes are reverse compliment 
                elif forbar == bioinfo.reverse_compliment(revbar) and forbar in knownbarcs_list:
                    writeto(file_obj_dict, forbar, currseq1, currseq2)
                    known+=1
                    permutedict[(forbar,forbar)]+=1


                # check if they are hopped 
                elif forbar in knownbarcs_list and bioinfo.reverse_compliment(revbar) in knownbarcs_list: 
                    writeto(file_obj_dict, "hopped", currseq1, currseq2)
                    hopped+=1
                    permutedict[(forbar,bioinfo.reverse_compliment(revbar))]+=1

                else: 
                    writeto(file_obj_dict, "unk", currseq1, currseq2)
                    unk+=1
                #     write seq out to R1_unk.fastq and R2_unk.fastq 
                #     unk+=1
                currbar1, currbar2, currseq1, currseq2 = [],[],[],[] # we are at the end of a record 

        # write out hopped, known, unk to tab separated file
        writestats(permutedict, file_obj_dict, hopped, known, unk, total)
        
        # close all the files 
        for key,file_lst in file_obj_dict.items():
            if isinstance(file_lst, tuple): 
                for file in file_lst: file.close()
            else: file_lst.close()


demultiplex(R1, R2, R3, R4)


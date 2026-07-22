#!/usr/bin/env python 

# Generate a per base distribution of quality scores for read1, read2, index1, and index2. 
# Average the quality scores at each position for all reads and generate a per nucleotide 
# mean distribution 

import bioinfo
import gzip

R1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
R2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
R3 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
R4 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"

# seqfiles = [R1, R4]
# for seqfile in files: 
#     sumlist = [0 for x in range(101)]
#     with gzip.open(file, 'wb') as f:
#         pass
        
with gzip.open(R1, 'rb') as f1, gzip.open(R2, 'rb') as f2, gzip.open(R3, 'rb') as f3, gzip.open(R4, 'rb') as f4: 
    for x in range(4): 
        line1, line2, line3, line4 = f1.readline().strip().decode('utf-8'), f2.readline().strip().decode('utf-8'), f3.readline().strip().decode('utf-8'), f4.readline().strip().decode('utf-8')
        if x == 1: 
            print(line2, line3)
            print(bioinfo.corrected_revcomp(line2, line3))
            print(bioinfo.reverse_compliment(bioinfo.corrected_revcomp(line2, line3)))
        if x == 3: 
            print(bioinfo.qual_score(line1))
            print(bioinfo.qual_score(line2))
            print(bioinfo.qual_score(line3))
            print(bioinfo.qual_score(line4))






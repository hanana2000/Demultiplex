#!/usr/bin/env python 

#SBATCH --account=bgmp                    # REQUIRED: which account to use
#SBATCH --partition=bgmp                  # REQUIRED: which partition to use
#SBATCH --cpus-per-task=8                 # optional: number of cpus, default is 1
#SBATCH --job-name=spades_k77             # optional: job name
#SBATCH --time=3:00:00                    # optional: time before timesout 

# Generate a per base distribution of quality scores for read1, read2, index1, and index2. 
# Average the quality scores at each position for all reads and generate a per nucleotide 
# mean distribution 
import sys 
import os 

sys.path.append(os.getcwd())

import bioinfo
import gzip
import matplotlib.pyplot as plt

R1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
R2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
R3 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
R4 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"

FILES = [R1, R2, R3, R4]

R1_test, R2_test, R3_test, R4_test = "../TEST-input_FASTQ/R1_test.fq","../TEST-input_FASTQ/R2_test.fq", "../TEST-input_FASTQ/R3_test.fq", "../TEST-input_FASTQ/R4_test.fq"

TESTFILES = [R1_test, R2_test, R3_test, R4_test]

def init_list(lst: list, length: int=101, value: float=0.0) -> list:
    """
    This function takes an empty list and will populate it with
    the value passed with length passed 
    """
    return [value for x in range(length)]

def populate_list(file: str, length: int=101) -> tuple[list, int]:
    """
    populates a list with running sums of phred scores for each base pair
    """
    scoresums: list = []
    scoresums = init_list(scoresums, length)

    with gzip.open(file, "rt", encoding='utf-8') as fh: 
    # with open(file, "r") as fh: 
        i = 0
        for line in fh:
            line = line.strip('\n')
            if i%4 == 3:
                for baseindex,base in enumerate(line): 
                    scoresums[baseindex]+=bioinfo.convert_phred(base)
            i+=1
    return scoresums,i

file = FILES[0] 
my_list, num_lines = populate_list(file)
print("num lines:", num_lines)
print(len(my_list))
for i,total in enumerate(my_list):
    my_list[i]=total/(num_lines/4) # turn the list of sums into a list of averages 
plt.bar([x for x in range(len(my_list))], my_list)
plt.title(f"Mean Quality Score for Base Positions of R{1}")
plt.xlabel("Base Position")
plt.ylabel("Mean Quality Score")
plt.savefig(f"distribution_R{1}.png")
plt.close()


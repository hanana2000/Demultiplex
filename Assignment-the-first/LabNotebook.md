# Assignment the First

1. Determine which files contain the indexes, and which contain the paired end reads containing the biological data of interest. Create a table and label each file with either read1, read2, index1, or index2.  
```bash 
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz | head -4

@K00337:83:HJKJNBBXX:8:1101:1265:1191 1:N:0:1
GNCTGGCATTCCCAGAGACATCAGTACCCAGTTGGTTCAGACAGTTCCTCTATTGGTTGACAAGGTCTTCATTTCTAGTGATATCAACACGGTGTCTACAA
+
A#A-<FJJJ<JJJJJJJJJJJJJJJJJFJJJJFFJJFJJJAJJJJ-AJJJJJJJFFJJJJJJFFA-7<AJJJFFAJJJJJF<F--JJJJJJF-A-F7JJJJ

zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | head -4

@K00337:83:HJKJNBBXX:8:1101:1265:1191 2:N:0:1
NCTTCGAC
+
#AA<FJJJ

zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | head -4

@K00337:83:HJKJNBBXX:8:1101:1265:1191 3:N:0:1
NTCGAAGA
+
#AAAAJJF

zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz | head -4

@K00337:83:HJKJNBBXX:8:1101:1265:1191 4:N:0:1
NTTTTGATTTACCTTTCAGCCAATGAGAAGGCCGTTCATGCAGACTTTTTTAATGATTTTGAAGACCTTTTTGATGATGATGATGTCCAGTGAGGCCTCCC
+
#AAFAFJJ-----F---7-<FA-F<AFFA-JJJ77<FJFJFJJJJJJJJJJAFJFFAJJJJJJJJFJF7-AFFJJ7F7JFJJFJ7FFF--A<A7<-A-7--

```

files 1 and 4 seem to contain biological reads, and files 2 and 3 are barcodes. 

2. Determine the length of the reads in each file.

```bash 
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz | head -2 | tail -1 | wc -c 
102

# but there is a newline char so it is 101 for the biological reads

zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | head -2 | tail -1 | wc -c
9

# without the newline the barcodes are 8 long  


```

3. Determine the phred encoding for these data.

```bash 
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz | head -100 | grep "#" | head -4
A#A-<FJJJ<JJJJJJJJJJJJJJJJJFJJJJFFJJFJJJAJJJJ-AJJJJJJJFFJJJJJJFFA-7<AJJJFFAJJJJJF<F--JJJJJJF-A-F7JJJJ
A#AAFJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJAJJJJJJJJJJJJJJFJJJJJFFFFJJJJJJJJJJJJJJJJJJ77F
A#AFFFJFJJFJJJJFJJJJJJJJAJJFJJJJJFJFJ7<FAFJJFJFJJFJFJJJFJAAJJJFJJJJJJJJJJJJJJJAJJJFAJJJJJFFJJJAJJJ<F-
A#<AAFJFJJJJFJJFJJ7JFJJJFJFAJJ<FF<<JJ<JJ<F<JJFAJJFFFJJJJJJA--77FJ--<<-AA<<AFJJJJJJFJJJFFFJ-<7--7-FFFA
```
There are #'s in the quality score lines, which would mean this is phred 33. 



notes for bash script: 

```bash
#SBATCH --account=bgmp                    # REQUIRED: which account to use
#SBATCH --partition=bgmp                  # REQUIRED: which partition to use
#SBATCH --cpus-per-task=8                 # optional: number of cpus, default is 1
#SBATCH --job-name=spades_k77             # optional: job name
#SBATCH --time=3:00:00                    # optional: time before timesout 

# scancel <jobid>
# sbatch <file name>
# squeue -u hankap
# history | tail -10
```
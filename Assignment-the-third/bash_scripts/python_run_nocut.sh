#!/bin/bash

#SBATCH --account=bgmp                    # REQUIRED: which account to use
#SBATCH --partition=bgmp                  # REQUIRED: which partition to use
#SBATCH --cpus-per-task=8                 # optional: number of cpus, default is 1
#SBATCH --job-name=demultiplex             # optional: job name
#SBATCH --time=3:00:00                    # optional: time before timesout 

# chmod 755 demultiplex.py
# /usr/bin/time -v ./demultiplex.py \
#     -R1 ../TEST-input_FASTQ/R1_test.fq \
#     -R2 ../TEST-input_FASTQ/R2_test.fq \
#     -R3 ../TEST-input_FASTQ/R3_test.fq \
#     -R4 ../TEST-input_FASTQ/R4_test.fq \
#     -o /projects/bgmp/hankap/bioinfo/Bi622/Demultiplex/Assignment-the-third \
#     -k ./known_barcodes.tsv

chmod 755 demultiplex.py
/usr/bin/time -v ./demultiplex.py \
    -R1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz \
    -R2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz \
    -R3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz \
    -R4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz \
    -o /scratch/bgmp/hankap/demux \
    -k ./known_barcodes.tsv \
    -c -500 \
    -f stats_demultiplex_nocut


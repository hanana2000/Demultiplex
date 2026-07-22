# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here:

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read 1 | 101 | Phred 33 |
| 1294_S1_L008_R2_001.fastq.gz | index 1 | 8 | Phred 33 |
| 1294_S1_L008_R3_001.fastq.gz | index 2 | 8 | Phred 33 |
| 1294_S1_L008_R4_001.fastq.gz | read 2 | 101 | Phred 33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. **YOUR ANSWER HERE**
    3. **YOUR ANSWER HERE**
    
## Part 2
1. Define the problem

There are 4 files, R1, 2, 3, and 4. R1 and R4 contain biological sequences, and R2 and R3 have indices. R1 sequences are paired with R2 indices, and R3 indeces are paired with R4 sequences. All 4 sequence names will match regardless, and are in the same order in each file. 

The problem is that most of the time, the indices will match with each other (majority of them will always be paired with the same other index, and they should be reverse complimentary to each other). However, sometimes the barcodes "hop" and do not match the reverse compliment of the paired index. 

There is also an issue of the barcodes containing 'N's, not matching our list of known barcodes, or ambiguity/ low quality, in which case they are not reliable. 

2. Describe output

The output should be R1 and R2 FASTQ files with sequences that had matching barcodes, 
two FASTQ files with non-matching barcodes (index hopping), with pairs of barcodes split into the two files, 
and two more FASTQ files with barcodes that were low quality or do not match the list of 24 known indices. 

3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
```
I just used the first entry in all 4 provided files 

ave quality scores for these lines: 

R1: 36.366336633663366
R2: 31.625
R3: 31.125
R4: 31.07920792079208

I have not yet chosen the phred score cutoff, but for this test I will assume they are above the cutoff. 

test barcodes: 
NCTTCGAC NTCGAAGA

corrected reverse compliment: 
TCTTCGAC

this is a known barcode, so it would go into a R1 R2 file
```
[input R1](../TEST-input_FASTQ/R1_test.fq)
[input R2](../TEST-input_FASTQ/R2_test.fq)
[input R3](../TEST-input_FASTQ/R3_test.fq)
[input R4](../TEST-input_FASTQ/R4_test.fq)

[output R1](../TEST-input_FASTQ/TCTTCGAC_GTCGAAGA_R1.fq)
[output R2](../TEST-input_FASTQ/TCTTCGAC_GTCGAAGA_R2.fq)  
[output hopped R1](../TEST-input_FASTQ/testout_hopped_R1.fq)
[output hopped R2](../TEST-input_FASTQ/testout_hopped_R2.fq)  
[output unk R1](../TEST-input_FASTQ/testout_unk_R1.fq)
[output unk R1](../TEST-input_FASTQ/testout_unk_R1.fq)


4. Pseudocode

```python 

# R2 and R3 are barcode files 
# R1 and R4 are sequence files
# No file names should have Ns in them 

list of known barcodes = []
dict of {barcodepair: 0, barcodepair: 0, unk: 0} # itertools to create all possible combinations of barcodes

with open (R2), with open (R3), with open(R1), with open(R4): 
    current record = [],[],[],[] #four lists for current R2, R3, R1, and R4, will have max 4 items at any time 
    while True: 
        readline for R2, R3, R1, and R4 and add to current record lists 
        add barcodes to ends of header lines 
        if len(current record) == 4: # we are at the end of a record 
            if bioinfo.qualscore(qual line of barcode) is above cutoff: # check quality cuttoff first for barcode (set by histogram results)
                if either barcode has 'N': # check if has N, then check if corrected compliment and in known barcodes
                    if corrected_revcomp(barcode1, barcode2) and in known barcodes: write seq to {barcodes}_R1.fastq and {barcodes}_R2.fastq and update dict
                    else: write seq to R1_unk.fastq and R2_unk.fastq and update dict
                # check if barcodes are reverse compliment 
                elif barcode1 == reverse_compliment(barcode2) and in known barcodes: write seq to {barcodes}_R1.fastq and {barcodes}_R2.fastq and update dict
                # check if they are hopped 
                elif barcode1 and barcode2 in known barcodes (but NOT revcomp): write seq to hopped_R1.fastq and hopped_R2.fastq and update dict
                else: write seq out to R1_unk.fastq and R2_unk.fastq and update dict
            else: # quality score is not above cutoff 
                write seq out to R1_unk.fastq and R2_unk.fastq and update dict
        if it is the end of the file (empty string for all four files): 
            break

```

5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement

''' 

```python
 
def reverse_compliment(seq: str) -> str: 
    """
    return the reverse compliment of a string 
    """
    comp = ""
    return comp

def corrected_revcomp(seq1: str, seq2: str) -> str: 
    """
    takes two sequences with 'N's and checks if they are reverse complimented 
    returns the corrected rev comp if they are 
    returns an empty string if they are not 
    """ 
    out = ""
    # iterate through each seq and if all the non-N bases match, it is reverse compliment
    return out 
     
```

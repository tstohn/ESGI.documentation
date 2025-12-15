---
title: Getting Started
tags: 
 - jekyll
 - github
description: Getting started with Docsy Jekyll
---

# Application


## Multimodal

Application example for multimodal data using **SIGNALseq**, which contains sequencing reads for both RNA and protein generated through combinatorial indexing. In this method, protein measurements are obtained using DNA barcoded antibodies that target the protein of interest. The DNA barcodes are sequenced alongside with the RNA, enabling simultaneous profiling of gene expression and protein abundances. 

Each modality has its own separate FASTQ files for both the forward and reverse reads. See below the instructions to download the HELA data:
```
---
wget "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR28056729/SRR28056729" -O SRR28056729.sra # Transcripts
wget "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR28056728/SRR28056728" -O SRR28056728.sra # Proteins

fastq-dump --split-files --gzip SRR28056729.sra
fastq-dump --split-files --gzip SRR28056728.sra
---
```

For each modality, we need to define the barcode structure, with it's possible barcode patterns for each position. Both modalities use a ten-pattern barcode design, where the first pattern encodes the feature identity. In the protein modality, this feature identifying sequence is DNA, whereas in the RNA modality is it RNA. The remaining nine barcode patterns are DNA sequences in both modalities, where the BC barcode patterns encode for the single-cell identity with either constant or random base sequences separating them. 

In both modalities, the forward read includes the feature encoding barcode, while the reverse read extends up to this barcode. In the protein modality, a polyA tail follows the feature-encoding barcode, causing overlap between the forward and reverse read. In contrast, in the RNA modality the forward read terminates at the end of the RNA sequence, and there is no overlap with the reverse read. 

Below are the barcode structures for the protein and RNA modalities, shown as ten bracket-enclosed sequence substrings. Each of the brackets corresponds to a specific position in the barcode and contains a comma-separates list of possible barcode sequences for that position. 
```
---
# Barcode structure for protein modality:
------------> <------------------------------ 
[Antibodies][*][BC1][22X][BC2][30X][BC2][10X]
|           |  |    └─── Sequence for twenty two random bases      
|           |  └─── List of possible sequences
|           └─── End of forward and reverse reads with overlap
└─── Antibody identity sequence

# Barcode structure for RNA modality:
-----> <----------------------------------------------------------------------------
[RNA][-][BC1][CCACAGTCTCAAGCACGTGGAT][BC2][AGTCGTACGCCGATGCGAAACATCGGCCAC][BC2][10X]
|    |       └─── Constant barcode
|    └─── End of forward and reverse reads without overlap            
└─── Sequence of the transcripts
---
```

For each barcode pattern, the allowed number of mismatches must be defined. The protein modality allows one mismatch in the feature-encoding barcode, whereas the RNA modality allows none. Both modalities accept one mismatch in each of the single-cell encoding barcodes (BC) and zero mismatches in the constant and random base sequences. 
```
---
# Mismatches protein modality
1,0,1,0,1,0,1,0

# Mismatches RNA modality
0,0,1,0,1,0,1,0
---
```




Set up required paths to run demultiplex and count for protein and RNA: 

```
---
# Set paths to raw data, tools, and analyses
Path_data = "/path/to/raw_data"
Path_tool= "/path/to/tools"
Path_analyses = "/path/to/analyses"

# Define log file
LOGFILE_AB="${"Path_analyses"}/output/ESGI_Protein/ESGI_PROTEIN_LOG.txt"
# Remove old log file 
rm -f $LOGFILE

#RUN WITH 1MM in SC-BARCODE and 1MM in AB-BARCODE: results have prefix A_
/usr/bin/time -v "${Path_tool"/bin/demultiplex" \
                -i "${Path_data}/SRR28056728_1.fastq.gz" \
                -r "${Path_data}/SRR28056728_2.fastq.gz" \
                -o "${Path_analyses}/output/ESGI_Protein" \
                -p "${Path_analyses}/background_data/ESGI_files/pattern_PROTEIN.txt" \
                -m "${Path_analyses}/background_data/ESGI_files/mismatches_PROTEIN_1MM.txt" \
                -n A \
                -t 70 -f 1 -q 1 2>> $LOGFILE_AB

/usr/bin/time -v "${Path_tool}/bin/count" \
                -i "${Path_analyses}/output/ESGI_Protein/A_PROTEIN.tsv" \
                -o "${Path_analyses}/output/ESGI_Protein/A_PROTEIN_Counts.tsv" \
                -t 70 -d "${Path_analyses}/background_data/ESGI_files" \
                -a "${Path_analyses}/background_data/ESGI_files/antibody_names_as_in_KITE.txt" \
                -c 1,3,5 -x 0 -u 6 -m 0 -s 1 2>> $LOGFILE_AB
---
```

###

## Spatial

## Multipattern

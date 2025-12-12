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

For each modality, we need to define the barcode structure, with it's possible barcode patterns for each location. Both modalities use a ten-pattern barcode desing, where the first pattern encoses the feature identity. In the protein modality, this feature identifying sequence is DNA, whereas in the RNA modality is it RNA. The remaining nine barcode patterns are DNA sequences in both modalities. 

```
---
# Barcode structure for protein modality:
------------> <------------------------------ 
[Antibodies][*][BC1][22X][BC2][30X][BC2][10X]
|           |  |    └─── Sequence for twenty two random bases      
|           |  └─── List of possible sequences
|           └─── Read stops here but forward and reverse do overlap
└─── Antibody identity sequence

# Barcode structure for RNA modality:
-----> <----------------------------------------------------------------------------
[RNA][-][BC1][CCACAGTCTCAAGCACGTGGAT][BC2][AGTCGTACGCCGATGCGAAACATCGGCCAC][BC2][10X]
|     |      └─── Constant barcode
|     └─── Hard stop for forward and reverse read without overlap            
└─── Sequence of the transcripts
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

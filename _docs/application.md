---
title: Application
tags: 
 - Multimodal
 - Multipattern
 - Spatial
description: Application examples for different types of data.
---

# Application


## Multimodal

Application example for multimodal data using **SIGNALseq**, which contains sequencing reads for both RNA and protein generated through combinatorial indexing. In this method, protein measurements are obtained using DNA barcoded antibodies that target the protein of interest. These DNA barcodes are sequenced alongside with the RNA, enabling simultaneous profiling of protein abundances and gene expression. 

Each modality has its own separate FASTQ files for both the forward and reverse reads. See below the instructions to download the HELA data:
```
wget "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR28056729/SRR28056729" -O SRR28056729.sra # Transcripts
wget "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR28056728/SRR28056728" -O SRR28056728.sra # Proteins

fastq-dump --split-files --gzip SRR28056729.sra
fastq-dump --split-files --gzip SRR28056728.sra
```
For the RNA-modality, specifically, you need to download the reference genome as described in [Reference Genome](software-overview#demultiplex)

To execute **ESGI** on both modalities, we must create two separate ESGI-initialization files. Both modalities use a barcode pattern consisting of eight positional elements, though they differ in how the forward read encodes feature identity. 

#### Forward read: positional element 0
The forward read contains the first pattern element and encodes the feature identity. 
* **Protein modality**: The first pattern element is a variable DNA barcode. A polyA tail follows this barcode, resulting in an overlap betweem the forward and reverse reads. 
* **RNA modality**: The first pattern element is a transcript sequence. The forward read terminates at the end of the RNA sequence, resulting in no overlap with the reverse read.

#### Read transition: positional element 1
* Protein modality uses `*` to indicate overlapping sequences.
* RNA modality uses `-` to indicate no overlap.

#### Reverse read: positional elements 2-7
For both modalities, the reverse read contains the remaining six pattern elements. These DNA-sequences are structured as following: 

| Element Index | Type | Encoding |
| --------- | ----------- | ------ | 
| 2, 4, 6 | Variable barcode element | Single-cell identities
| 3, 5 | Constant pattern element | anchors or linkers
|  7 |  Random pattern element | Unique Molecule Identifier (UMI)


Below an illustration of the barcode patterns for the two modalities shown as eight bracket-enclosed sequence substrings. Each bracket corresponds to a pattern element containing a comma-separated list of possible barcode sequences for that position. 
>```
>          ---------------> <----------------------------------------
>PROTEIN:[Antibodies.txt][*][BC1.txt][22X][BC2.txt][30X[BC2.txt][10X]
>
>      ----> <-----------------------------------------
>RNA:[RNA][-][BC1.txt][22X][BC2.txt][30X][BC2.txt][10X]
>
>```

| File name | Element length (bases) | Number of subsequences | Sequence example |
| --------- | ----------- | ------ | ------ |  
| Antibodies.txt | 15 | 23 | AAGGCAGACGGTGCA,GGCTGCGCACCGCCT,CGTCCTAGGACATAT 
| BC1.txt | 8 | 96 | TTACGAGT,TATCGTTT,CGAGGTAA
| BC2.txt | 8 | 96 | ATCACGTT,CGATGTTT,TTAGGCAT

Example of a barcode-aligned sequence for the protein modality: 
>```
>[Antibodies.txt] [BC1.txt]           [22X]          [BC2.txt]              [30X]              [BC2.txt]   [10X]  
>AGACAGTGATGTCCG  CCGATCCC   ATCCACGTGCTTGAGACTGTGG  TTAGGCAT  GTGGCCGATGTTTCGCATCGGCGTACGACT  TAACGCTG  TAAAGGAAGT
>```

Next, for each pattern element, we define the allowed number of mismatches. The protein modality allows one mismatch in the feature-encoding barcode, whereas the RNA modality allows none. Both modalities accept one mismatch for the variable barcode elements and zero for the constant and pattern elements mismatches in the constant and random base sequences. 

```
# Mismatches protein modality
1,0,1,0,1,0,1,0

# Mismatches RNA modality
0,0,1,0,1,0,1,0
```

Now we have all information to create the ESGI-initialization files:
```
Path_data = "/path/to/raw_data"
# Includes FASTQ files of forward and reverse reads

Path_background_data = "/path/to/background_data"
# .txt files for barcode pattern, mismatches, annotation

Path_output = "/path/to/output"

# Forward and reverse reads:
forward="${Path_data}/SRR28056728_1.fastq.gz"
reverse="${Path_data}/SRR28056728_2.fastq.gz"

pattern="${Path_background_data}/pattern_PROTEIN.txt"
mismatches="${Path_background_data}/mismatches_PROTEIN_1MM.txt"

# Indexing for elements encoding: feature, single-cell ID and UMI:
FEATURE_ID=0
SC_ID=2,4,6
UMI_ID=7

FEATURE_NAMES=/home/n.vd.brug/Projects/Analysis-EZGI/SIGNALseq_Analysis/background_data/ESGI_files/antibody_names_as_in_KITE.txt
ANNOTATION_IDs=/home/n.vd.brug/Projects/Analysis-EZGI/SIGNALseq_Analysis/background_data/ESGI_files/BC1.txt

threads=10
prefix=MYEXPERIMENT
```

Run the workflow for the protein modality:
**demultiplex** → **count** 
```
---
LOGFILE_AB="${"Path_output"}/ESGI_PROTEIN_LOG.txt"
rm -f $LOGFILE

/usr/bin/time -v "${Path_tool"/bin/**demultiplex**" \
                -i "${Path_data}/SRR28056728_1.fastq.gz" \
                -r "${Path_data}/SRR28056728_2.fastq.gz" \
                -o "${Path_output}/ESGI_Protein" \
                -p "${Path_background_data}/pattern_PROTEIN.txt" \
                -m "${Path_background_data}/mismatches_PROTEIN.txt" \
                -n A \
                -t 70 -f 1 -q 1 2>> $LOGFILE_AB

/usr/bin/time -v "${Path_tool}/bin/**count**" \
                -i "${Path_output}/A_PROTEIN.tsv" \
                -o "${Path_output}/A_PROTEIN_Counts.tsv" \
                -t 70 -d "${Path_background_data}/background_data/ESGI_files" \
                -a "${Path_background_data}/background_data/ESGI_files/antibody_names_as_in_KITE.txt" \
                -c 2,4,6 -x 0 -u 7 -m 0 -s 1 2>> $LOGFILE_AB

---
```

Run the workflow for the RNA modality:
**demultiplex** → **STAR** → **annotate** → **count** 
```
---
LOGFILE="${"Path_output"}/ESGI_RNA_LOG.txt"
rm $LOGFILE

/usr/bin/time -v "${Path_data}/SRR28056729_1.fastq" \
              -r "${Path_data}/SRR28056729_2.fastq" \
              -o "${Path_output}/ESGI_RNA" \
              -p "${Path_background_data}/pattern_RNA.txt" \
              -m "${Path_background_data}/mismatches_RNA.txt" \
              -t 70 -f 1 -q 1 2>> $LOGFILE
              
/usr/bin/time -v STAR --runThreadN 70 \
     --genomeDir "${Path_reference_genome}/GRCh38/GRCh38_STAR_index" \
     --readFilesIn "${Path_output}/RNA.fastq" \
     --outFileNamePrefix "${Path_output}/RNA_" \
     --sjdbGTFfile "${Path_reference_genome}/GRCh38/gencode.v43.annotation.gtf" \
     --sjdbOverhang 73 \
     --outSAMtype BAM Unsorted \
     --outSAMattributes NH HI AS nM GX GN \
     --quantMode TranscriptomeSAM \
     --outFilterMultimapNmax 50 \
     --outSAMmultNmax 1 --outSAMunmapped Within \
     --limitOutSJcollapsed 2000000 \
     --twopassMode Basic 2>> $LOGFILE
     
/usr/bin/time -v "${Path_tool}/bin/annotate \
              -i "${Path_output}/ESGI_RNA/RNA.tsv" \
              -b "${Path_output}/RNA_Aligned.out.bam" \
              -f GX 2>> $LOGFILE

/usr/bin/time -v "${Path_tool}/bin/count
              -i "${Path_output}/RNA_annotated.tsv" \
              -o "${Path_output}/RNA_Counts_umi0.tsv" -t 70 \
              -d "${Path_background_data}/" \
              -c 2,4,6 -x 0 -u 7 -m 1 -s 1 \
              -w "${Path_background_data}/bc_sharing_revComp.tsv" 2>> $LOGFILE
---
```

###

## Spatial

## Multipattern

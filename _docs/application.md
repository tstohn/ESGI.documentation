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

Application example for multimodal data using **SIGNALseq**, which produces separate sequencing libraries for both RNA and protein modalities via combinatorial indexing. In this method, protein levels are captured using DNA barcoded antibodies and transcripts through mRNA reverse transcription. Both modalities are processed into distinct FASTQ files, but share a common barcode pattern: the forward read encodes the feature identity and the reverse read encodes the single-cell identity and Unique Molecular Identifier (UMI). 

Download the SRA files for both modalities and extract the paired-end FASTQ reads:
```
# Download SRA files for Protein and Transcript
wget "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR28056728/SRR28056728" -O SRR28056728.sra
wget "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR28056729/SRR28056729" -O SRR28056729.sra

# Extract paired-end FASTQ files
fastq-dump --split-files --gzip SRR28056728.sra
fastq-dump --split-files --gzip SRR28056729.sra
```

### Set up:
To process both modalities, separate ESGI-initialization files must be created. Both modalities share a pattern structure consisting of eight positional elements, but they differ in how the feature identity is encoded in the forward read.

#### Forward read: positional element 0
The forward read contains the first pattern element and encodes the feature identity. 
* **Protein modality**: The first element is a variable antibody derived DNA barcode. This is followed by a poly-A tail, which creates a sequence overlap between the forward and reverse reads. 
* **RNA modality**: The first element is the genomic transcript sequence. The forward read terminates at the end of the RNA sequence, resulting in no overlap with the reverse read.

#### Read transition: positional element 1
A specific symbol is used to define the relationship between the forward and reverse reads:
* **Protein modality:** uses the `*` symbol to indicate overlapping sequences.
* **RNA modality:** uses the `-` symbol to indicate no overlap.

#### Reverse read: positional elements 2-7
For both modalities, the reverse read contains the remaining six pattern elements encoding the single-cell identity and UMI.

| Element Index | Type | Encoding |
| --------- | ----------- | ------ | 
| 2, 4, 6 | Variable barcode element | Single-cell identities
| 3, 5 | Constant pattern element | anchors or linkers
|  7 |  Random pattern element | Unique Molecule Identifier (UMI)

Below an illustration of the patterns for both modalities shown as eight bracket-enclosed sequence substrings. Each bracket corresponds to a pattern element containing a comma-separated list of possible barcode sequences for that position. 
          
>```
>PROTEIN:[Antibodies.txt][*][BC1.txt][22X][BC2.txt][30X[BC2.txt][10X]
>```
       
>```
>RNA:[RNA][-][BC1.txt][22X][BC2.txt][30X][BC2.txt][10X]
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

Next, for each pattern element, we define the allowed number of mismatches. The protein modality allows one mismatch in the feature-encoding barcode, whereas the RNA modality allows none. Both modalities accept one mismatch for the variable barcode elements and zero for the constant and random pattern elements. 

Mismatches protein modality:
>```
>1,0,1,0,1,0,1,0
>```

Mismatches RNA modality
>```
>0,0,1,0,1,0,1,0
>```

Now, we have defined all structural parameters to create the ESGI-initialization files for both modalities:`myExperiment_PROTEIN.ini` and `myExperiment_RNA.ini`.

Both configuration files use the element indexes for feature identities, single-cell IDs, and UMIs from the table, along with the patterns and mismatch information illustrated in the example blocks to define `patterns.txt` and `mismatches.txt`.

`myExperiment_PROTEIN.ini`:
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

FEATURE_NAMES="${Path_background_data}/antibody_names.txt"
ANNOTATION_IDs="${Path_background_data}/BC1.txt"

threads=10
prefix=MYEXPERIMENT
```

For the RNA modality, genomic sequences are aligned to the human reference genome (GRCh38) using the **STAR** aligner. 

First, download the GRCh38 primary assembly and the corresponding GENCODE annotation files.
>```
>mkdir GRCh38
>cd GRCh38
>
>wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
>wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz
>gunzip *.gz
>```

Next, generate the STAR genome index: 
>```
>STAR --runThreadN 70 \
>     --runMode genomeGenerate \
>     --genomeDir GRCh38_STAR_index \
>     --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
>     --sjdbGTFfile gencode.v43.annotation.gtf \
>     --sjdbOverhang 73
>```

To integrate the reference genome with **ESGI**, the initialization file requires an additional parameter defining the directory path to the STAR genome index. 

`myExperiment_RNA.ini`:
```
Path_data = "/path/to/raw_data"
# Includes FASTQ files of forward and reverse reads

Path_background_data = "/path/to/background_data"
# .txt files for barcode pattern, mismatches, annotation

Path_output = "/path/to/output"

# Forward and reverse reads:
forward="${Path_data}/SRR28056729_1.fastq.gz"
reverse="${Path_data}/SRR28056729_2.fastq.gz"

pattern="${Path_background_data}/pattern_RNA.txt"
mismatches="${Path_background_data}/mismatches_RNA_1MM.txt"

# Indexing for elements encoding: feature, single-cell ID and UMI:
FEATURE_ID=0
SC_ID=2,4,6
UMI_ID=7

genomeDir="/path/to/GRCh38_STAR_index"/

threads=10
prefix=MYEXPERIMENT
```

Once the configuration files are ready, you can initiate the process for each modality by running **ESGI** from the terminal.

For the protein modality, execute the following command:
```
./bin/esgi myExperiment_PROTEIN.ini
```

For the RNA modality, run this command:
```
./bin/esgi myExperiment_RNA.ini
```

## Multipattern
Application example for multipattern data, using **scIDseq**. This technology quantifies intracellular protein abundances using antibodies conjugated to unique DNA barcodes. The approach uses a multipattern design where linkers of varying lengths, in combination with variable barcode elements, encode the specific protein identity. The remaining of the pattern encodes the single cell identity and UMI.

### Set up:
**ESGI** can demultiplex reads for multiple barcode patterns simultaneously. This specific dataset includes eight barcode patterns defined by different linker lengths, ranging from one to eight bases. To maintain the same total sequence length, the final element in the six-element patterns varies inversely with the linker length. The exception is the longest linker length pattern, which contains only five elements. 

The barcode patterns are contained entirely within the forward read and consist of five to six positional elements. The table below outlines the  element indexes and what they encode for. 

| Element Index | Type | Encoding |
| --------- | ----------- | ------ | 
| 0,3,(5) | Constant pattern element | Linker or anchor
| 1 | Random pattern element | Unique Molecule Identifier (UMI)
| 2 | Variable pattern element | Feature identities
| 4 | Variable pattern element | Well-plate position

The multipattern design contains specific and shared barcodes elements:
* Feature identity (element 2): each pattern has a unique linker, length and sequence, and is associated with a unique set of barcode sequences that define the protein identity. 
* Well plate position (element 4), universal barcode set to ensure consistent assignent of well-plate positions.  

Each pattern is assigned a unique name, with its pattern elements represented as a series of bracket enclosed substrings. Within each set of brackets is a comma-separated list of all possible barcode sequences for that specific position. Accordingly, patterns 1 through 7 are defined by six bracketed substrings, while pattern 8 consists of five. 

>```
>---
>PATTERN_1:[][][][][][]
>PATTERN_2:[][][][][][]
>PATTERN_3:[][][][][][]
>PATTERN_4:[][][][][][]
>PATTERN_5:[][][][][][]
>PATTERN_6:[][][][][][]
>PATTERN_7:[][][][][][]
>PATTERN_8:[][][][][]
>---
>```

For each pattern, the maximum number of allowed mismatches per element is defined using a comma-separated list. Each integer in the list corresponds to a specific positional element. Patterns 1 throug 7 require six mismatch values, while pattern 8 requires five. 

>```
>0,1,1,1,1,2
>0,1,1,1,1,2
>1,1,1,1,1,1
>1,1,1,1,1,1
>1,1,1,1,1,1
>1,1,1,1,1,1
>1,1,1,1,1,1
>1,1,1,1,1
>```

With all structural patameters defined, you can now create the ESGI-initialization file, `myExperiment.ini`. This file links your raw forward reads to the pattern and mismatch definitions established in the previous sections. 

The configuration uses the element indexes for feature identities, single-cell IDs, and UMIs as defined in the table, alongside the pattern and mismatch information as illustrates in the example blocks to define the `patterns.txt` and `mismatches.txt` files. 

```
Path_data = "/path/to/raw_data"
# Includes FASTQ files of forward reads for all plates

Path_background_data = "/path/to/background_data"
# .txt files for barcode pattern, mismatches, annotation

Path_output = "/path/to/output"

# Forward read:
forward="${Path_data}/plate.fastq.gz"

pattern="${Path_background_data}/patterns.txt"
mismatches="${Path_background_data}/mismatches.txt"

# Indexing for elements encoding: feature, single-cell ID and UMI:
FEATURE_ID=2
SC_ID=4
UMI_ID=1

threads=10
prefix=MYEXPERIMENT
```

Execute **ESGI** by running the following command in your terminal:  
```
./bin/esgi myExperiment.ini
```

## Spatial
Application example for spatial data using **Multiplexed Deterministic Barcoding in Tissue (xDBiT)**. This technology uses microfluidic-based deterministic barcoding with DNA oligonucleotides to encode transcriptomes alongside their spatial coordinates for multiple tissue sections in parallel. 

### Set up:
The **xDBiT** barcode pattern consists of eight positional pattern elements. The data is generated as two independent forward reads. The first forward read captures the transcript and the second forward read encodes the UMI and spatial (x,y) coordinates. 

#### Forward read I: positional element 0
The read contains the first pattern element, a genomic DNA sequence encoding the transcript. 

#### Read transition: positional element 1
A discrete transition separating the two forward reads without sequence overlap. By including the `[-]` symbol in the pattern and enabling the `independent` flag, the tool treats both reads as two distinct sequences in the `5'→3'` direction.

#### Forward read II: positional elements 2-6
The read contains the remaining six pattern elements, which collectively encode the (x,y) spatial coordinates and UMI. 

| Element Index | Type | Encoding |
| --------- | ----------- | ------ | 
| 0 | Genomic sequence | Transcript identity
| 2 | Random sequence element of 10 bases | UMI
| 3,5,7 | Variable pattern element | (x,y) spatial coordinates
| 4,6 | Constant sequence element of 30 bases | Linkers or anchors

The barcode pattern is represented as eight bracket-enclosed substrings. Each bracket identifies a specific positional element and contains a comma-separated list of possible barcode sequences for that position. 
          
>```
>SPATIAL:[DNA][-][10X][coordinate_barcode.txt][30X][coordinate_barcode.txt][30X][coordinate_barcode.txt]
>```

The `coordinate_barcode.txt` file defines the (x,y) spatial coordinates using an 8x12 matrix of 96 unique 8-base long barcode sequences.

Maximum allowed mismatched for each of the eight pattern elements:
>```
>0,0,0,1,0,1,0,1
>```

## Refences
1 SIGNALseq
2 scIDseq
3 xDBiT

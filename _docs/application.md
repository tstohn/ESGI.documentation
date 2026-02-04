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

Also, download the GRCh38 human reference genome and annotation files as described in the [STAR](software-overview#demultiplex) software overview.

### Set up:
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

          ---------------> <----------------------------------------
          
>```
>PROTEIN:[Antibodies.txt][*][BC1.txt][22X][BC2.txt][30X[BC2.txt][10X]
>```
>
       ----> <-----------------------------------------
       
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

Next, for each pattern element, we define the allowed number of mismatches. The protein modality allows one mismatch in the feature-encoding barcode, whereas the RNA modality allows none. Both modalities accept one mismatch for the variable barcode elements and zero for the constant and pattern elements mismatches in the constant and random base sequences. 

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

While both files follow a similar structure, the RNA modality requires an additional parameter to define the directory path to the reference genome. 

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
Application example for multipattern data, using **single cell Immuno-Detection by sequencing (scIDseq)** data. This technology quantifies intracellular protein abundances using antibodies conjugated to unique DNA barcodes. The approach uses a multipattern design where linkers of varying lengths, in combination with variable barcode elements, encode the specific protein identity. 

### Set up:
**ESGI** can demultiplex reads for multiple barcode patterns simultaneously. This specific dataset includes eight barcode patterns defined by different linker lengths, ranging from one to eight bases. To maintain the same total sequence length, the final element in the six-element patterns varies inversely with the linker length. The exception is the longest linker length pattern, which contains only five elements. 

The barcode patterns are contained entirely within the forward read and consist of five to six positional elements. The table below outlines the  element indexes and what they encode for. 

| Element Index | Type | Encoding |
| --------- | ----------- | ------ | 
| 0,3,(5) | Constant pattern element | Linker or anchor
| 1 | Random pattern element | Unique Molecule Identifier (UMI)
| 2 | Variable pattern element | Feature identities
| 4 | Variable pattern element | Single-cell identities

Across the patterns, there are specific and shared elements: 
* Variable pattern element 2, encoding the feature identity, consists of a defined list of barcode sequences unique to each linker.
* Variable pattern element 4, encoding the single-cell identity, is shared across all patterns to assign well-plate positions. 

Each pattern is defined by a unique name, with its pattern elements represented as a series of bracket enclosed substring. Within each set of brackets is a comma-separated list of all possible barcode sequences for that specific position. Accordingly, patterns 1 through 7 are defined by six bracketed substrings, while pattern 8 consists of five. 

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


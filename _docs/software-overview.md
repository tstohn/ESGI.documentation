---
title: Versioning your Docs
description: How to turn on and use versioning
---

# Software Overview

**ESGI** is a debarcoding tool for single-cell sequencing data, consisting of two submodules: **demultiplex** and **count**. They can be executed directly via **ESGI** or called individually. For more details, run demultiplex, count, or ESGI with the  `--help` flag. 

## Demultiplex
Demultiplexing maps sequencing reads to barcode patterns, where positional pattern elements are used to encode experiment-specific information, like cell identities, molecular modalities, and experimental conditions. The tool handles simultaneous mapping to multiple barcode patterns, supports pattern elements of varying lengths, and allows for mismatches in the pattern elements arising from insertions, deletions, and substitutions. 

### Input:
To demultiplex barcode sequences, provide the experimental FASTQ files along with text files defining the barcode patterns and elements. Also, you must specify the number of allowed mismatches for each pattern element. 

Required input parameters:

| Option | Description | Type |
| --------- | ----------- | ------ | 
| `--input`, `-i` | Single-end or forward read file | fastq(.gz) 
| `--reverse`, `-r` | Reverse read file (optional) | fastq(.gz)
|  `--output`, `-o` |  Output directory | Directory 
| `--BarcodePatternsFile`, `-p` | Description of the barcode patterns, specifying the pattern name followed by its positional pattern elements in bracket-enclosed sequence substrings. Each bracket contains a comma separated list of possible barcode sequences for that position, and these barcodes may vary in length | (.txt)
|  `--mismatchFile`, `-m` |  A comma-separated list of integers, one for each pattern in the barcode scheme, specifying the number of mismatches allowed for in each bracket-enclosed substring | (.txt) 

Example of a barcode pattern called PATTERN consisting of five positional pattern elements:

>```
>         -----------------> <-----------------------
>PATTERN:[Ab_barcodes.txt][-][BC1.txt][AGCTCATC][10X]
>         |                |  |        |         └─── Random pattern element
>         |                |  |        └─── Constant pattern element
>         |                |  └─── Variable barcode element encoding single-cell ID
>         |                └─── Read separator for forward and reverse read
>         └─── Variable barcode element encoding feature identity
>```

Variable barcode elements contain a list of possible barcode sequences for that position
>```
>CTGATC,CGTTGA,GTAGCG,CATCGT
>```

Define the number of allowed mismatched per pattern element as a comma-separated list.
>``` 
>1,0,1,1,0
>```

Optional parameters:

| Option | Description | Default |
| --------- | ----------- | ------ | 
| `--independent`, `-d` | Treat the forward read as two separate sequences in the 5'-->3' direction. Use together with the read separator [-] in the `--BarcodePatternsFile` to indicate where one read ends and the other begins.  | Disabled |
|  `--namePrefix`, `-n` | Prefix for file names. | None |
|  `--writeStats`, `-q` | Creates output files for statistics of pattern element alignment. | Disabled |
|  `--writeFailedLined`, `-f` | Generate output file for reads with failed pattern alignment. | 0 |
|  `--hamming`, `-H` | Use Hamming distance instead of Levenshtein for variable barcodes. Only supported when allowing for one mismatch per pattern element. | Levenshtein |


### Output: 
When **demultiplex** completes, it generates output files containing demultiplexed reads and quality metrics. The output type and subsequent downstream steps depend on the feature modality encoded by the pattern element. 
* **Generic barcode sequence**: Reads are split and aligned to the barcode pattern and saved as a TSV file.
* **Genomic sequence**, Gene or transcript sequences are extracted into a FASTQ file for alignment with **STAR** and subsequent annotation. Remaining read elements are aligned to the pattern and saved in a TSV file.

| Name | Description | File type | 
| ---------- | ------------- | ----------- |
| `Quality_numberMM`| Reports the frequency of mismatches (0, 1, 2) per pattern element. |TXT 
| `Quality_typeMM` | Counts per pattern element the according type of mismatch; insertion, deletion or subsitution. | TXT
| `FailedLines` | Contains all reads that were not aligned to a pattern | TXT
| `PATTERN`  | For every barcode pattern, a TSV file is created with split and aligned reads. The first column contains the read name, followed by the pattern elements. | TSV 
| `Genomic Sequence` | Extracted gene of transcript sequences | FASTQ 

Also, it reports the alignment performance into three match type categories, expressed as percentages of the total sequencing reads:
* **Perfect matches:** Reads whose barcodes match completely
* **Moderate matches:** Reads that match within the allowed number of mismatches
* **Mismatches:** Reads that cannot be matched given the number of allowed mismatches

{% include alert.html type="info" title="Generic barcode sequences follow demultiplex → count, whereas genomic sequences follow demultiplex → STAR → annotate → count." %}

### STAR
Following demultiplexing, the genomic sequences within the FASTQ files are mapped to a reference genome using the **STAR** (Spliced Transcripts Alignment to a Reference) aligner. 

Instructions for downloading reference genomes - like for human (GRCh38) or mice (GRCm38) - along with instructions for generating the required STAR genome indices, are provided in the [Multimodal](application#multimodal) and [Spatial](application#spatial) application sections, respectively. 

Once the genome index is constructed, STAR can perform the allignment to generate the output files for annotation. 

| File type | Description |
| --------- | ----------- |
| aligned.out.bam | Binary file containing the sequence reads aligned to the reference genome coordinates.
| TSV | Tab separated table of detected splice junctions and exon-intron boundaries.

### Annotate
The **annotate** step uses the **STAR**-derived genomic coordinates to assign gene annotations to the reads and merges this information with the final demultiplexed TSV output file. 

{% include alert.html type="info" title="Finally, both workflows output one TSV file per barcode pattern; containing all barcode-element-aligned reads." %}

## Count
The **count** submodule groups the demultiplexed reads by single-cell and feature barcode. By default, identical UMI-tagged entries are collapsed to produce the final counts for each unique cell- feature combination. 

### Input:
As input, **count** takes barcode-aligned reads (TSV-file) along with the barcode information, including barcode scheme, patterns and allowed mismatches. It also requires indices of the barcode patterns encoding the single-cell, feature, and UMI identities. Counting starts at zero and includes every pattern element defined within square brackets. Note when defining these indexes, the TSV output file of **demultiplex** lists the read name in de first column, followed by the positional barcode patterns. 

Required parameters:

| Option | Description | Type |
| --------- | ----------- | ------ | 
| `--input`, `-i` | Barcode-aligned reads for features in single cells |TSV
| `--barcodeDir`, `-d` | Directory containing all barcode-related files, including barcode scheme, patterns, and mismatches | Directory
| `--FeatureIndex`, `-x` | Indices referring to the positions of the barcode patterns that encode feature identities | Comma-separated list
| `--singleCellIndices`, `-c` | Indices referring to the positions of the barcode patterns that encode the single-cell identities | Comma-separated list
| `--umiIndex`, `-u` | Indices referring to the positions of the barcode patterns that encode the UMI | Comma-separated list

Optionally, feature names and single-cell annotations, like experimental conditions, can be added by providing them in the same order as corresponding barcodes in the barcode pattern files. 

Optional parameters for annotation:

| Option | Description | File type | Default |
| ---------- | ------------- | ----------- | -------------- |
| `--featureNames`, `-a` | List of all feature names, in same order as the feature-barcodes in the barcode file | (.txt) | Nucleotide sequences
|  `--annotationFiles`, `-g` | List of paths to files, where every file contains a list of single-cell annotations | (.txt) | None
|  `--annotationIdxs`, `-y` | List of space separated indices used to annotate cells  | (.txt) | None
|  `--scIdAsString`, `-s` | Stores the single-cell ID as the sequence string instead of barcode ID | None | 1

By default, UMI collapsing is enabled using Hamming distance, allowing a single mismatch. Using the parameters below, you can disable UMI collapsing or change the distance measure and the number of allowed mismatches.

Optional parameters for UMI collapsing:

| Option | Description | Default |
| ---------- | ------------- | -------------- |
|  `--mismatches`, `-m` | Number of mismatches allowed during UMI collapsing | 1
|  `--hamming`, `-H` | Use Hamming instead of Levenshtein distance during UMI collapsing (allowing for insertions and deletions) | 0
|  `--umiThreshold`, `-f` | Threshold to filter UMIs before collapsing | 0
|  `--umiRemoval`, `-z` | UMIs are collapsed | 1

### Output
After **count** completes, it creates several output files. The most important are the collapsed and uncollapsed count matrices. 
- **UMI uncollapsed** : Reports the total number of reads per UMI, providing insight into UMI amplification during the PCR experiment.
- **UMI collapsed**: All reads sharing the same UMI are collapsed to create the final count matrix for each unique cell- feature combination.

| File name | Description | Type |
| ---------- | ------------- | -------------- |
|  `LOG` | A summary of total processed reads and the number of detected UMI mismatches.  | tsv
|  `COUNTDATA` | The final countmatrix with counts for unique cell-feature combinations (UMI collapsed). | tsv
|  `UMIDATA` | Matrix with counts for unique cell-feature-UMI combinations (UMI uncollapsed). | tsv
|  `UMISTAT` | Describes UMI amplification per feature and its occurrence | tsv

## ESGI
Instead of running **demultiplex** and **count** separately, you can execute them together using **ESGI**. To do this, you must provide an *initialization-file (.ini)* containing the essential input, including paths to all barcode files, as well as indexing information. Including feature names and single-cell annotations is optional.

Example of ESGI Initialization-file:

>```
># Input directory
>forward="/path/to/forward_reads.fastq.gz"
># Reverse reads (optional)
>reverse="/path/to/reverse_reads.fastq.gz"
>
># Output directory
>output="/path/to/output"
>
># Paths to barcode scheme and mismatch settings
>pattern="/path/to/barcode_pattern.txt"
>mismatches="/path/to/mismatches.txt"
>
># Indexes for single-cell, feature, and UMI assignment
>SC_ID=1,5
>FEATURE_ID=3
>UMI_ID=4
>
># Optional: feature names and single-cell annotations
>FEATURE_NAMES="/path/to/feature_names.txt"
>
>ANNOTATION_IDs="/path/to/annotation_ID.txt"
>ANNOTATION_NAMES="/path/to/annotation_names.txt"
>
>threads=10
>prefix=MYEXPERIMENT
>```

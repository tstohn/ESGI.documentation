---
title: Versioning your Docs
description: How to turn on and use versioning
---

# Software Overview

**ESGI** is a flexible debarcoding tool for single-cell sequencing data. It can be split into two separate steps: demultiplexing the reads and counting the single-cell features. ESGI can be run as a wrapper of these tools, or **demultiplex** and **count** can be executed separately for more control. 

During the demultiplexing step, the FASTQ reads are split into the elements as defined in the input pattern. Elements are corrected when they fall within the allowed mismatch threshold. If genomic sequences are included, ESGI can call **STAR** to map these sequences to a reference genome and add the resulting annotation to the demultiplexed reads. The counting step collapses UMIs, generates a single-cell feature matrix, and assigns additional barcoded attributes, like treatment conditions. Alongside the final matrix, ESGI can report qulaity control metrics for both the demultiplexing and counting step.

For more details, run demultiplex, count, or ESGI with the  `--help` flag. 

## Demultiplex
Demultiplex maps sequencing reads to barcode patterns, where positional pattern elements are used to encode experiment-specific information, like single cell identities, molecular modalities, and experimental conditions. The tool handles simultaneous mapping to multiple barcode patterns, supports pattern elements of varying lengths, and allows for mismatches in the pattern elements arising from insertions, deletions, and substitutions. 

### Input:
To perform demultiplexing, input yout experimental FASTQ files and define barcode patterns en elements using text files. Also, specify the allowed mismatches for each pattern element as a comma-separated list of integers (e.g. `0,1,1`). 

**Required input parameters:**
        
| Option | Description | Type |
| --------- | ----------- | ------ |
| `--input`, `-i` | Path to the single-end or forward read FASTQ file. | fastq(.gz)
| `--reverse`, `-r` | Path to the reverse read FASTQ file (optional). | fastq(.gz)
|  `--output`, `-o` |  Directory for all output files. | Directory
| `--BarcodePatternsFile`, `-p` | A text file defining barcode patterns. Format: patternname followed by bracketed `[A,B,C]` lists of barcodes for each positional element. | (.txt)
|  `--mismatchFile`, `-m` |  A text file containing a comma-separated list of integers. Each integer defines the allowed number of mismatches for that bracketed element in the patterns file. | (.txt)

**Pattern element options:**
* Barcodes `[filename.txt]`: Path to text file containing a list of barcode sequences. 
* Constant region `[ATCG]`: Fixed nucleotide sequence.
* UMI `[nX]`: A unique molecule identifier consisting of `n` random bases. 
* Forward/reverse separator `[-]`: Use when there is no overlap between forward and reverse reads and you want them to be mapped independently.
* Stop mapping element `[*]`: stops mapping at this point. This allows you to ignore long, error-prone, or non-informative regions.

The pattern expects the forward read to be in the `5′→3′` direction, while the revere read corresponds to the reverse complement of the pattern. Forward and reverse reads may overlap. If the forward and reverse reads are non-overlapping, the read separator `[-]` can be applied together with the `independent` flag, as described in the optional parameters table below. 

*Example of a barcode pattern called `PATTERN` consisting of five positional pattern elements:*

```
        -----------------> <-----------------------
PATTERN:[Ab_barcodes.txt][-][BC1.txt][AGCTCATC][10X]
         |                |  |        |         └─── Random element (ten bases)
         |                |  |        └─── Constant element
         |                |  └─── Barcode element encoding single-cell ID
         |                └─── Read separator for forward and reverse read
         └─── Barcode element encoding feature identity
```

For barcode elements, define the possible sequences for that specific positional element as a comma-separated list.

>```
>CTGATC,CGTTGA,GTAGCG,CATCGT
>```

Specify the maximum allowed mismatched for each pattern element using a comma-separated list of integers.
>``` 
>1,0,1,1,0
>```

**Optional parameters:**

| Option | Description | Default |
| --------- | ----------- | ------ |
| `--independent`, `-d` | Treats the forward read as two distinct sequences in the `5′→3′` direction. Requires the separator `[-]` in the `--BarcodePatternsFile` to indicate where one read ends and the other begins.  | Disabled
|  `--namePrefix`, `-n` | Sets a prefix for all output file names. | None
|  `--writeStats`, `-q` | Generates statistics files regarding the alignment of individual pattern elements. | Disabled
|  `--writeFailedLined`, `-f` | Generates an output file containing reads that failed to align with any pattern. | 0
|  `--hamming`, `-H` | Uses Hamming distance instead of Levenshtein for variable barcodes. Only supported when allowing for one mismatch per pattern element. | Levenshtein

### Output: 

For every pattern, **demultiplex** outputs a tab-delimited (TSV) file containing the demultiplexed reads. In this file, every successfully demutliplexed read appears as a row split and aligned according to the elements specified in the pattern. If the pattern includes a DNA element, an additional FASTQ file is generated alongside the TSV file, containing the extracted genomic sequences. These sequences can later be aligned using **STAR**. In addition, several quality metrics are generated. 

**Output files:**

| Name | Description | File type |
| ---------- | ------------- | ----------- |
| `Quality_numberMM`| Reports the frequency of mismatches (0, 1, 2) detected per pattern element. |TXT
| `Quality_typeMM` | Counts of mismatch types (insertion, deletion, subtitution) per pattern element. | TXT
| `FailedLines` | Contains all reads that failed to aligned to a pattern | TXT
| `PATTERN`  | A TSV file per pattern where reads are split and aligned by element. The first column is the read name, followed by the aligned elements. | TSV
| `Genomic Sequence` | Extracted gene or transcript sequences | FASTQ

The tool also reports the alignment performance into three match type categories, expressed as percentages of the total sequencing reads:
* **Perfect matches:** Reads whose barcodes match completely.
* **Moderate matches:** Reads that match within the maximum allowed mismatches.
* **Mismatches:** Reads that cannot be matched given the maximum allowed mismatches.

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
The **count** module takes as input the pattern-aligned TSV file generated during demultiplexing, along with indexing information. It requires indices of the pattern elements encoding the single-cell, feature, and UMI identities. 

Indices follow zero-based counting and correspond to the columns in the pattern-aligned TSV file. In this file, the first column (index 0) is always the read name, folowed by the pattern elements (starting at index 1).

**Required parameters:**

| Option | Description | Type |
| --------- | ----------- | ------ |
| `--input`, `-i` | The TSV file of pattern-aligned reads generated by the demultiplex module. |TSV
| `--barcodeDir`, `-d` | Path to the directory containing all pattern-related files. | Directory
| `--FeatureIndex`, `-x` | Index of the pattern element encoding the feature identity. | Comma-separated list
| `--singleCellIndices`, `-c` | Indices of the pattern elements encoding the single-cell identities. | Comma-separated list
| `--umiIndex`, `-u` | Index of pattern element encodeding the UMI. It can be a single element or a list of elements. If it is a list all elements are concatenated and treated as a single long UMI. This feature is only available when running count individually. | Comma-separated list

**Optional parameters for annotation:**

Barcode sequences within a pattern element can be replaced by labels, such as feature names or single-cell annotations (e.g. experimental conditions). To implement this, provide an annotation file where the entries follow the exact order of the barcode sequences defined in the pattern element.

*Example:*

Barcode sequences for a pattern element encoding feature identity:
>```
>CTGATC,CGTTGA,GTAGCG,CATCGT
>```

Feature names in same order as corresponding barcode sequences:
>```
>feature1,feature2,feature3,feature4
>```

| Option | Description | File type | Default |
| ---------- | ------------- | ----------- | -------------- |
| `--featureNames`, `-a` | A list of feature names provided in same order as the barcode sequences in the pattern element. | (.txt) | Nucleotide sequences
|  `--annotationIdxs`, `-y` | A space separated list of indices for the pattern elements you want to annotate. | (.txt) | None
|  `--annotationFiles`, `-g` | A space-separated list of file path containing the annotations. The first file maps to the first index and contains labels for the corresponding barcode sequences.  | (.txt) | None
|  `--scIdAsString`, `-s` | Stores the single-cell ID as the sequence string instead of barcode ID. | None | 1

**Optional parameters for UMI collapsing:**

By default, UMI collapsing is enabled using Hamming distance, allowing a single mismatch. Using the parameters below, you can adjust the distance metric, mismatch threshold, or disable collapsing entirely.

| Option | Description | Default |
| ---------- | ------------- | -------------- |
|  `--mismatches`, `-m` | Maximum number of mismatches allowed to consider two UMIs identical. | 1
|  `--hamming`, `-H` | If set to 0, use Levenshtein distance, allowing for insertions and deletions. If 1, use Hamming distance, allowing for substitutions only. | 0
|  `--umiThreshold`, `-f` | Minimum read count required for a UMI to be considered valid before collapsing. | 0
|  `--umiRemoval`, `-z` | UMIs are collapsed. If set to 0, every unique read is counted separately. | 1

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
Instead of running **demultiplex** and **count** separately, you can execute them together using **ESGI**. To do this, you must provide an `initialization-file (.ini)` containing the essential input, including paths to all pattern files, as well as indexing information. Including feature names and single-cell annotations is optional.

*Example of ESGI Initialization-file:*

```
# Input directory
forward="/path/to/forward_reads.fastq.gz"
# Reverse reads (optional)
reverse="/path/to/reverse_reads.fastq.gz"

# Output directory
output="/path/to/output"

# Paths to barcode scheme and mismatch settings
pattern="/path/to/barcode_pattern.txt"
mismatches="/path/to/mismatches.txt"

# Indexes for single-cell, feature, and UMI assignment
SC_ID=1,5
FEATURE_ID=3
UMI_ID=4

# Optional: feature names and single-cell annotations
FEATURE_NAMES="/path/to/feature_names.txt"

ANNOTATION_IDs="/path/to/annotation_ID.txt"
ANNOTATION_NAMES="/path/to/annotation_names.txt"

threads=10
prefix=MYEXPERIMENT
```

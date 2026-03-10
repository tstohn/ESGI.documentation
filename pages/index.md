---
layout: page
title: ESGI
permalink: /
---

# Demultiplexing barcoded single-cell sequencing data with ESGI

ESGI - Efficient Splitting of Generic Indices - makes demultiplexing single-cell data easy. Download the tool from [Release](https://github.com/tstohn/ESGI/releases/tag/v1.0.0) and start demultiplexing your own data. 

ESGI enables flexible demultiplexing and counting of various barcoding-based single-cell data types, including multimodal, combinatorial and spatial formats. As input, ESGI uses user-defined barcode patterns as a template to processes raw FASTQ files. After demultiplexing, the tool counts single-cell features and collapses UMIs to generate a final single-cell feature matrix. In addition to this final data matrix, ESGI generates quality metrics for botch the demultiplexing and count steps. 

![assets/img/overview.png](assets/img/overview.png)

## Features

ESGI works on any multiplexed single-cell data. Its generic input pattern allows you to demultiplex almost any barcoded single-cell data.
Here are some features and use-cases that ESGI can handle:

 - ESGI can demultiplex barcoded data like 10X, Combinatorial-indexing experiments, spatial sequencing data, etc.
 - ESGI takes a pattern as input, which encodes the barcode-structure. It allows barcodes with fixed bases, variable barcodes that can be given by a file, RNA sequences or UMIs (15X for 15 nucleotide long UMI), e.g.: [GACGCATACGT][15X][SINGLE_CELL_BARCODES_FILE.txt][PROTEIN_BARCODES_FILE.txt]
 - ESGI can map the RNA part in the read, but can also count other modalities, like CITE-seq or intra-cellular protein measurement technologies.
 - For RNA mapping we recommend to run ESGI (together with a mapper like STAR) on a server
 - For simpler experiments (like counting protein tags) ou can easily run the App on your local laptop
 - ESGI allows also barcodes on different length
 - For flexability ESGI can also take several patterns at once as input, in case you pooled different modalities, or read-structures


For features, getting started with development, see the {% include doc.html name="Getting Started" path="getting-started" %} page. Would you like to request a feature or contribute?
[Open an issue]({{ site.repo }}/issues)

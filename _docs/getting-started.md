---
title: Getting Started
tags: 
 - jekyll
 - github
description: Getting started with Docsy Jekyll
---

# Getting Started


## Getting the Tool

The tool can either be downloaded or build from scratch. 

### Download Tool

The precompiled binaries are available on the GitHub releases: 
[https://github.com/tstohn/ESGI/releases/tag/v1.0.0](https://github.com/tstohn/ESGI/releases/tag/v1.0.0) 

Go to latest release, find the assets for your operating system, then download and extract the ZIP file. 

### Build Tool

To compile the tool on your own system, clone the GitHub repository, install the dependencies, and build **ESGI**. 

>```
># Clone repository
>git clone https://github.com/tstohn/ESGI.git
>cd ESGI/
>
># Install dependencies
>make install
>
># Build ESGI
>make esgi
>```

After building, the executables **demultiplex**, **count**, **annotate** and **ESGI**, will be available in the *bin/ directory*. You can run each submodule individually, or use ESGI to execute the complete workflow. 


## Reference Genome

The **STAR** alignment tool is used to map barcode patterns that contain gene or transcript sequences to a reference genome.  
Below are the instructions for downloading the reference genome and annotation files.

#### Download the Reference genome:
>```
>mkdir GRCh38
>cd GRCh38
>
>wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
>wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz
>gunzip *.gz
>```

#### Generate STAR genome index and annotation file:
>```
>STAR --runThreadN 70 \
>     --runMode genomeGenerate \
>     --genomeDir GRCh38_STAR_index \
>     --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
>     --sjdbGTFfile gencode.v43.annotation.gtf \
>     --sjdbOverhang 73
>```

STAR outputs genomic coordinates for aligned reads, which can be annotated using bedtools. To enable this, the gene annotation file must first be converted from *.gtf* format to *.bed* format.
>```
># Convert GTF to BED (filtering by gene_name)
>gtf2bed --attribute-key=gene_name < gencode.v43.annotation.gtf > genes.bed
>
># Keep only relevant columns and remove duplicates
>cut -f1-4 genes.bed | sort -u > genesFiltered.bed
>
># Create a simplified BED file with chr, start, end, and gene_id
>awk '{match($0, /gene_id "([^"]+)"/, arr); print $1"\t"$2"\t"$3"\t"arr[1];}' genes.bed > genes_simple.bed
>```

#### Next, read how to execute **ESGI** and its submodules in the [Software Overview](software-overview#demultiplex)


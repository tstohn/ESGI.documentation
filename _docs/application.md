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

| Positional element | Type | Encoding |
| --------- | ----------- | ------ | 
| 2, 4, 6 | Variable barcode element | Single-cell identities
| 3, 5 | Constant pattern element | anchors or linkers
|  7 |  Random pattern element | Unique Molecule Identifier (UMI)


Below an illustration of the barcode patterns for the two modalities shown as eight bracket-enclosed sequence substrings. Each bracket corresponds to a pattern element containing a comma-separated list of possible barcode sequences for that position. 
>```
>          ---------------> <--------------------------------------------------------------------------------------
>PROTEIN:[Antibodies.txt][*][BC1.txt][CCACAGTCTCAAGCACGTGGAT][BC2.txt][AGTCGTACGCCGATGCGAAACATCGGCCAC][BC2.txt][10X]
>
>      ----> <--------------------------------------------------------------------------------------
>RNA:[RNA][-][BC1.txt][CCACAGTCTCAAGCACGTGGAT][BC2.txt][AGTCGTACGCCGATGCGAAACATCGGCCAC][BC2.txt][10X]
>```

For each barcode pattern, the allowed number of mismatches must be defined. The protein modality allows one mismatch in the feature-encoding barcode, whereas the RNA modality allows none. Both modalities accept one mismatch in each of the single-cell encoding barcodes (BC) and zero mismatches in the constant and random base sequences. 
```
---
# Mismatches protein modality
1,0,1,0,1,0,1,0

# Mismatches RNA modality
0,0,1,0,1,0,1,0
---
```

Next, we need to specify which barcode patterns encode the single-cell identity, feature identity and Unique molecule Identifier (UMI). This is done by providing a list of indices matching to their positions in the barcode structure, with counting starting at 0. Here, indices are derived from the output *TSV* file of **demultiplex**, in which the first column is  the readname followed by a column for each barcode pattern. In both modalities, the first barcode pattern encodes the feature identity, followed by pattern 3, 5, and 7 encoding the single-cell identity, and the last pattern the UMI.
```
---
# Single-cell indices (-c)
2, 4, 6 

# Feature index (-x)
0

# UMI index (-u)
7
---
```

Now we can set up the required paths to execute the full workflow: 
```
---
Path_data = "/path/to/raw_data"
- Includes FASTQ files of forward and reverse reads

Path_background_data = "/path/to/background_data"
- .txt files for barcode pattern, mismatches, annotation

Path_STAR = "/path/to/star"
Path_reference_genome = "/path/to/reference_genome"

Path_tool= "/path/to/tools"

Path_output = "/path/to/output"
---
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

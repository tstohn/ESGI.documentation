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
>         |                |  └─── Variable barcode element encoding single-cell identity
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
For genomic sequences, demultiplexed FASTQ reads are aligned to a reference genome using the **STAR** (Spliced Transcripts Alignment to a Reference) aligner. 

Use the commands below to download the human reference genome (GRCh38) and the corresponding annotation files.

>```
>mkdir GRCh38
>cd GRCh38
>
>wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
>wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz
>gunzip *.gz
>```

Next, generate a STAR genome index: 
>```
>STAR --runThreadN 70 \
>     --runMode genomeGenerate \
>     --genomeDir GRCh38_STAR_index \
>     --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
>     --sjdbGTFfile gencode.v43.annotation.gtf \
>     --sjdbOverhang 73
>```

Once the index is build, run **STAR** for alignment to produce the output files required for annotation.

| File type | Description |
| --------- | ----------- |
| aligned.out.bam | Sequence reads aligned to the reference genome.
| TSV | Table of detected splice junctions and exon-intron boundaries.

### Annotate
The **annotate** step uses the **STAR**-derived genomic coordinates to assign gene annotations to the reads and merges this information with the final demultiplexed TSV output file. 

{% include alert.html type="info" title="Finally, both workflows output one TSV file per barcode scheme; containing all barcode-scheme-aligned reads." %}

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
After **count** completes, it creates the following output files: 

| File name | Description | Type |
| ---------- | ------------- | -------------- |
|  `LOG` | A summary of total processed reads and the number of detected UMI mismatches.  | tsv
|  `COUNTDATA` | The final countmatrix with counts for unique cell-feature combinations (UMI collapsed). | tsv
|  `UMIDATA` | Matrix with counts for unique cell-feature-UMI combinations (UMI uncollapsed). | tsv
|  `UMISTAT` | Describes UMI amplification per feature and its occurrence | tsv

- **UMI uncollapsed** : Reports the number of reads per UMI, providing insight into UMI amplification during the PCR experiment.
- **UMI collapsed**: All reads sharing the same UMI are collapsed to create the final count matrix for each unique cell- feature combination.

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









### Version Keywords
There are six special keywords that can be used in the versions config.
* main
* alpha
* beta
* rc
* pre
* *current*

The value **main**, is used to display the latest development version of the documentation, a banner is used to denote that the code and documentation maybe unstable and that they are not guaranteed to be correct or up to date and could change at any time.

The value of **alpha**, **beta**, **rc** and **pre** are used to denote that the documentation is a pre-production version and as such a banner is displayed announcing that the code/documentation maybe unstable and the content in the pages are still work in progress so may not be completely up to date or correct.

If you wish to utilize this feature then you must create folders within your _docs/*Archive* directory with the names of those you wish to use and place your documentation inside it/them.

The value of **current** is used to align with the latest version of your documentation in your base _docs folder so there **doesn't** need to be a specific folder created and named 'current' to serve them. The name 'current', for the value of 'latest' in the previous example, can be changed to whatever you wish to name the current release but be aware that the version set in the _config.yml file also needs to match it! See the example below where we have called the latest release v3.2 rather than 'current' and we also have a version exactly named v3.2 too, as mentioned this version will point to your base _docs directory. 

```yml
version_params:
  version_menu: "Release"
  version_dir: Archive
  tocversion_dir: versions
  versioning: true
  latest: v3.2
  versions:
    - main
    - v3.2
    - beta
    - v3.1
```

### Configuring Version Alerts

Depending on the version keywords you may use Docsy Jekyll tries to determine if the documentation is historical (e.g. out of date) or is a upcoming release of information and will show one of the following banners to alert users to the fact that they are not reviewing the current version of the documentation. 

|![Historical Versions]({{ site.baseurl }}/assets/img/versionalertoutdated.png)|
|:--:| 
| *This banner is shown for historical versions of your documentation* |

|![Beta Versions]({{ site.baseurl }}/assets/img/versionalertbeta.png)|
|:--:| 
| *This banner is shown for pre-production versions of your documentation* |

|![Development Versions]({{ site.baseurl }}/assets/img/versionalertmain.png)|
|:--:| 
| *This banner is shown for development versions of your documentation* |

These banners can be customised to have the wording or imagery you require by altering the HTML in the _includes/versionalert.html file.

## Table of Contents Handling (TOC)

As with the documents, your TOC contents also need to be updated so they correctly point to the document versions you are viewing.

You will need to update the toc file as the structure of your documentation/toc changes over time but as we have taken a snapshot of our documentation so we need to do the same thing for our toc so it reflects the document structure at that point in time too. You will see that for this site we have the versions of Current and Previous and you will see that the Current version has all the information about versioning which didn't exist in the prior version.

To enable the TOC to change based on the version, we have to create a new toc file and name it {version}-toc.yml and it must be placed in a subdirectory of _data whose name is specified by the tocversion_dir parameter in the _config.yml file (for this site it can be found in the _data/versions subdirectory). We have created previous-toc.yml in the _data/versions directory that has the TOC structure for the prior version where we did not have the versioning.md file present. 


To create this version of the toc file file we can use the following commands
```shell
# Our last release is version "Previous" so we want to create a TOC for that inside our versions directory
mkdir -p _data/versions
# copy the existing TOC file to our version subdirectory and name it as {version}-toc.yml
cp  _data/toc.yml _data/versions/previous-toc.yml
```


{% include alert.html type="info" title="Note: The toc.yml file will always be our current table of contents pointing to the most recent version of the documentation in the base of our _docs folder. When adding additional items to your toc.yml file you do not need to worry about the version identifier as that is taken care of by Docsy Jekyll when viewing alternative versions, continue to have the urls pointing to the files in the base docs directory." %}

Here is what our current toc.yml file in our _data directory looks like, the previous-toc.yml file will look identical apart from omitting the 'Versioning' title and link which was added for this release, and of course will live in the versions subdirectory of _data as well.

```yml
title: Documentation
  url: docs
  links:
    - title: "Getting Started"
      url: "docs/getting-started"
      children:
        - title: Features
          url: "docs/getting-started#getting-started"
        - title: Development
          url: "docs/getting-started#development"
        - title: Customization
          url: "docs/getting-started#customization"
    - title: Versioning
      url: "docs/versioning"
    - title: "About"
      url: "about"
    - title: "News"
      url: "news"
- title: "Extras"
  url: "docs/extras"
  links:
    - title: Quizzes
      url: "docs/extras/example-quiz"
    - title: Tags Page
      url: "tags"
```

Now not every version of your documentation will have a change in structure and therefore require an update to toc.yml. To remove this overhead a mapping file is used to tell Docsy Jekyll which TOC it should use for a particular version,
the benefit of this is that if the structure of your site does not change then you do not need to create an additional toc.yml file.

You will have a toc-mapping.yml file in your _data directory and the contents will contain the code below to map a particular version with a particular TOC file to use. You need to update this file for each version you create even if you do not need to create a new TOC file.

```yml
# This file can be used to explicitly map a release to a specific table-of-contents
# (TOC). You'll want to use this after any revamps to information architecture, to ensure
# that the navigation for older versions still works.

Current: toc
Previous: previous-toc
```
**You do not need to prefix the name of the new TOC filename with the folder (in our case versions) as Docsy Jekyll will retrieve that from the _config.yml file and add that to the path.**  

{% include alert.html type="info" title="Note: If a version cannot be found in this mapping file then the standard toc.yml in the _data directory will be used instead." %}

Because we added the versioning.md markdown file to this site we needed a new toc so that this new content could be accessed from the TOC menu. If we had only made an update to the content of the files that already existed, then we could have pointed Previous to toc in the yml above to use the same toc.yml file.

## Versioning and Search

As you will have seen in the [Versioning Options](#versioning-options) you have the ability to configure Docsy Jekyll to perform it's searches over the entire site's documentation, or only certain versions which you can specify. If you choose not to enable it in _config.yml (allow_search = false) then the searches that get performed will only be over the most current version of the documentation in your base _docs folder. If you prefer to have this enabled then you will need to provide a list of search versions (search_version in _config.yml) that you wish to be displayed if a match is found. The following image shows what your users would be shown if you have this enabled. You can see that the search was for the term 'start' and we have several results returned but two both for the 'Getting Started' topic. The reason for this is due to a result being returned within the current version as well as a result being returned for the previous version as denoted by the badges of 'Current' and 'Previous' below the titles. 


| *Searching with Versioning enabled*|
|:--:|
|![Search with Versions]({{ site.baseurl }}/assets/img/docsy-jekyll-version-search.png)


## Issues with Permalinks

Please do not use permalinks in your archive documentation Front Matter if you wish to use versioning as this will break the versioning system. If you have permalinks in your current documentation you will need to remove them when you place the documents into a subdirectory so that Jekyll can serve them correctly. 

{% include alert.html type="warning" title="Using permalinks in your archive documentation Front Matter will cause issues with versioning as this effectively 'hardcodes' the url of that page and will break the versioning links. Please be aware!" %}

## 404 Errors

There will be times where users will navigate around the documentation and could end up pointing to an invalid page that doesn't exist. For instance while you are reading this page on the **Current** version if you were to change this to the Previous version, via the dropdown, then you will get a page not found as this versioning page did not exist in the prior version. We have created a 404.md page in the base of the site so that the user still gets to navigate even though the page was not found. You can update this page with your own information should you so wish.

{% include alert.html type="info" title="When pages are not found and the 404 page is displayed the TOC menus default back to links pointing to the current release." %}

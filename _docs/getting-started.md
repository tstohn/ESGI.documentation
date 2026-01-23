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

##### Next, read how to execute **ESGI** and its submodules in the [Software Overview](software-overview#demultiplex)


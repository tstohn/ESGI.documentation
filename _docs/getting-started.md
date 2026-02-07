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

Go to latest release and download the ZIP file compatible with your operating system. Once extracted, you can find all executable files in the `bin/ directory`.

### Build Tool

To build **ESGI** locally, clone the GitHub repository, install the dependencies, and run the build command. 

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

Once the build is complete, the `bin/ directory` will contain the **ESGI**, **demultiplex**, **count**, and **annotate** executables. You can run the submodules individually or use ESGI to execute the full workflow. 

For detailed usage, visit the [Software Overview](software-overview#demultiplex). 

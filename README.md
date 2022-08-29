# bulk-RNAseq

## Introduction

This repo contains in-house tutorials and outside links primarily related to bulk RNA-seqencing data processing and analysis. Basic pipelines and default options are primarily used throughout the tutorials. Additional links to package and software documentation are provided in each tutorial - please be sure to review  these documents when completing your own analysis. A separate repo for sc/sn-RNAseq will be developed in the future. 

## Getting Started

The Allen Lab has a server dedicated to bioinformatics - this server can be used to store and process your sequencing data. 

### How to Access

- Mac OsX users can use the Terminal to access and command the server
- Windows users, please download [putty](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html) or setup Windows Subsystem for Linus (WSL, e.g. Ubuntu) 

1)	Sign-in using Username (@mnla1.salk.edu) and password
2) When you login for the first time you will be asked to create a new password. If you need to reset your password, contact JB
3) In your command prompt, you will use the secure command and type: ssh username (including mnla1.salk.edu)
4) Note: putty users you will select ssh client and then enter the address where prompted
5) Note: You need to be connected to the Salk VPN to access the server 
6) R/Rstudio can be downloaded and used locally on your computer or through the lab server. For help accessing or downloading R/Rstudio, feel free to reach out to JB.
7) Download [FileZilla](https://filezilla-project.org/) to access Server files without using the command line

### A few rules

- Backups of all raw data generated in the lab should be kept on the server.
- The default directory is /home/username. Space is **extremely** limited in the /home directory. Please save all data and analysis files on /data
- If you experience any issues on the lab server, please reach out to JB or LLB asap. 


## In-house tutorials/notes
> For information on setting up and using R/RStudio, please check out our other [R and Statistics: The Basics - Repo](https://github.com/MNL-A/r-statistics-basics) or [R and Statistics: The Basics - Webpage Verion](https://mnl-a.github.io/r-statistics-basics/)

| Link | Description |
| ----------- | ----------- |
|[Mapping with Salmon](bulk-RNAseq-tutorials/Salmon_DESeq2_Pipeline/RNAseq_Salmon.html)| Bash Commands and Scripts|
| [RNAseq Pipeline w/ Salmon + DESeq2](bulk-RNAseq-tutorials/Salmon_DESeq2_Pipeline/AgingAstrocyteTranscritptome_Tutorial.html) | RNA-seq analysis in R after Salmon quasi-mapping |
|[PDF: RNAseq Data Processing STAR & FeatureCounts](/RNAseq_STARFeatureCounts_Tau12mRibotagISH.pdf) | Command line notes/scripts for RNA-seq data processing with STAR and FeatureCounts using ISH data |
|[HTML: RNAseq Data Processing STAR & FeatureCounts](/RNAseq_STARFeatureCounts_Tau12mRibotagISH.html) | HTML version of the above notes | 

## Below are some additional resources that I've found particularly useful

### Linux Command Line

[Cheatsheet 1](https://phoenixnap.com/kb/linux-commands-cheat-sheet#linux-commands-cheat-sheet-pdf)

[Cheatsheet 2](https://www.guru99.com/linux-commands-cheat-sheet.html)

### Salmon for RNA-seq

[Official Salmon Docs](https://salmon.readthedocs.io/en/latest/salmon.html)

[Example workflow with good explanations](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/08_salmon.html)
 
### General RNA-seq pipelines

[Comparison of methods](https://hbctraining.github.io/Training-modules/planning_successful_rnaseq/slides/RNAseq-strategies_mm.pdf)

### R Studio

[Cheatsheets](https://www.rstudio.com/resources/cheatsheets/)

### Markdown

[Offical Markdown Guide](https://www.markdownguide.org/basic-syntax/)

#### 
Reach out to Jillybeth if you have suggestions, questions or comments related to the information provided above.


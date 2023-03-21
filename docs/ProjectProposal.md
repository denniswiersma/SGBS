# Project proposal
Dennis Wiersma - 378279  
Project proposal for final dataprocessing assignment.  
Academic year 2022-2023

Pipeline to reproduce [Fast-GBS](https://bmcbioinformatics.biomedcentral.com/counter/pdf/10.1186/s12859-016-1431-9.pdf) pipeline.

## Fast-GBS
Fast-Genotyping-by-sequencing is a bioinformatics pipline which makes use of GBS, which itself is a genotyping approach which makes use of NGS to scan a genome.
This technique allows for the simultaneous discovery and genotyping of countless SNPs across a wide range of species.
Fast-GBS attempts to to simplify the process of performing genotyping on the huge amount of data produced by NGS, while simultaneously speeding up this process.
It is capable of handling data from different sequencing platforms and can automatically detect different types of variants.

### Input
Fast-GBS takes a FASTQ file containing sequenced DNA fragments from any restriction enzyme–based GBS protocol as its main input. 
Next to this it also uses a configuration file (referred to as parameter file) which contains processing options as well as file paths.
These paths can refer to the input FASTQ file, a reference genome file, or a file containing barcodes.

### Steps
Fast-GBS uses publicly available packages for its main analysis pipeline, while using internally developed scripts to guide the process and create directories.
The following is a list of used packages and their purpose:
- Sabre -> demultiplexing
- Cutadapt -> trimming & cleaning
- BWA -> alignment
- Platypus -> post-processing of mapping, haplotype construction, and variant calling

### Output
The main output of Fast-GBS is a Variant Call Format (VCF) file containing data on all found variants.
Accompanying the VCF file are a text file containing only the genotypic data and a log file.

## Research question
How much does remaking Fast-GBS using snakemake benefit the pipeline's performance and ease of use?

## Testing
Testing of the pipline, while it is still under development, will be performed by feeding a (yet to be determined) rather small FASTQ file into the pipeline.
Using this smaller FASTQ file will ensure a rapid pace of development.  
Testing of performance of the final product will be done by reproducing one or multiple of the analysis procedures described by the original paper. 
An analysis like this may take hours to complete, so the amount of benchmarking done will depend on the available amount of time left before the deadline.

## Visualising results
Benchmarking of the pipelines performance may be visualised by producing a report containing one or multiple bar graphs comparing said performance.
Apart from this, a DAG file visualising the pipeline will be included.
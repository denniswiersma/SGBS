# SGBS
Snake Genotyping By Sequencing. A reimplementation of Fast-GBS using snakemake.

## Usage

First make a conda environment using the `environment.yaml` file:  
`conda env create -f environment.yml`

Make sure you have:
- Your paired ended data as `{sample_name}_1.fq.gz` and `{sample_name}_2.fq.gz`.
- A reference genome.
- An adapter fasta file.
- A barcodes fasta file.

Then run the pipeline using the following command:

`snakemake -c10 --use-conda`

## Divergence from Fast-GBS

- Uses [Flexbar](https://github.com/seqan/flexbar) instead of Sabre.
- Uses [Flexbar](https://github.com/seqan/flexbar) instead of Cutadapt.
- Uses [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) instead of Platypus.

## Gathering test data
The original paper, which can be found in this repository at `/docs/Fast-GBSPaper.pdf`, says the following about the data used to test the pipeline:  
"To test the performance of Fast-GBS, we used existing sequence datasets for panels of 24 unrelated accessions / clones for three species covering a range of genomic situations: soybean [22], barley [Abed et al., unpublished], and potato [Bastien et al., unpublished]."

Reference [22] links to [a paper](https://dx.doi.org/10.1371/journal.pone.0131533) which lists the data as available under study accession `SRP059747` in the NCBI SRA database. However, querying the database using this accession number returns 324 results. Way more than the 24 mentioned in the paper, so which ones to use?
Luckily this same accession number nets four results on PubMed Central, [the first of which](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6652137/) links to an [xlsx file](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6652137/bin/12859_2019_2859_MOESM2_ESM.xlsx) containing the relevant accession numbers. These are:

```
SRR2073085
SRR2073084
SRR2073083
SRR2073082
SRR2073081
SRR2073080
SRR2073079
SRR2073078
SRR2073077
SRR2073076
SRR2073075
SRR2073074
SRR2073073
SRR2073072
SRR2073071
SRR2073070
SRR2073069
SRR2073068
SRR2073067
SRR2073066
SRR2073065
SRR2073064
SRR2073063
```

This data can be downloaded using the [NCBI SRA Toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software) and using the `fasterq-dump SRR2073085` command.
Alternatively, the data can first be downloaded as `SRR2073085.sra` using the `prefetch SRR2073085` command. This will download the data to the directory set up during the SRA Toolkit installation. The data can then be converted to fastq format using the `fasterq-dump SRR2073085` command.
This data should be placed in the `resources/data/` directory, although this can be changed in the `config.yaml` file. It should be in the `fq.gz` format.

Additionally, a reference genome is needed. For this a reference genome for Glycine max (soybean) was used, which can be downloaded from [here](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000004515.6/).
This reference genome should be placed in the `resources` directory, although this can be changed in the `config.yaml` file.

Last but not least, adapters and barcodes should be provided in fasta format. The path to these files can again be set in the `config.yaml` file.

## Dependencies
- [bwa](https://github.com/lh3/bwa)
- [samtools](http://www.htslib.org/)
- [flexbar](https://github.com/seqan/flexbar)
- [gatk](https://gatk.broadinstitute.org/hc/en-us)
- [picard](https://broadinstitute.github.io/picard/)

## DAG
[DAG of pipeline](dag.png)
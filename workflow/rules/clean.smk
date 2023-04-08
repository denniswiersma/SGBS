import pathlib
import re
from Bio import SeqIO

configfile: "config/config.yaml"

def get_barcodes():
    barcode_file = config["barcodes"]
    barcodes = []
    for record in SeqIO.parse(barcode_file, "fasta"):
        barcodes.append(record.id)
    return barcodes

def get_sample_metadata(path: str) -> dict[str, list[str]]:
    """
    Get sample metadata from a directory of fastq files
    :param path: path to directory of fastq files
    :return: a dictionary of sample names and their associated fastq files
    """
    # Initialise dictionary to store sample metadata
    samples = {}
    # Initialise pathlib object
    file_path = pathlib.Path(path)
    # Check if path is a directory
    if file_path.is_dir():

        # Iterate over files in directory
        for file in file_path.iterdir():
            # Check if file is a fastq file
            if file.suffixes in [[".fastq", ".gz"], [".fq", ".gz"], [".fastq"], [".fq"]]:
                # Split file name on first period
                name_parts = file.name.split(".", 1)
                # Remove _1 or _2 from sample name
                sample_name = re.sub(r"_[12]", "", name_parts[0])
                # Add sample name and file name to dictionary
                if sample_name in samples.keys():
                    samples[sample_name].append(file.name)
                else:
                    samples[sample_name] = [file.name]
    return samples


rule demultiplex:
    """
    Demultiplex fastq files using flexbar"""
    input:
        forward_read="resources/data/{sample}_1.fq.gz",
        reverse_read="resources/data/{sample}_2.fq.gz"
    output:
        forward_demultiplexed="resources/data/demultiplexed/{sample}_barcode_{barcode}_1.fastq.gz",
        reverse_demultiplexed="resources/data/demultiplexed/{sample}_barcode_{barcode}_2.fastq.gz",
        forward_unassigned="resources/data/demultiplexed/{sample}_barcode_unassigned_{barcode}_1.fastq.gz",
        reverse_unassigned="resources/data/demultiplexed/{sample}_barcode_unassigned_{barcode}_2.fastq.gz",

    message: "Demultiplexing {input}"
    log: "logs/demultiplexing_{sample}_{barcode}.log"
    shell:
        "(flexbar -r {input.forward_read} -p {input.reverse_read} --barcodes {config[barcodes]} --target resources/data/demultiplexed/{wildcards.sample} --barcode-unassigned --zip-output GZ) >{log} 2>&1"


rule trim_adapters:
    """
    Trim adapters from fastq files"""
    input:
        forward_read="resources/data/demultiplexed/{sample}_barcode_{barcode}_1.fastq.gz",
        reverse_read="resources/data/demultiplexed/{sample}_barcode_{barcode}_2.fastq.gz"
    output:
        forward_read="resources/data/trimmed/{sample}_{barcode}_1.fastq.gz",
        reverse_read="resources/data/trimmed/{sample}_{barcode}_2.fastq.gz"
    message: "Trimming adapters from {input}"
    log: "logs/adapter_trimming_{sample}_{barcode}.log"
    shell:
        "(flexbar -r {input.forward_read} -p {input.reverse_read} -a {config[adapter]} -t resources/data/trimmed/{wildcards.sample}_{wildcards.barcode} --zip-output GZ) > {log} 2>&1"

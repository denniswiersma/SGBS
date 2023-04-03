import pathlib
import re

configfile: "config/config.yaml"

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
        forward_read = expand("resources/data/{sample[0]}", sample=get_sample_metadata(config["data"]).values()),
        reverse_read = expand("resources/data/{sample[1]}", sample=get_sample_metadata(config["data"]).values())
    output:
        samples = "resources/data/demultiplexed/{barcode}.fastq.gz"
    message: "Demultiplexing {input}"
    log: "logs/demultiplexing_{barcode}.log"
    shell:
        "(flexbar -r {input.forward_read} -p {input.reverse_read} --barcodes {config[barcodes]} --target resources/data/demultiplexed/flexbarOut --barcode-unassigned --zip-output GZ) >{log} 2>&1"

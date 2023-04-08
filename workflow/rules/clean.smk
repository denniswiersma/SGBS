configfile: "config/config.yaml"

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
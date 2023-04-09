configfile: "config/config.yaml"

rule demultiplex:
    """
    Demultiplex fastq files using flexbar
    """
    input:
        forward_read="resources/data/{sample}_1.fq.gz",
        reverse_read="resources/data/{sample}_2.fq.gz"
    output:
        forward_demultiplexed="resources/data/demultiplexed/{sample}_barcode_{barcode}_1.fastq.gz",
        reverse_demultiplexed="resources/data/demultiplexed/{sample}_barcode_{barcode}_2.fastq.gz",

    message: "Demultiplexing {input}"
    log: "logs/flexbar/demultiplexing_{sample}_{barcode}.log"
    threads: config["flexbar_demultiplex_threads"]
    shell:
        "(flexbar -n {threads} -r {input.forward_read} -p {input.reverse_read} --barcodes {config[barcodes]} "
        "--target resources/data/demultiplexed/{wildcards.sample} --zip-output GZ) >{log} 2>&1"


rule trim_adapters:
    """
    Trim adapters from fastq files
    """
    input:
        forward_read="resources/data/demultiplexed/{sample}_barcode_{barcode}_1.fastq.gz",
        reverse_read="resources/data/demultiplexed/{sample}_barcode_{barcode}_2.fastq.gz"
    output:
        forward_read="resources/data/trimmed/{sample}_{barcode}_1.fastq.gz",
        reverse_read="resources/data/trimmed/{sample}_{barcode}_2.fastq.gz"
    message: "Trimming adapters from {input}"
    log: "logs/flexbar/adapter_trimming_{sample}_{barcode}.log"
    threads: config["flexbar_trim_threads"]
    shell:
        "(flexbar -n {threads} -r {input.forward_read} -p {input.reverse_read} -a {config[adapter]} "
        "-t resources/data/trimmed/{wildcards.sample}_{wildcards.barcode} --zip-output GZ) > {log} 2>&1"
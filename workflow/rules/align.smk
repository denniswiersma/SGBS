rule align:
    """
    Aligns trimmed reads to reference genome using BWA-MEM
    """
    input:
        forward_read = "resources/data/trimmed/{sample}_{barcode}_1.fastq.gz",
        reverse_read = "resources/data/trimmed/{sample}_{barcode}_2.fastq.gz",
        bwa_amb = config["reference"] + ".amb",
        bwa_ann = config["reference"] + ".ann",
        bwa_bwt = config["reference"] + ".bwt",
        bwa_pac = config["reference"] + ".pac",
        bwa_sa = config["reference"] + ".sa",
        sam_fai = config["reference"] + ".fai"
    output:
        "resources/data/aligned/{sample}_{barcode}.sam"
    log:
        "logs/align_{sample}_{barcode}.log"
    message: "Aligning {input.forward_read} and {input.reverse_read} to reference genome"
    shell:
        "bwa mem {config[reference]} {input.forward_read} {input.reverse_read} > {output} 2> {log}"
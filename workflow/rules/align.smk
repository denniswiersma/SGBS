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

rule post_processing_view:
    input:
        "resources/data/aligned/{sample}_{barcode}.sam"
    output:
        "resources/data/processed/{sample}_{barcode}.bam"
    log:
        "logs/post_processing_view_{sample}_{barcode}.log"
    message: "Post-processing {input} with samtools view"
    shell:
        "(samtools view -bSh {input} > {output}) >{log} 2>&1"

rule post_processing_sort:
    input:
        "resources/data/processed/{sample}_{barcode}.bam"
    output:
        "resources/data/processed/{sample}_{barcode}_sorted.bam"
    log:
        "logs/post_processing_sort_{sample}_{barcode}.log"
    message: "Post-processing {input} with samtools sort"
    shell:
        "(samtools sort {input} -o {output}) >{log} 2>&1"

rule post_processing_index:
    input:
        "resources/data/processed/{sample}_{barcode}_sorted.bam"
    output:
        "resources/data/processed/{sample}_{barcode}_sorted.bam.bai"
    log:
        "logs/post_processing_index_{sample}_{barcode}.log"
    message: "Post-processing {input} with samtools index"
    shell:
        "(samtools index {input}) >{log} 2>&1"
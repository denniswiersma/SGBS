rule bwa_index:
    input:
        reference = config["reference"]
    output:
        config["reference"] + ".amb",
        config["reference"] + ".ann",
        config["reference"] + ".bwt",
        config["reference"] + ".pac",
        config["reference"] + ".sa"
    message: "Indexing reference genome with BWA"
    shell:
        "bwa index -a bwtsw {input.reference}"

rule samtools_index:
    input:
        reference = config["reference"]
    output:
        config["reference"] + ".fai"
    message: "Indexing reference genome with samtools"
    shell:
        "samtools faidx {input.reference}"
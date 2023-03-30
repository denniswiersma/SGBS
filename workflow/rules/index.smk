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
    log:
        "logs/bwa_index.log"
    shell:
        "(bwa index -a bwtsw {input.reference}) >{log} 2>&1"

rule samtools_index:
    input:
        reference = config["reference"]
    output:
        config["reference"] + ".fai"
    message: "Indexing reference genome with samtools"
    log:
        "logs/samtools_index.log"
    shell:
        "(samtools faidx {input.reference}) >{log} 2>&1"
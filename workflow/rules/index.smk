configfile: "config/config.yaml"

rule bwa_index:
    """
    Index reference genome with BWA.
    """
    input:
        reference = config["reference"]
    output:
        config["reference"] + ".amb",
        config["reference"] + ".ann",
        config["reference"] + ".bwt",
        config["reference"] + ".pac",
        config["reference"] + ".sa"
    message: "Indexing reference genome with BWA. This may take a few minutes."
    log:
        "logs/bwa/bwa_index.log"
    shell:
        "(bwa index -a bwtsw {input.reference}) >{log} 2>&1"

rule samtools_index:
    """
    Index reference genome with samtools.
    """
    input:
        reference = config["reference"]
    output:
        config["reference"] + ".fai"
    message: "Indexing reference genome with samtools."
    log:
        "logs/samtools/samtools_index.log"
    shell:
        "(samtools faidx {input.reference}) >{log} 2>&1"
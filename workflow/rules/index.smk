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
    message: "Indexing reference {input} with BWA."
    conda: "SGBS"
    log:
        "logs/bwa/bwa_index.log"
    params:
        batch_size = config["bwa_index_batch_size"]
    shell:
        "(bwa index -a bwtsw {input.reference} -b {params.batch_size}) >{log} 2>&1"

rule samtools_index:
    """
    Index reference genome with samtools.
    """
    input:
        reference = config["reference"]
    output:
        config["reference"] + ".fai"
    message: "Indexing reference {input} with samtools."
    conda: "SGBS"
    log:
        "logs/samtools/samtools_index.log"
    shell:
        "(samtools faidx {input.reference}) >{log} 2>&1"
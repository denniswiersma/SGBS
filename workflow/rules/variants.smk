ref_name = config["reference"].split(".")[0]

rule create_sequence_dict:
    input:
        reference = config["reference"]
    output:
        ref_name + ".dict"
    log:
        "logs/gatk/create_sequence_dict.log"
    shell:
        "(picard CreateSequenceDictionary R={input.reference} O={output}) >{log} 2>&1"

rule variant_calling:
    """
    HaplotypeCaller is a tool that performs local re-assembly of haplotypes in an active region to discover SNPs and indels.
    """
    input:
        bam = "resources/data/processed/{sample}_{barcode}_sorted.bam",
        bai = "resources/data/processed/{sample}_{barcode}_sorted.bam.bai",
        reference = config["reference"],
        dict = ref_name + ".dict"
    output:
        "resources/data/variants/{sample}_{barcode}_variants.g.vcf"
    log:
        "logs/gatk/{sample}_{barcode}_haplotypecaller.log"
    message: "Calling variants for {input.bam}"
    conda: "SGBS"
    params:
        min_base_quality = config["gatk_min_baseq"],
    shell:
        "(gatk HaplotypeCaller -R {input.reference} -I {input.bam} -O {output} "
        "--min-base-quality-score {params.min_base_quality}) >{log} 2>&1"

rule post_variant_calling:
    """
    GenotypeGVCFs is a tool that performs joint genotyping on one or more samples pre-called with HaplotypeCaller in GVCF mode.
    """
    input:
        gvcf = "resources/data/variants/{sample}_{barcode}_variants.g.vcf",
        reference = config["reference"]
    output:
        "resources/data/variants/{sample}_{barcode}_variants.vcf"
    log:
        "logs/gatk/{sample}_{barcode}_genotypegvcfs.log"
    message: "Post-processing variants for {input.gvcf}"
    conda: "SGBS"
    shell:
        "(gatk GenotypeGVCFs -R {input.reference} -V {input.gvcf} -O {output}) >{log} 2>&1"
# MODULE 4: 

rule kat_com:
    input:
        fq_1 = config["results"]["trimmed"] + "/{Sample_ID}_P1.fastq.gz",
        fq_2 = config["results"]["trimmed"] + "/{Sample_ID}_P2.fastq.gz",
        ref  = config["ref_genome"]
    output:
        comp = config["results"]["kat"] + "/{Sample_ID}/{Sample_ID}_vs_ref.mx.spectra-cn.png"
    conda:
        "../envs/kat.yaml"
    shell:
        """
        mkdir -p {config[results][kat]}/{wildcards.Sample_ID}
        kat comp -o {config[results][kat]}/{wildcards.Sample_ID}/{wildcards.Sample_ID}_vs_ref {input.fq_1} {input.ref}  > {config[results][kat]}/{wildcards.Sample_ID}/kat_log.txt 2>&1
        """

rule kat_hist:
    input:
        fq_1 = config["results"]["trimmed"] + "/{Sample_ID}_P1.fastq.gz",
        fq_2 = config["results"]["trimmed"] + "/{Sample_ID}_P2.fastq.gz"
    output:
        hist = config["results"]["kat"] + "/{Sample_ID}/{Sample_ID}_hist.png"
    conda:
        "../envs/kat.yaml"
    shell:
        """
        mkdir -p {config[results][kat]}/{wildcards.Sample_ID}
        kat hist -o {config[results][kat]}/{wildcards.Sample_ID}/{wildcards.Sample_ID}_hist {input.fq_1} {input.fq_2}  > {config[results][kat]}/{wildcards.Sample_ID}/kat_log.txt 2>&1
        """

rule kat_gc:
    input:
        fq_1 = config["results"]["trimmed"] + "/{Sample_ID}_P1.fastq.gz",
        fq_2 = config["results"]["trimmed"] + "/{Sample_ID}_P2.fastq.gz"
    output:
        gc = config["results"]["kat"] + "/{Sample_ID}/{Sample_ID}_gc.png"
    conda:
        "../envs/kat.yaml"
    shell:
        """
        mkdir -p {config[results][kat]}/{wildcards.Sample_ID}
        kat gcp -o {config[results][kat]}/{wildcards.Sample_ID}/{wildcards.Sample_ID}_gc {input.fq_1} {input.fq_2}  > {config[results][kat]}/{wildcards.Sample_ID}/kat_log.txt 2>&1
        """
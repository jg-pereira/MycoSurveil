# MODULE 1: Renaming, QC, Trimming, 
#########################################
# workflow/rules/1_rename_files
#########################################

rule rename_files:
    wildcard_constraints:
        Sample_ID="|".join(SAMPLES)

    input:
        r1=lambda wc: config["input_dir"] + "/" + short_inverse[wc.Sample_ID][0] + short_inverse[wc.Sample_ID][1] + short_inverse[wc.Sample_ID][3],
        r2=lambda wc: config["input_dir"] + "/" + short_inverse[wc.Sample_ID][0] + short_inverse[wc.Sample_ID][2] + short_inverse[wc.Sample_ID][3]

    output:
        r1=config["input_dir"] + "/{Sample_ID}_1.fastq.gz",
        r2=config["input_dir"] + "/{Sample_ID}_2.fastq.gz"

    conda:
        "../envs/project.yaml"
        
    shell:
        """
        mv {input.r1} {output.r1}
        mv {input.r2} {output.r2}
        """

#########################################
# workflow/rules/2_fastqc
#########################################

rule fastqc_raw:
    input:
        r1 = config["input_dir"] + "/{Sample_ID}_1.fastq.gz",
        r2 = config["input_dir"] + "/{Sample_ID}_2.fastq.gz",
    output:
        zip = config["results"]["fastqc_raw"] + "/{Sample_ID}_1_fastqc.zip",
        html = config["results"]["fastqc_raw"] + "/{Sample_ID}_1_fastqc.html"
    params:
        outdir_qc = config["results"]["fastqc_raw"] 
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        mkdir -p {params.outdir_qc}
        fastqc {input.r1} {input.r2} -o {params.outdir_qc}
        """

rule unzip_fastqc:
    input:
        zip = config["results"]["fastqc_raw"] + "/{Sample_ID}_1_fastqc.zip",
    output:
        config["results"]["fastqc_raw"] + "/{Sample_ID}_1_fastqc/fastqc_data.txt"
    params:
        outdir_qc = config["results"]["fastqc_raw"] 
    conda:
        "../envs/unzip.yaml"
    shell:
        """
        unzip -d {params.outdir_qc} {input.zip}
        """

#########################################
# workflow/rules/4_trimmomatic
#########################################

rule trimmomatic:
    input:
        r1 = config["input_dir"] + "/{Sample_ID}_1.fastq.gz",
        r2 = config["input_dir"] + "/{Sample_ID}_2.fastq.gz",
    output:
        r1_paired = temp(config["results"]["trimmed"] + "/{Sample_ID}_P1.fastq.gz"),
        r2_paired = temp(config["results"]["trimmed"] + "/{Sample_ID}_P2.fastq.gz"), 

        r1_unpaired = temp(config["results"]["trimmed"] + "/{Sample_ID}_U1.fastq.gz"),
        r2_unpaired = temp(config["results"]["trimmed"] + "/{Sample_ID}_U2.fastq.gz")
    params:
        adapters = "TruSeq3-PE.fa",
        threads = 4,
        outdir_trimmed = config["results"]["trimmed"]+ "/{Sample_ID}"
    conda:
        "../envs/trimmomatic.yaml"
    shell:
        """
        echo "DEBUG: trimmed = {params.outdir_trimmed}"
        mkdir -p {params.outdir_trimmed}
        trimmomatic PE \
            -threads {params.threads} \
            {input.r1} {input.r2} \
            {output.r1_paired} {output.r1_unpaired} \
            {output.r2_paired} {output.r2_unpaired} \
            ILLUMINACLIP:{params.adapters}:2:30:10 \
            LEADING:10 TRAILING:10 \
            SLIDINGWINDOW:4:15 MINLEN:31
        """

#########################################
# workflow/rules/5_fastqc_trimmed
#########################################

rule fastqc_trimmed:
    input:
        r1_paired = config["results"]["trimmed"] + "/{Sample_ID}_P1.fastq.gz",
        r2_paired = config["results"]["trimmed"] + "/{Sample_ID}_P2.fastq.gz"
    output:
        zip = config["results"]["fastqc_trimmed"] + "/{Sample_ID}_P1_fastqc.zip",
        html = config["results"]["fastqc_trimmed"] + "/{Sample_ID}_P1_fastqc.html"
    params:
        outdir_qc = config["results"]["fastqc_trimmed"] 
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        mkdir -p {params.outdir_qc}
        fastqc {input.r1_paired} {input.r2_paired} -o {params.outdir_qc}
        """

rule unzip_fastqc_trimmed:
    input:
        zip = config["results"]["fastqc_trimmed"] + "/{Sample_ID}_P1_fastqc.zip",
    output:
        config["results"]["fastqc_trimmed"] + "/{Sample_ID}_P1_fastqc/fastqc_data.txt"
    conda:
        "../envs/unzip.yaml"
    shell:
        """
        unzip -d {config[results][fastqc_trimmed]} {input.zip}
        """
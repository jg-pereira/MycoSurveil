# MODULE 2: 

#########################################
# workflow/rules/1_index_ref_genome
#########################################
rule bwa_index:
    input:
        REF
    output:
        REF + ".bwt",
        REF + ".pac",
        REF + ".ann",
        REF + ".amb",
        REF + ".sa"
    conda:
        "../envs/mapping.yaml"
    shell:
        "bwa index {input}"

rule faidx:
    input:
        REF
    output:
        REF + ".fai"
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools faidx {input}"

rule gatk_dict:
    input:
        REF
    output:
        REF_NO_EXT + ".dict"
        #"ref_genome/genome_albicans.dict"
    conda:
        "../envs/mapping.yaml"
    shell:
        "gatk CreateSequenceDictionary -R {input} -O {output}"

#########################################
# workflow/rules/2_MAPPING WITH BWA MEM
#########################################

rule bwa_mem:
    input:
        REF = REF,
        R1 = config["results"]["trimmed"] + "/{Sample_ID}_P1.fastq.gz", 
        R2 = config["results"]["trimmed"] + "/{Sample_ID}_P2.fastq.gz"
    output:
        temp(config["results"]["bams"] + "/{Sample_ID}/{Sample_ID}.sam")
    conda:
        "../envs/mapping.yaml"
    params:
        outdir = config["results"]["bams"] 
    threads: 8
    shell:
        """
        mkdir -p {params.outdir}/{wildcards.Sample_ID}
        bwa mem -t {threads} \
        -R "@RG\\tID:{wildcards.Sample_ID}\\tPL:ILLUMINA\\tSM:{wildcards.Sample_ID}" \
        {input.REF} {input.R1} {input.R2} > {output}
        """

# #########################################
# workflow/rules/3_SORTED BY COORDINATES
# #########################################

rule sort_bam:
    input:
        config["results"]["bams"] + "/{Sample_ID}/{Sample_ID}.sam"
    output:
        temp(config["results"]["bams"] + "/{Sample_ID}/{Sample_ID}.sorted.bam")
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        gatk SortSam \
            -I {input} \
            -O {output} \
            -SO coordinate
        """

# #########################################
# workflow/rules/4_MARK DUPLICATES
# #########################################

rule mark_duplicates:
    input:
        config["results"]["bams"] + "/{Sample_ID}/{Sample_ID}.sorted.bam"
    output:
        bam = config["results"]["bams"] + "/{Sample_ID}/{Sample_ID}.sorted.markdup.bam",
        metrics = config["results"]["bams"] + "/{Sample_ID}/{Sample_ID}.dup_metrics.txt"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        gatk MarkDuplicates \
            -I {input} \
            -O {output.bam} \
            -M {output.metrics}
        """

# #########################################
# workflow/rules/5_BAM INDEX
# #########################################

rule index_bam:
    input:
        bam = config["results"]["bams"] + "/{Sample_ID}/{Sample_ID}.sorted.markdup.bam"
    output:
        bai = config["results"]["bams"] + "/{Sample_ID}/{Sample_ID}.sorted.markdup.bam.bai"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        gatk BuildBamIndex \
            -I {input.bam} \
            -O {output.bai}
        """

# #########################################
# workflow/rules/6_STATISTICS AND COVERAGE
# #########################################

rule coverage_stats:
    input:
        bam = config["results"]["bams"] + "/{Sample_ID}/{Sample_ID}.sorted.markdup.bam",
        bai = config["results"]["bams"] + "/{Sample_ID}/{Sample_ID}.sorted.markdup.bam.bai" 
    output:
        coverage = config["results"]["stats"] + "/{Sample_ID}/{Sample_ID}.coverage.txt",
        flagstat = config["results"]["stats"] + "/{Sample_ID}/{Sample_ID}.mapping_flagstats.txt",
        stats    = config["results"]["stats"] + "/{Sample_ID}/{Sample_ID}.mapping_stats.txt"
    conda:
        "../envs/samtools_env.yaml"
    params:
        outdir = config["results"]["stats"]
    shell:
        """
        mkdir -p {params.outdir}/{wildcards.Sample_ID}
        
        samtools flagstat {input.bam} > {output.flagstat}
        samtools stats {input.bam} > {output.stats}
        samtools coverage {input.bam} > {output.coverage}
        """

#########################################
# workflow/rules/7_Ploidy ESTIMATION
#########################################

rule nquire_create:
    input:
        bam = config["results"]["bams"] + "/{Sample_ID}/{Sample_ID}.sorted.markdup.bam"
    output:
        bin = config["results"]["nquire"] + "/{Sample_ID}.bin"
    params:
        exe = config["paths"]["nquire_bin"],
        out_dir = config["results"]["nquire"]
    shell:
        """
        mkdir -p {params.out_dir}
        
        {params.exe} create -b {input.bam} -o {params.out_dir}/{wildcards.Sample_ID} -q 20
        """

rule nquire_denoise:
    input:
        bin = config["results"]["nquire"] + "/{Sample_ID}.bin"
    output:
        bin_clean = config["results"]["nquire"] + "/{Sample_ID}_denoised.bin"
    params:
        exe = config["paths"]["nquire_bin"],
        out_dir = config["results"]["nquire"]
    shell:
        """
        {params.exe} denoise {input.bin} -o {params.out_dir}/{wildcards.Sample_ID}_denoised
        """

rule nquire_lrdmodel:
    input:
        bin = config["results"]["nquire"] + "/{Sample_ID}_denoised.bin"
    output:
        txt = config["results"]["nquire"] + "/{Sample_ID}_lrdmodel.txt"
    params:
        exe = config["paths"]["nquire_bin"]
    shell:
        """
        {params.exe} lrdmodel {input.bin} > {output.txt}
        """

rule nquire_histotest:
    input:
        bin = config["results"]["nquire"] + "/{Sample_ID}_denoised.bin"
    output:
        txt = config["results"]["nquire"] + "/{Sample_ID}_histotest.txt"
    params:
        exe = config["paths"]["nquire_bin"]
    shell:
        """
        {params.exe} histotest {input.bin} > {output.txt}
        """
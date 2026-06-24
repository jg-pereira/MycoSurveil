# MODULE 3: 

VCF_DIR = config["results"]["vcf"]
BAM_DIR = config["results"]["bams"]
REF_GENOME = config["ref_genome"] 
PYLO  = config["results"]["phylo"]


########################################
#workflow/rules/1_Calling_variants
########################################
rule haplotypecaller:
    input:
        bam = BAM_DIR + "/{Sample_ID}/{Sample_ID}.sorted.markdup.bam",
        bai = BAM_DIR + "/{Sample_ID}/{Sample_ID}.sorted.markdup.bam.bai",
        ref = REF_GENOME
    output:
        vcf = temp(VCF_DIR + "/{Sample_ID}/{Sample_ID}.raw.vcf"),
    conda:
        "../envs/mapping.yaml"
    params:
        outdir = VCF_DIR,
    shell:
        """
        mkdir -p {params.outdir}/{wildcards.Sample_ID}        
        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.vcf} \
            --standard-min-confidence-threshold-for-calling 30 \
            --ploidy 2 \
        """

#########################################
#workflow/rules/2_Variant_Filtration SNP
#########################################
rule variantfiltration:
    input:
        vcf = VCF_DIR + "/{Sample_ID}/{Sample_ID}.raw.vcf",
        ref = REF_GENOME
    output:
        vcf = VCF_DIR + "/{Sample_ID}/{Sample_ID}.filtered.vcf",
    conda:
        "../envs/mapping.yaml"
    shell:
        """       
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.vcf} \
            -O {output.vcf} \
            -G-filter-name "heterozygous" \
            -G-filter "isHet == 1" \
            --filter-name "BadDepthofQualityFilter" \
            --filter "DP <= 20 || QD < 2.0 || MQ < 40.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
            --cluster-size 5 \
            --cluster-window-size 20
        """

# #########################################
# workflow/rules/3_Select SNPs and INDELs
# #########################################

rule select_snps:
    input:
        vcf = VCF_DIR + "/{Sample_ID}/{Sample_ID}.filtered.vcf",
        ref = REF_GENOME
    output:
        vcf = VCF_DIR + "/{Sample_ID}/{Sample_ID}.snps.vcf"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            --select-type-to-include SNP \
            -O {output.vcf}
        """

rule select_indels:
    input:
        vcf = VCF_DIR + "/{Sample_ID}/{Sample_ID}.filtered.vcf",
        ref = REF_GENOME
    output:
        vcf = VCF_DIR + "/{Sample_ID}/{Sample_ID}.indels.vcf"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            --select-type-to-include INDEL \
            -O {output.vcf}
        """

rule compress_vcfs:
    input: 
        vcf_snps = config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.snps.vcf",
        vcf_indels = config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.indels.vcf",
        vcf_all = config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.filtered.vcf"
    output:
        vcf_snps_gz = config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.snps.vcf.gz",
        vcf_indels_gz = config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.indels.vcf.gz",
        vcf_all_gz = config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.filtered.vcf.gz"
    conda:
        "../envs/tabix.yaml"
    shell:
        """
        bgzip -c {input.vcf_snps} > {output.vcf_snps_gz}
        bgzip -c {input.vcf_indels} > {output.vcf_indels_gz}
        bgzip -c {input.vcf_all} > {output.vcf_all_gz}
        """

# #########################################
# workflow/rules/4_Variant_statistics_metadata
# #########################################

rule  nr_pass_hom_het:
    input: 
        vcf = config["results"] ["vcf"]+ "/{Sample_ID}/{Sample_ID}.snps.vcf.gz"
    output:
        counts_homo = config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.pass_hom.txt",
        counts_het = config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.pass_het.txt",
    conda: 
        "../envs/bcf.yaml"
    shell:
        """
        bcftools view -v snps -f PASS -g hom {input.vcf} -H | wc -l > {output.counts_homo}
        bcftools view -v snps -f PASS -g het {input.vcf} -H | wc -l > {output.counts_het}
        """

rule count_variants:
    input:
        vcf = config["results"] ["vcf"]+ "/{Sample_ID}/{Sample_ID}.filtered.vcf"
    output:
        counts = config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.variant_counts.txt",
        stats = config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.variant_stats.txt"
    conda:
        "../envs/bcf.yaml"
    shell:
        """
        bcftools view -H {input.vcf} | wc -l > {output.counts}
        bcftools stats {input.vcf} > {output.stats}
        """

rule total_snps:
    input:
        vcf_snp = config["results"] ["vcf"]+ "/{Sample_ID}/{Sample_ID}.snps.vcf.gz"
    output:
        counts_snps = config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.total_snps.txt"
    conda:
        "../envs/bcf.yaml"
    shell:
        """
        bcftools view -v snps {input.vcf_snp} -H | wc -l > {output.counts_snps}
        """

rule total_indels:
    input:
        vcf_indel = config["results"] ["vcf"]+ "/{Sample_ID}/{Sample_ID}.indels.vcf.gz"
    output:
        counts_indels = config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.total_indels.txt"
    conda:
        "../envs/bcf.yaml"
    shell:
        """
        bcftools view -v indels {input.vcf_indel} -H | wc -l > {output.counts_indels}
        """

rule update_metadata:
    input:
        metadata = config["metadata"],
        mapping_flagstats = expand(config["results"]["stats"] + "/{Sample_ID}/{Sample_ID}.mapping_flagstats.txt", Sample_ID=SAMPLES),
        coverage = expand(config["results"]["stats"] + "/{Sample_ID}/{Sample_ID}.coverage.txt", Sample_ID=SAMPLES),
        ploidy = expand(config["results"]["nquire"] + "/{Sample_ID}_histotest.txt", Sample_ID=SAMPLES),
        fastqc_file_raw = expand(config["results"]["fastqc_raw"] + "/{Sample_ID}_1_fastqc/fastqc_data.txt", Sample_ID=SAMPLES),
        fastqc_file_trimmed = expand(config["results"]["fastqc_trimmed"] + "/{Sample_ID}_P1_fastqc/fastqc_data.txt", Sample_ID=SAMPLES),
        nr_variants = expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.variant_counts.txt", Sample_ID=SAMPLES),
        sample_snps = expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.snps.vcf", Sample_ID=SAMPLES),
        sample_indels = expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.indels.vcf", Sample_ID=SAMPLES),
        pass_homo = expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.pass_hom.txt", Sample_ID=SAMPLES),
        pass_hetero = expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.pass_het.txt", Sample_ID=SAMPLES),
    output:
        updated_metadata = config["results"]["base"] + "/updated_metadata.xlsx"
    params:
        outdir = config["results"]["base"]
    shell:
        """
        python workflow/scripts/update_metadata.py \
            --metadata {config[metadata]} \
            --mapping_flagstats {input.mapping_flagstats} \
            --coverage {input.coverage} \
            --ploidy {input.ploidy} \
            --fastqc_raw {input.fastqc_file_raw} \
            --fastqc_trimmed {input.fastqc_file_trimmed} \
            --nr_variants {input.nr_variants} \
            --snps {input.sample_snps} \
            --indels {input.sample_indels} \
            --output_metadata {params.outdir} \
            --pass_homo {input.pass_homo} \
            --pass_hetero {input.pass_hetero}
        """

#**********************************************
#         PHYLOGENETIC ANALYSIS 
#**********************************************

###############################################
# workflow/rules/5_VCFs TO SNPs only PASS + HOM
###############################################

rule vcf_pass_hom:
    input:
        vcf = VCF_DIR + "/{Sample_ID}/{Sample_ID}.snps.vcf.gz"
    output:
        vcf_pass_hom = VCF_DIR + "/{Sample_ID}/{Sample_ID}.snps.pass.hom.vcf.gz"
    conda:
        "../envs/bcf.yaml"
    shell:
        """
        bcftools view -v snps -f PASS -g hom {input.vcf} -Oz -o {output.vcf_pass_hom}
        """

rule vcf_pass_hom_index:
    input:
        index_samples = config["results"] ["vcf"] + "/{Sample_ID}/{Sample_ID}.snps.pass.hom.vcf.gz"
    output:
        tbi = config["results"] ["vcf"] + "/{Sample_ID}/{Sample_ID}.snps.pass.hom.vcf.gz.tbi"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        gatk IndexFeatureFile -I {input.index_samples}
        """

rule vcf_pass_hom_TO_fasta:
    input:
        ref = REF_GENOME,
        vcf = config["results"] ["vcf"] + "/{Sample_ID}/{Sample_ID}.snps.pass.hom.vcf.gz",
        tbi = config["results"] ["vcf"] + "/{Sample_ID}/{Sample_ID}.snps.pass.hom.vcf.gz.tbi"
    output:
        fasta = temp(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.genome.fasta")
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        gatk FastaAlternateReferenceMaker \
            -R {input.ref} \
            -V {input.vcf} \
            -O {output.fasta}
        """

rule rename_chromossomes_in_fasta:
    input:
        fasta = config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.genome.fasta"
    output:
        fasta_renamed = temp(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.genome.renamed.fasta")
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        sed -E 's/^>[^ ]+ ([^:]+):.*/>\\1/' {input.fasta} > {output.fasta_renamed}
        """

##########################################
# workflow/rules/6_BED of reference genome
##########################################

rule bed_genome_reference:
    input:
        index_ref = REF + ".fai"
    output:
        bed_ref = config["results"]["vcf"] + "/sites.reference.genome.bed"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        awk '{{print $1"\t0\t"$2}}' {input.index_ref} > {output.bed_ref}
        """

#######################################################
# workflow/rules/7_Exclude all SNPs+PASS+hom from VCFs
#######################################################

rule vcf_exclude_pass:
    input:
        vcf_all = config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.filtered.vcf.gz",  #I=10920
        vcf_pass_hom = config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.snps.pass.hom.vcf.gz"  #I=37
    output:
        not_pass_vcf = config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.ind.snps.bad.vcf.gz"  #I=10883
    conda:  
        "../envs/bcf.yaml"
    shell:
        """
        bcftools index {input.vcf_all} 
        bcftools index {input.vcf_pass_hom}
        bcftools isec -C -w1 {input.vcf_all} {input.vcf_pass_hom} -Oz -o {output.not_pass_vcf}
        """

#############################################################################################
# workflow/rules/8_Posições do VCF que não são SNPs+PASS+hom convertidas para BED (bad sites)
#############################################################################################

rule vcf_to_bed:
    input:
        not_pass_vcf = config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.ind.snps.bad.vcf.gz"
    output:
        bad_sites_bed = config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.ind.snps.bad.bed"
    conda:
        "../envs/bcf.yaml"
    shell:
        """
        bcftools query -f '%CHROM\t%POS0\t%POS\n' {input.not_pass_vcf} > {output.bad_sites_bed}
        """

###############################################################################################################################################
# workflow/rules/9_Global merge all bad sites, and subtract to reference genome to get good sites for phylogenetic analysis
###############################################################################################################################################

rule global_good_site:
    input:
        bad_sites = expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.ind.snps.bad.bed", Sample_ID=SAMPLES),
        bed_ref = config["results"]["vcf"] + "/sites.reference.genome.bed"
    output:
        good_sites_bed = config["results"]["vcf"] + "/good_sites.global.bed",
        all_bad_sites_tmp = config["results"]["vcf"] + "/all_bad_sites.bed"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        cat {input.bad_sites} | sort -k1,1 -k2,2n | uniq > {output.all_bad_sites_tmp}
        bedtools subtract -a {input.bed_ref} -b {output.all_bad_sites_tmp} > {output.good_sites_bed}
        """

########################################################################
# workflow/rules/10_good sites in fasta format for phylogenetic analysis
########################################################################

rule good_sites:
    input:
        fasta = config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.genome.renamed.fasta",
        good_sites_bed = config["results"]["vcf"] + "/good_sites.global.bed",
    output:
        fasta_clean = temp(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.genome.clean.fasta"),
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools getfasta -fi {input.fasta} -bed {input.good_sites_bed} -fo {output.fasta_clean}
        """

rule good_sites_ref:
    input:
        good_sites_bed = config["results"]["vcf"] + "/good_sites.global.bed",
        ref_seq = REF_GENOME
    output:
        Refseq = config["results"]["vcf"] + "/ref.seq.fasta"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools getfasta -fi {input.ref_seq} -bed {input.good_sites_bed} -fo {output.Refseq}
        """
########################################################################################
# worflow/rules/11_concatenate all reagions of the fasta to get the full genome sequence 
########################################################################################

rule make_genome_sample_full:
    input:
        clean_fasta = config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.genome.clean.fasta"
    output:
        full_fasta = temp(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.genome_full.fasta")
    shell:
        """
        echo ">{wildcards.Sample_ID}" > {output.full_fasta}
        grep -v ">" {input.clean_fasta} | tr -d '\\n' >> {output.full_fasta}
        echo "" >> {output.full_fasta}
        """

rule make_genome_ref_full:
    input:
        Refseq = config["results"]["vcf"] + "/ref.seq.fasta"
    output:
        full_Refseq = temp(config["results"]["vcf"] + "/ref.seq.full.fasta")
    shell:
        """
        echo ">Reference_genome" > {output.full_Refseq}
        grep -v ">" {input.Refseq} | tr -d '\\n' >> {output.full_Refseq}
        echo "" >> {output.full_Refseq}
        """

###############################################################################################
# workflow/rules/12_concatenate all samples to get the multiple genome sequence for all samples 
###############################################################################################

rule join_all_samples:
    input:
        full_fastas = expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.genome_full.fasta", Sample_ID=SAMPLES),
        full_Refseq = config["results"]["vcf"] + "/ref.seq.full.fasta"
    output:
        all_full_fasta = config["results"]["vcf"] + "/Allsamples.genome_full.fasta"
    shell:
        """
        cat {input.full_fastas} {input.full_Refseq} > {output.all_full_fasta}
        """
######################################################
# workflow/rules/13_Phylogenetic analysis with IQ-TREE
######################################################

rule iqtree:
    input: 
        fasta = config["results"]["vcf"] + "/Allsamples.genome_full.fasta"
    output:
        tree   = config["results"]["phylo"] + "/Allsamples_genome.contree",
        log    = config["results"]["phylo"] + "/Allsamples_genome.log",
        iqtree = config["results"]["phylo"] + "/Allsamples_genome.iqtree"
    conda:
        "../envs/iqtree.yaml"
    params:
        outdir = config["results"]["phylo"]
    shell:
        """
        mkdir -p {params.outdir}
        iqtree -s {input.fasta} -m GTR -B 1000 -alrt 1000 -pre {params.outdir}/Allsamples_genome
        """

rule chanege_iqtree_to_nwk:
    input:
        tree = config["results"]["phylo"] + "/Allsamples_genome.contree"
    output:
        nwk = config["results"]["phylo"] + "/Allsamples_genome.nwk"
    shell:
        """
        cp {input.tree} {output.nwk}
        """
######################################################
# SNP-DIST
######################################################

rule distance_genetic_matrix:
    input:
        fasta = config["results"]["vcf"] + "/Allsamples.genome_full.fasta"
    output:
        dist_matrix = config["results"]["snp_dists"] + "/Allsamples_genome_snp_dist.tsv"
    conda:
        "../envs/snp-dists.yaml"
    params:
        outdir = config["results"]["snp_dists"]
    shell:
        """
        mkdir -p {params.outdir}
        snp-dists {input.fasta} > {output.dist_matrix} 
        """


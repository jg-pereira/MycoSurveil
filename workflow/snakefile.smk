#########################################
# workflow/Snakefile
#########################################
configfile: "/mnt/c/Users/joana.gomes/MycoSurveil/configs/config.yaml"
import os
from scripts.utils_samples import load_metadata, make_samples_list, mapping_samples_list

# --- 1.SAMPLES LIST ---

mapping_metadata = load_metadata(config["metadata"])
samples_list_NGS = make_samples_list(config["input_dir"])
final_mapping, short_inverse = mapping_samples_list(samples_list_NGS, mapping_metadata)
SAMPLES = list(short_inverse.keys())

print("Samples:", SAMPLES)
SAMPLES = [str(s) for s in SAMPLES]

# --- 2. REFERENCE SETTING---

REF = config["ref_genome"]
REF_NO_EXT = os.path.splitext(REF)[0]

include: "rules/1_module.smk" 
include: "rules/2_module.smk"
include: "rules/3_module.smk"
#include: "rules/4_module.smk"  


rule all:
    input:
        # --- 1. QUALITY CONTROL ---
        # Raw 
        expand(config["input_dir"] + "/{Sample_ID}_1.fastq.gz", Sample_ID=SAMPLES),
        expand(config["input_dir"] + "/{Sample_ID}_2.fastq.gz", Sample_ID=SAMPLES),
        # FastQC raw
        expand(config["results"]["fastqc_raw"] + "/{Sample_ID}_1_fastqc.zip", Sample_ID=SAMPLES),
        expand(config["results"]["fastqc_raw"] + "/{Sample_ID}_1_fastqc.html", Sample_ID=SAMPLES),
        expand(config["results"]["fastqc_raw"] + "/{Sample_ID}_1_fastqc/fastqc_data.txt", Sample_ID=SAMPLES),
        # Trimmed 
        expand(config["results"]["trimmed"] + "/{Sample_ID}_P1.fastq.gz", Sample_ID=SAMPLES),
        expand(config["results"]["trimmed"] + "/{Sample_ID}_P2.fastq.gz", Sample_ID=SAMPLES),
        # FastQC trimmed
        expand(config["results"]["fastqc_trimmed"] + "/{Sample_ID}_P1_fastqc.zip", Sample_ID=SAMPLES),
        expand(config["results"]["fastqc_trimmed"] + "/{Sample_ID}_P1_fastqc.html", Sample_ID=SAMPLES),
        expand(config["results"]["fastqc_trimmed"] + "/{Sample_ID}_P1_fastqc/fastqc_data.txt", Sample_ID=SAMPLES),

        # --- 2. ALIGMENT  ---
        # INDEXING: BWA, Samtools faidx (.fna.fai), GATK/Picard Dictionary (.dict)
        REF + ".bwt",
        REF + ".pac",
        REF + ".ann",
        REF + ".amb",
        REF + ".sa",
        REF + ".fai",
        REF_NO_EXT + ".dict",

        # # MAPPING 
        #expand(config["results"]["bams"] + "/{Sample_ID}/{Sample_ID}.sam", Sample_ID=SAMPLES),
        #expand(config["results"]["bams"] + "/{Sample_ID}/{Sample_ID}.sorted.bam", Sample_ID=SAMPLES),
        expand(config["results"]["bams"] + "/{Sample_ID}/{Sample_ID}.sorted.markdup.bam", Sample_ID=SAMPLES),
        expand(config["results"]["bams"] + "/{Sample_ID}/{Sample_ID}.dup_metrics.txt", Sample_ID=SAMPLES),
        expand(config["results"]["bams"] + "/{Sample_ID}/{Sample_ID}.sorted.markdup.bam.bai", Sample_ID=SAMPLES),
        expand(config["results"]["stats"] + "/{Sample_ID}/{Sample_ID}.coverage.txt", Sample_ID=SAMPLES),
        expand(config["results"]["stats"] + "/{Sample_ID}/{Sample_ID}.mapping_flagstats.txt", Sample_ID=SAMPLES),

        # # PLOIDY ESTIMATION 
        expand(config["results"]["nquire"] + "/{Sample_ID}.bin", Sample_ID=SAMPLES),
        expand(config["results"]["nquire"] + "/{Sample_ID}_denoised.bin", Sample_ID=SAMPLES),
        expand(config["results"]["nquire"] + "/{Sample_ID}_lrdmodel.txt", Sample_ID=SAMPLES),
        expand(config["results"]["nquire"] + "/{Sample_ID}_histotest.txt", Sample_ID=SAMPLES),
        
        # --- 3. VARIANT CALLING 
        #expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.raw.vcf", Sample_ID=SAMPLES),
        expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.filtered.vcf", Sample_ID=SAMPLES),
        expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.snps.vcf.gz",Sample_ID=SAMPLES),
        expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.indels.vcf.gz", Sample_ID=SAMPLES),
        expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.filtered.vcf.gz", Sample_ID=SAMPLES),
        expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.snps.pass.hom.vcf.gz", Sample_ID=SAMPLES),
        expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.snps.pass.hom.vcf.gz.tbi",Sample_ID=SAMPLES),
        #expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.genome.fasta", Sample_ID=SAMPLES),
        #expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.genome.renamed.fasta", Sample_ID=SAMPLES),
        config["results"]["vcf"] + "/sites.reference.genome.bed",

        expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.ind.snps.bad.vcf.gz", Sample_ID=SAMPLES),
        expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.ind.snps.bad.bed", Sample_ID=SAMPLES),
        config["results"]["vcf"] + "/good_sites.global.bed",
        config["results"]["vcf"] + "/all_bad_sites.bed",
        config["results"]["vcf"] + "/ref.seq.fasta",
        #expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.genome.clean.fasta", Sample_ID=SAMPLES),
        #expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.genome_full.fasta", Sample_ID=SAMPLES),
        config["results"]["vcf"] + "/Allsamples.genome_full.fasta",
        config["results"]["phylo"] + "/Allsamples_genome.contree",
        config["results"]["phylo"] + "/Allsamples_genome.log",
        config["results"]["phylo"] + "/Allsamples_genome.iqtree",
        config["results"]["phylo"] + "/Allsamples_genome.nwk",
        
        # 4 - Variant statistics
        expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.pass_hom.txt", Sample_ID=SAMPLES),
        expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.pass_het.txt", Sample_ID=SAMPLES),
        expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.variant_counts.txt", Sample_ID=SAMPLES),
        expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.variant_stats.txt", Sample_ID=SAMPLES),
        expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.total_snps.txt", Sample_ID=SAMPLES),
        expand(config["results"]["vcf"] + "/{Sample_ID}/{Sample_ID}.total_indels.txt", Sample_ID=SAMPLES),
        config["results"]["base"] + "/updated_metadata.xlsx",

        config["results"]["snp_dists"] + "/Allsamples_genome_snp_dist.tsv"




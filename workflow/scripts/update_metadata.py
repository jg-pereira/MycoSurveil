#update_metadata.py

import os
import pandas as pd
import argparse

def mapping_flagstats_reading(mapping_flagstats):
    
    total_reads_trimmed = {}
    total_mapped_line = {}
    for file in mapping_flagstats:
        with open(file, "r") as f:
            lines = f.readlines()

            reads_trimmed = lines[1].strip().split(" + ")[0]
            mapped_line = lines[6].strip().split("mapped (")[1].split(" : N/A)")[0]
            parent_dir = os.path.abspath(file)
            sample_folder_name = os.path.basename(parent_dir)
            sample_name = sample_folder_name.split('.mapping_flagstats.txt')[0] 
            total_reads_trimmed[sample_name]= reads_trimmed
            total_mapped_line[sample_name]= mapped_line

    
    return total_reads_trimmed, total_mapped_line

# mapping_flagstats=["/mnt/c/Users/joana.gomes/results/stats/F/F.mapping_flagstats.txt" ]
# total_reads_trimmed, total_mapped_line = mapping_flagstats_reading(mapping_flagstats)
# print(f"samples_mapping_data: {total_mapped_line}")
# print(f"samples_reads_trimmed: {total_reads_trimmed}")

def flagstats_coverage(coverage):

    depth_coverage={}
    coverage_data={}
    for file in coverage:
        df=pd.read_csv(file, sep="\t")
        df_filtered = df.iloc[:-1]
        mean_coverage = df_filtered['coverage'].mean()
        #print(df_filtered['coverage'])
        meandepth = df_filtered['meandepth'].mean()
        #print(df_filtered['meandepth'])

        parent_dir = os.path.abspath(file)
        sample_folder_name = os.path.basename(parent_dir)
        sample_name = sample_folder_name.split('.coverage.txt')[0] 
        depth_coverage[sample_name]= str(meandepth)
        coverage_data[sample_name]= str(mean_coverage)
    
    return depth_coverage, coverage_data

#coverage = ["/mnt/c/Users/joana.gomes/MycoSurveil/results/stats/chr_A_CN10_2x/chr_A_CN10_2x.coverage.txt"]

#depth_coverage, coverage_data = flagstats_coverage(coverage)
#print(f"samples_depth_coverage: {depth_coverage}")
#print(f"samples_coverage_data: {coverage_data}")

def nquire_histotest_reading(histotest):
    ploidy_samples = {}
    rquire = {}

    for file in histotest:

        with open(file, "r") as f:
            lines = [line.strip() for line in f]

            diploid_raw = lines[3].split("r^2: ")[-1]
            triploid_raw = lines[8].split("r^2: ")[-1]
            tetraploid_raw = lines[13].split("r^2: ")[-1]
    
            diploid_value = float(diploid_raw)
            triploid_value = float(triploid_raw)
            tretraploid_value = float(tetraploid_raw)

            ploidy = "NA"
            r_quire = "0"

            if diploid_value > triploid_value and diploid_value > tretraploid_value:
                ploidy = "2"
                r_quire = diploid_value
            elif triploid_value > diploid_value and triploid_value > tretraploid_value:
                ploidy = "3"
                r_quire = triploid_value
            elif tretraploid_value > diploid_value and tretraploid_value > triploid_value:
                ploidy = "4"
                r_quire = tretraploid_value

            parent_dir = os.path.abspath(file)
            sample_folder_name = os.path.basename(parent_dir)
            sample_name = sample_folder_name.split('_histotest.txt')[0] 
            ploidy_samples[sample_name]= ploidy
            rquire[sample_name]= str(r_quire)
    
    return ploidy_samples, rquire

# histotest=["/mnt/c/Users/joana.gomes/results/nquire/F_histotest.txt"]
# ploidy_samples, rquire = nquire_histotest_reading(histotest)
# print(f"samples_ploidy: {ploidy_samples}")
# print(f"samples_rquire_data: {rquire}")

def fastqc_reading(fastqc_file):
    fastqc_raw_data = {}
    qc_raw_data = {}
    fastqc_trimmed_data = {}
    qc_trimmed_data = {}
    
    for file in fastqc_file:
        with open(file, "r") as f:
            lines = f.readlines()
            line7 = lines[6].strip()
            gc = lines[10].strip()
            
            int_gc = gc.split("\t")[-1]
            parts = line7.split("\t")
            nr_reads = int(parts[-1]) * 2 # Assuming paired-end data
            
            parent_dir = os.path.dirname(file)
            sample_folder_name = os.path.basename(parent_dir)

            if "_trimmed" in parent_dir:
                sample_name = sample_folder_name.replace("_P1_fastqc", "")
                fastqc_trimmed_data[sample_name] = nr_reads
                qc_trimmed_data[sample_name] = int_gc
            else:
                sample_name = sample_folder_name.replace("_1_fastqc", "") 
                fastqc_raw_data[sample_name] = nr_reads
                qc_raw_data[sample_name] = int_gc

    return fastqc_raw_data, qc_raw_data, fastqc_trimmed_data, qc_trimmed_data


# fastqc_file_raw = ["/mnt/c/Users/joana.gomes/results/fastqc_raw/F_1_fastqc/fastqc_data.txt"]
# fastqc_file_trimmed = ["/mnt/c/Users/joana.gomes/results/fastqc_trimmed/F_P1_fastqc/fastqc_data.txt"]
# raw_counts, raw_gc, _, _ = fastqc_reading(fastqc_file_raw)
# _, _, trimmed_counts, trimmed_gc = fastqc_reading(fastqc_file_trimmed)
# print(f"fastqc_raw_data: {raw_counts}")
# print(f"qc_raw_data: {raw_gc}")
# print(f"fastqc_trimmed_data: {trimmed_counts}")
# print(f"qc_trimmed_data: {trimmed_gc}")

def variant_counts(lista_arquivos):
    samples_vcf_data = {}
   
    for file in lista_arquivos:
        with open(file, "r") as f:
            lines=f.readlines()
            num_variants = lines[0]

        parent_dir = os.path.dirname(file)
        sample_folder_name = os.path.basename(parent_dir)
        sample_name= sample_folder_name.split('.variant_counts.txt')[0]
        samples_vcf_data[sample_name] = num_variants
    
    return samples_vcf_data

# sample_filter = ["/mnt/c/Users/joana.gomes/results/vcf/F/F.variant_counts.txt"]
# samples_vcf_data = variant_counts(sample_filter)
# print(f"Nr variants: {samples_vcf_data}")


def nr_homozigotes(files):
    samples_homozygoty_data = {}

    for file in files:
        with open(file, "r") as f:
            lines=f.readlines()
            nr_homozygoty = lines[0].strip()

        parent_dir = os.path.dirname(file)
        sample_folder_name = os.path.basename(parent_dir)
        sample_name= sample_folder_name.split('.pass_homo.txt')[0]
        samples_homozygoty_data[sample_name] = nr_homozygoty

    return samples_homozygoty_data
# files = ["/mnt/c/Users/joana.gomes/results/vcf/F/F.pass_homo.txt"]
# samples_homozygoty_data= nr_homozigotes(files)
# print(f"Nr pass homo: {samples_homozygoty_data}")

def nr_heterozigotes(files):
    samples_heterozigoty_data = {}

    for file in files:
        with open(file, "r") as f:
            lines = f.readlines()
            nr_heterozygoty = lines[0].strip() 

        parent_dir = os.path.dirname(file)
        sample_folder_name = os.path.basename(parent_dir) 
        sample_name = sample_folder_name.split('.pass_hetero.txt')[0]
        samples_heterozigoty_data[sample_name] = nr_heterozygoty

    return samples_heterozigoty_data

# files = ["/mnt/c/Users/joana.gomes/results/vcf/F/F.pass_hetero.txt"]
# samples_heterozigotes_data = nr_heterozigotes(files)
# print(f"Nr pass hetero: {samples_heterozigotes_data}")

def snps_reading(sample_snps):

    samples_snps_data = {}
    
    for file in sample_snps:
        nr = 0
        nr_pass = 0
        
        with open(file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue  
                nr += 1
                if "PASS" in line:
                    nr_pass += 1
        parent_dir = os.path.abspath(file)
        sample_folder_name = os.path.basename(parent_dir)
        sample_name = sample_folder_name.split('.snps.vcf')[0]
        samples_snps_data[sample_name]= (nr, nr_pass)
    
    return samples_snps_data

# sample_snps = ["/mnt/c/Users/joana.gomes/results/vcf/F/F.snps.vcf"]
# samples_snps_data = snps_reading(sample_snps)
# print(f"Nr snps: {samples_snps_data}")

def indels_reading(sample_indels):
    samples_indels_data = {}
    
    for file in sample_indels:
        nr = 0
        nr_pass = 0
        
        with open(file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue  
                nr += 1
                if "PASS" in line:
                    nr_pass += 1
        parent_dir = os.path.abspath(file)
        sample_folder_name = os.path.basename(parent_dir)
        sample_name = sample_folder_name.split('.indels.vcf')[0]
        samples_indels_data[sample_name]= (nr, nr_pass)
    
    return samples_indels_data

# sample_indels = ["/mnt/c/Users/joana.gomes/results/vcf/F/F.indel.vcf"]
# samples_indels_data = indels_reading(sample_indels)
# print(f"Nr indels: {samples_indels_data}")

def update_metadata(metadata,output_dir,raw_counts, raw_gc,total_reads_trimmed,trimmed_counts, trimmed_gc,total_mapped_line,ploidy_samples, rquire, depth_coverage, coverage_data, samples_vcf_data, samples_snps_data, samples_indels_data, samples_homozygoty_data, samples_heterozigotes_data):

    df = pd.read_excel(metadata)

    df['Raw_Reads'] = df['Sample_ID'].map(raw_counts)
    df['Raw_GC_Content'] = df['Sample_ID'].map(raw_gc)
    df['Trimmed_Reads'] = df['Sample_ID'].map(total_reads_trimmed)
    df['Trimmed_Reads'] = df['Sample_ID'].map(trimmed_counts)
    df['Trimmed_GC_Content'] = df['Sample_ID'].map(trimmed_gc)
    df['Mapping_Rate'] = df['Sample_ID'].map(total_mapped_line)
   
    df['Ploidy'] = df['Sample_ID'].map(ploidy_samples)
    df['R^2_Nquire'] = df['Sample_ID'].map(rquire)
    df['Depth_Coverage'] = df['Sample_ID'].map(depth_coverage)
    df['Coverage'] = df['Sample_ID'].map(coverage_data)
    df['Nr_Variants'] = df['Sample_ID'].map(samples_vcf_data)

    mapa_snps = df['Sample_ID'].map(samples_snps_data) 
    df['SNPs_Total'] = mapa_snps.str[0]
    df['SNPs_Filtered'] = mapa_snps.str[1]

    mapa_indels = df['Sample_ID'].map(samples_indels_data)
    df['Indels_Total']    = mapa_indels.str[0]
    df['Indels_Filtered'] = mapa_indels.str[1]

    df['High-Quality Homozygous Variants (Pass_Homo_SNP)'] = df['Sample_ID'].map(samples_homozygoty_data)
    df['High-Quality Heterozygous Variants (Pass_Hetero_SNP)'] = df['Sample_ID'].map(samples_heterozigotes_data)

    output_metadata_path = os.path.join(output_dir, "updated_metadata.xlsx")
    df.to_excel(output_metadata_path, index=False)
    print(f"Updated metadata saved to: {output_metadata_path}")

# metadata = "/mnt/c/Users/joana.gomes/resources/metadata.xlsx"
# output_dir = "/mnt/c/Users/joana.gomes/results/"


def main():

    parser = argparse.ArgumentParser(description = "Updates metadata with FastQC results.")
    parser.add_argument("--metadata", required = True, help = "Path to the input metadata Excel file.")
    parser.add_argument("--mapping_flagstats", required = True, nargs = '+', help = "Paths to the mapping flagstats files.")
    parser.add_argument("--coverage", required = True, nargs = '+', help = "Paths to the coverage files.")
    parser.add_argument("--ploidy", required = True, nargs = '+', help = "Paths to the histotest files.")
    parser.add_argument("--fastqc_raw", required = True, nargs = '+', help = "Paths to the FastQC raw data files (fastqc_data.txt).")
    parser.add_argument("--fastqc_trimmed", required = True, nargs = '+', help = "Paths to the FastQC trimmed data files (fastqc_data.txt).")
    parser.add_argument("--nr_variants", required = True, nargs = '+', help = "Paths to the filtered VCF files.")
    parser.add_argument("--snps", required = True, nargs = '+', help = "Paths to the SNPs VCF files.")
    parser.add_argument("--indels", required = True, nargs = '+', help = "Paths to the Indels VCF files.")
    parser.add_argument("--output_metadata", required = True, help = "Path to the ouput folder to save the results.")
    parser.add_argument("--pass_homo", required = True, nargs = '+', help = "Paths to the pass homo files.")
    parser.add_argument("--pass_hetero", required = True, nargs = '+', help = "Paths to the pass hetero files.")

    args = parser.parse_args()

    metadata = args.metadata
    output_dir = args.output_metadata
    fastqc_file_raw = args.fastqc_raw
    fastqc_file_trimmed = args.fastqc_trimmed
    mapping_flagstats = args.mapping_flagstats
    histotest = args.ploidy
    coverage = args.coverage
    sample_filter = args.nr_variants
    sample_snps = args.snps
    sample_indels = args.indels
    files_homo = args.pass_homo
    files_hetero = args.pass_hetero

    total_reads_trimmed, total_mapped_line = mapping_flagstats_reading(mapping_flagstats)
    depth_coverage, coverage_data = flagstats_coverage(coverage)
    ploidy_samples, rquire = nquire_histotest_reading(histotest)
    raw_counts, raw_gc, _, _ = fastqc_reading(fastqc_file_raw)
    _, _, trimmed_counts, trimmed_gc = fastqc_reading(fastqc_file_trimmed)  
    samples_vcf_data = variant_counts(sample_filter)
    samples_snps_data = snps_reading(sample_snps)
    samples_indels_data = indels_reading(sample_indels)
    samples_homozygoty_data= nr_homozigotes(files_homo)
    samples_heterozigotes_data = nr_heterozigotes(files_hetero)
    update_metadata(metadata,output_dir,raw_counts, raw_gc,total_reads_trimmed,trimmed_counts, trimmed_gc,total_mapped_line,ploidy_samples, rquire, depth_coverage, coverage_data, samples_vcf_data, samples_snps_data, samples_indels_data, samples_homozygoty_data, samples_heterozigotes_data)

if __name__ == "__main__":
    main()
# MycoSurveil (under development)

**MycoSurveil** is an open-source, start-to-end bioinformatics framework tailored for the genomic surveillance of fungal pathogens. Developed to address the critical gaps in fungal genomics, the project provides automated solutions to detect, characterize, and track public health threats caused by *Candida*.

This project is hosted and developed at the **National Health Institute Dr. Ricardo Jorge (INSA)**.

## Public Health Impact
By integrating clinical, epidemiological, and genomic data, MycoSurveil aims to:

- Enhance the preparedness of reference laboratories to tackle fungal outbreaks.
- Accelerate the technological transition toward genomics-based fungal surveillance.
- Provide an immediate impact on national and global capacity to fight fungal priority pathogens (WHO).

## Implementation
The MycoSurveil pipeline is a command-line tool, implemented using Snakemake, with modules tailored for fungal genomic surveillance:

1. **Quality Control & Trimming:** Automated raw data processing.
2. **Read Mapping and Ploidy Estimation:** Robust assessment of genome ploidy levels and aneuploidies.
3. **Variant Calling and Phylogenetic Analysis:** Accurate detection of Single Nucleotide Polymorphisms (SNPs) and Indels against reference fungal genomes. Reconstruction of evolutionary relationships and high-resolution automated cluster detection to identify outbreaks and transmission chains.

## Input
To run this pipeline, an input folder with a specific structure is required.
For instance, you can create a folder named `resources` containing your data. Crucially, this folder must be located at the same directory level as the cloned MycoSurveil repository.

### Required Directory Structure:
```
parent_directory/
├── MycoSurveil/               # Cloned GitHub repository
└── resources/                 # Input folder (at the same level as MycoSurveil)
    ├── rawdata_fastq/         # Subfolder containing reads in FASTQ format
    ├── ref_genome/            # Subfolder containing the reference genome (*Candida*) in .fna format
    └── metadata.xlsx          # Excel file containing two mandatory columns: Sample_ID and NGS_ID
```
(Note: A ready-to-use toy dataset following this exact structure is provided in the `examples/` folder for testing and demonstration purposes).

## Output
- A metadata `.xlsx` file containing the following variables: 
  `Sample_ID`, `NGS_ID`, `Raw_Reads`, `Raw_GC_Content`, `Trimmed_Reads`, `Trimmed_GC_Content`, `Mapping_Rate`, `Ploidy`, `R^2_Nquire`, `Depth_Coverage`, `Coverage`, `Nr_Variants`, `SNPs_Total`, `SNPs_Filtered`, `Indels_Total`, `Indels_Filtered`, `High-Quality Homozygous Variants (Pass_Homo_SNP)`, and `High-Quality Heterozygous Variants (Pass_Hetero_SNP)`.
- A phylogenetic tree file (newick format).


### Installation with conda

1. **Clone the repository**

```bash
git clone https://github.com/jg-pereira/MycoSurveil.git
```

2. **Create the Conda environment**

```bash
conda create -n mycosurveil -c conda-forge -c bioconda snakemake=9.13.7
```

3. **Activate the environment**

```bash
conda activate mycosurveil
cd MycoSurveil
```

### Examples of command-line usage
To run the pipeline using Snakemake with Conda integration (utilizing 4 cores as an example):
```bash
snakemake -s workflow/snakefile.smk -j 4 -p --use-conda
```


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
- A folder with two subfolders:
  - `rawdata_fastq`: containing reads in FASTQ format.
  - `ref_genome`: containing the reference genome of *Candida* in fna format.
- An Excel metadata file (.xlsx) containing two mandatory columns: Sample_ID and NGS_ID

## Output
- A metadata `.xlsx` file containing the following variables: 
  `Sample_ID`, `NGS_ID`, `Raw_Reads`, `Raw_GC_Content`, `Trimmed_Reads`, `Trimmed_GC_Content`, `Mapping_Rate`, `Ploidy`, `R^2_Nquire`, `Depth_Coverage`, `Coverage`, `Nr_Variants`, `SNPs_Total`, `SNPs_Filtered`, `Indels_Total`, `Indels_Filtered`, `High-Quality Homozygous Variants (Pass_Homo_SNP)`, and `High-Quality Heterozygous Variants (Pass_Hetero_SNP)`.
- A phylogenetic tree file (newick format).


### Installation with conda

#### 1. Manual installation from GitHub repository

```bash
git clone https://github.com/joana.gomes/MycoSurveil.git
```

### Examples of command-line usage
```bash
snakemake -s workflow/snakefile.smk -j 4 -p --use-conda
```


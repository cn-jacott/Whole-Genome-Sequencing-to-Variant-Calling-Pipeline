# Step 6. Chromosome-level filtering of high-quality variants

## Description

In this step, jointly genotyped chromosome-level VCF files are filtered to retain only **high-quality variants**. Filtering is performed **per chromosome** to reduce memory usage and improve scalability on HPC systems.

Because **SNPs** and **INDELs** require different filtering thresholds, the script allows selecting the variant type at the beginning and applies the appropriate filtering criteria automatically.

The workflow extracts the selected variant type, applies **GATK VariantFiltration**, and retains only variants with `PASS` status.

## Script summary

1.  Defines the chromosome to process  
2.  Selects the variant type (**SNPs** or **INDELs**) 
3.  Extracts the selected variant type from the input VCF
4.  Applies **GATK VariantFiltration** using variant-type-specific thresholds  
5.  Keeps only variants with `PASS` status
6.  Indexes the final filtered VCF
7.  Reports variant counts before and after filtering
    

## Inputs

-   Chromosome-level genotyped VCF  
    `../../working/vcfs/gmax/<GM>_genotyped_chunks_renamed_gmax.vcf.gz`
   
-   Reference genome  
    `../../ref/glyma.Wm82.gnm6.S97D.genome_main.fna`
    

## Outputs

Written to
`../../working/vcfs/basic_gatk_filtered/`

Output files depend on the selected variant type:
SNP output: `<GM>_snps_PASS.vcf.gz`
INDEL output:`<GM>_indels_PASS.vcf.gz`

Example for chromosome 1:
`Gm01_snps_PASS.vcf.gz`  
`Gm01_indels_PASS.vcf.gz`

## Filtering criteria

Filtering thresholds depend on the selected variant type. Only variants with `PASS` status are retained in the final output VCF.

### SNP filtering
Variants are labeled using: `SNP_filter`
with the following expression: `QD < 26.0 || FS > 60.0 || MQ < 40.0`

### INDEL filtering

Variants are labeled using: `INDEL_filter`
with the following expression: `QD < 26.0 || FS > 200.0 || MQ < 40.0`

Only variants with `PASS` status are retained in the final output VCF.

## Tools / modules

-   **GATK v4.5.0.0**   
-   **BCFtools v1.18** 
-   **HTSlib v1.18**
-   **SLURM workload manager**
    
Modules loaded within the script:
`module load GATK/4.5.0.0-GCCcore-12.3.0-Java-17`  
`module load BCFtools/1.18-GCC-12.3.0`  
`module load HTSlib/1.18-GCC-12.3.0`

## How to run on the cluster (SLURM)

### 1) Set the chromosome and variant type
```text
GM="Gm01"  
VARIANT_TYPE="snps"  # options: snps or indels
```
### 2) Submit the job

`sbatch filter_variants_Gm01.sh`

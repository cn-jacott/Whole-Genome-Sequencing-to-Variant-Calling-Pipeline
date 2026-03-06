# Step 7. Concatenation of chromosome-level VCF files

## Description

In this step, chromosome-level VCF files are concatenated into a single genome-wide VCF file using **BCFtools concat**. This step merges the filtered variant datasets from all chromosomes into one file suitable for downstream analyses.

The script reads chromosome-specific VCF files in numeric order (`Gm01–Gm20`) and concatenates them into a compressed and indexed genome-wide VCF.

Concatenation is performed after chromosome-level filtering to improve performance and reduce memory usage.

## Script summary

1.  Defines the input directory containing chromosome-level VCF files
2.  Checks that all expected chromosome VCFs (`Gm01–Gm20`) are present
3.  Concatenates all chromosome VCF files using **bcftools concat**
4.  Compresses and indexes the final genome-wide VCF
5.  Reports the number of variants in the final file
    

## Inputs

Chromosome-level VCF files:
```text
../../working/vcfs/Gm01_snps_PASS.vcf.gz 
../../working/vcfs/Gm02_snps_PASS.vcf.gz  
...  
../../working/vcfs/Gm20_snps_PASS.vcf.gz 
```

## Outputs

Genome-wide VCF file:`../../working/vcfs/FINAL_genome_level.vcf.gz`

The final VCF is **bgzip-compressed and indexed**.

## Tools / modules

-   **BCFtools v1.18** 
-   **HTSlib v1.18**
-   **SLURM workload manager**
    

Modules loaded within the script:

`module load BCFtools/1.18-GCC-12.3.0 ` 

`module load HTSlib/1.18-GCC-12.3.0`

## How to run on the cluster (SLURM)

Submit the job:
`sbatch Step7_Merge_genome_level_VCF.sh`


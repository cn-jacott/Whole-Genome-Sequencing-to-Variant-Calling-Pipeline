## Step 5. Joint genotyping with GenotypeGVCFs

## Description

In this step, each **GenomicsDB chunk** generated in **Step 4b** is genotyped with **GATK GenotypeGVCFs** to produce one VCF per interval. These chunk-level VCFs are then concatenated in numeric order to generate a final chromosome-level VCF.

The script performs two sequential jobs:

-   a **SLURM array job** to genotype each GenomicsDB interval independently
-   a **dependent gather job** to concatenate all chunk VCFs into a single chromosome level VCF

## Script summary

1.  Defines the chromosome and chunk size
2.  Reads the interval list generated in Step 4a  
3.  Submits a **GenotypeGVCFs** SLURM array job, with one task per interval
4.  Writes one compressed VCF per chunk 
5.  Submits a dependent job to concatenate all chunk VCFs in numeric order
6.  Indexes the final chromosome VCF  
7.  Reports per-chunk and final variant counts

## Inputs

-   GenomicsDB workspaces from Step 4b `../working/genomicsDB/genomicsdb_<GM>/gendb_<GM>_<START>_<END>/`
    
-   Interval list from Step 4a  
`../working/genomicsDB/genomicsdb_<GM>/<GM>_<CHUNK>.intervals.list`
    
-   Reference genome  
`../ref/glyma.Wm82.gnm6.S97D.genome_main.fna`

## Outputs

### Chunk-level VCFs

Written to:
`../working/joint_genotype/genotyped_chunks_<GM>_<CHUNK>/`

Example:
```text
Gm01_chunk_1.vcf.gz  
Gm01_chunk_2.vcf.gz  
Gm01_chunk_3.vcf.gz
```
### Final chromosome VCF

Written to:
`../working/joint_genotype/<GM>_genotyped_chunks.vcf.gz`

Example:
`../working/joint_genotype/Gm18_genotyped_chunks.vcf.gz`

## Tools / modules

-   **GATK v4.5.0.0**
-   **BCFtools v1.18**
-   **SLURM workload manager**

Modules loaded within the submitted jobs:

`module load GATK/4.5.0.0-GCCcore-12.3.0-Java-17`  
`module load BCFtools/1.18-GCC-12.3.0`

## How to run on the cluster (SLURM)

### 1) Edit chromosome-specific settings
```text
GM="Gm01"  
CHUNK=10000000  
MAX_CONCURRENT=6
```
### 2) Submit the wrapper script
`bash genotype_genDB_chunks_order.sh`

This script automatically submits:
-   one **array job** for `GenotypeGVCFs`    
-   one **dependent gather job** for `bcftools concat`
    

## Notes

-   Each array task processes **one interval** from the interval list
-   Output chunk VCFs are renamed in simple numeric order
-   The final concatenation job checks that the sum of variants across chunks matches the final merged VCF

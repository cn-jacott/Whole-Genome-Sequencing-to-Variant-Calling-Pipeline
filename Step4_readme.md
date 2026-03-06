
## ## Step 4. GenomicsDBImport per chromosome

## Description

In this step, per-sample gVCF files generated in the previous step are imported into a **GenomicsDB workspace** using **GATK GenomicsDBImport**. This is a intermediate step before joint genotyping with `GATK GenotypeGVCFs`.
To reduce memory usage and improve parallelization on HPC systems, the import is performed **per chromosome** and split into **fixed-size chunks**. For this reason, Step 4 is divided into two parts:

-   **Step 4a** – prepare the **sample map** and **interval list**
-   **Step 4b** – run `GenomicsDBImport` per interval using a **SLURM array job**

---

##  ## Step 4a. Prepare sample map and interval list

## Description

This step prepares the input files required by `GenomicsDBImport` for one chromosome.

### Script summary

1.  Defines the chromosome to process

3.  Defines the chunk size
    
4.  Builds a **sample map** from chromosome-specific gVCF files
    
5.  Retrieves chromosome length from the reference `.fai`
    
6.  Generates a list of genomic **intervals**
    

### Inputs
-   Chromosome-specific gVCFs   `../working/<GM>/*.g.vcf.gz`
    
-   Reference genome    `../ref/glyma.Wm82.gnm6.S97D.genome_main.fna`
    

### Outputs

Written to:
`../working/genomicsDB/genomicsdb_<GM>/`
- sample map
`sample_map_<GM>.txt`

- intervals list
`<GM>_<CHUNK>.intervals.list`
    

### Sample map format

`sample_name    /path/to/sample.g.vcf.gz`

Example:

```text
sample_1    ../working/Gm18/sample_1.g.vcf.gz
sample_2    ../working/Gm18/sample_2.g.vcf.gz
sample_3    ../working/Gm18/sample_3.g.vcf.gz
```

## Tools / modules

-   **Python** – interval generation
    
-   **SLURM workload manager**

## How to run on the cluster (SLURM)

`sbatch prep_gendb_Gm18.sh`

---
## Step 4b 
### Script summary

1.  Defines the chromosome and chunk size
    
2.  Reads the **sample map**
    
3.  Reads the **interval list**
    
4.  Assigns one interval per SLURM task
    
5.  Runs `GenomicsDBImport`
6.  Creates one GenomicsDB workspace per interval

### Inputs
-   Sample map and intervals list from Step 4a  
    `sample_map_<GM>.txt`
    `<GM>_<CHUNK>.intervals.list`


-   gVCFs from Step 3 
-   Reference genome
      `../ref/glyma.Wm82.gnm6.S97D.genome_main.fna`

### Outputs

Written to:
`../working/genomicsDB/genomicsdb_<GM>/`

Example workspace:
`gendb_Gm18_1_10000000/`

## Tools / modules

-   **GATK v4.5.0.0**
-   **SLURM workload manager**
    

Modules loaded in the script:
`module load GATK/4.5.0.0-GCCcore-12.3.0-Java-17`

## How to run on the cluster (SLURM)

### 1) Set the array size

The number of tasks must match the number of intervals:
`#SBATCH --array=1-6`

### 2) Submit the job 
`sbatch genomicsdb_import_Gm18.sh`

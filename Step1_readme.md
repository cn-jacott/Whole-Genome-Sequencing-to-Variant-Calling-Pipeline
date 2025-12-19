
# Step 1. Quality control of raw sequencing reads (FastQC)

## Description

This step performs quality control (QC) on raw paired-end Illumina sequencing reads using **FastQC**. All `*.fq.gz` files found within the input data directory are processed in parallel to assess read quality, base composition, GC content, and potential adapter contamination.

The script is designed for execution on an HPC cluster using **SLURM** and uses **GNU Parallel** to efficiently process large numbers of FASTQ files within a single job.

FastQC output files are renamed to include the parent sample folder name, ensuring unique and traceable filenames when multiple samples or sequencing lanes are present.

## Input data structure

The script assumes the following directory structure:

<pre><code>working/
|-- data/
|   |-- 1/
|   |   |-- sample1_lane1_1.fq.gz
|   |   |-- sample1_lane1_2.fq.gz
|   |-- 2/
|   |   |-- sample2_lane1_1.fq.gz
|   |   |-- sample2_lane1_2.fq.gz
|   |   |-- sample2_lane2_1.fq.gz
|   |   |-- sample2_lane2_2.fq.gz
|   |-- ...
</code></pre>

Each subdirectory corresponds to a single sample and may contain one or more sequencing lanes.



## Output

All FastQC results are written to:

`working/fastqc/` 

For each FASTQ file, the following outputs are generated:

-   `*_fastqc.html`
-   `*_fastqc.zip`
    
Files are renamed using the format:
`{sample_folder}_{fastq_filename}_fastqc.{html,zip}` 

Example: 
`1_sample1_lane1_1_fastqc.html` 



## Tools and modules
-   FastQC v0.12.1
-   GNU Parallel
-   SLURM workload manager
   
Modules are loaded within the script:
`module load parallel/20230722-GCCcore-12.2.0`
`module load FastQC/0.12.1-Java-11`

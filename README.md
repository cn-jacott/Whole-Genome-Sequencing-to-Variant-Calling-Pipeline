
# Whole-Genome Sequencing to Variant Calling Pipeline (FASTQC → VCF)

##  Overview
This repository describes a reproducible pipeline for whole-genome sequencing data processing and variant calling, starting from raw Illumina FASTQ files and ending with jointly genotyped VCF files.
The workflow follows GATK best practices and is designed for execution on high-performance computing (HPC) systems using SLURM.

The pipeline supports large sequencing panels by performing per-sample variant calling in gVCF mode, followed by chromosome-wise joint genotyping. It is suitable for plant genomes and has been tested on soybean (*Glycine max* and *Glycine soja*) whole-genome data.

##  Pipeline overview
**FASTQ → BAM → gVCF → GenomicsDB → VCF**

An overview of each step is provided below. Detailed commands and execution scripts for each step are available in the corresponding Markdown files.
Before starting the pipeline, the reference genome for the species of interest must be indexed. Instructions for reference genome indexing are provided in the Markdown file `Reference_genome_indexing.md`

##  Step 1. Quality control of raw reads
Raw paired-end FASTQ files are inspected to assess read quality, adapter contamination, and base composition.

**Input:** Paired-end read files (`*.fastq.gz`)
**Output:** QC reports (FastQC, optionally MultiQC)
**Tools:** FastQC

## Step 2. Read alignment to the reference genome
Paired-end reads are aligned to the reference genome, merged across lanes per sample (if necessary), deduplicated, indexed to produce final BAM files.

**Input:** Paired-end read files (`*.fastq.gz`), indexed reference genome (`*.fna`)
**Output:** Sorted and indexed deduplicated BAM files (`*.bam`, `*.bai`)
**Tools:** BWA-MEM, SAMtools, Picard

## Step 3. HaplotypeCaller gVCF per sample per chromosome
Variants are called independently for each sample using GATK HaplotypeCaller in gVCF mode. Variant calling is performed per chromosome to improve scalability on HPC systems.

**Input:** Sorted and indexed deduplicated BAM files (`*.bam`, `*.bai`), indexed reference genome (`*.fna`)  
**Output:** Per-sample, per-chromosome gVCF files (`*.g.vcf.gz`, `*.tbi`)  
**Tools:** GATK `HaplotypeCaller`

## Step 4. GenomicsDBImport per chromosome
Per-sample gVCF files are imported into a GenomicsDB workspace on a per-chromosome basis. This step enables efficient joint genotyping of large cohorts.

**Input:** Per-sample gVCF files, sample map file  
**Output:** GenomicsDB workspace (one per chromosome)  
**Tools:** GATK `GenomicsDBImport`


## Step 5. Joint genotyping with GenotypeGVCFs

Joint genotyping is performed on each chromosome using the corresponding GenomicsDB workspace to produce multisample VCF files.

**Input:** GenomicsDB workspace, indexed reference genome  
**Output:** Jointly genotyped multisample VCF files (`*.vcf.gz`, `*.tbi`)  
**Tools:** GATK `GenotypeGVCFs`


## Step 6. VCF merging and variant filtering

Jointly genotyped VCF files from individual chromosomes are merged to generate genome-wide VCF files. Variant-level filtering is applied to remove low-quality sites using hard filters or downstream filtering criteria, depending on the study design.

**Input:** Per-chromosome joint VCF files   (`*.vcf.gz`, `*.tbi`)  
**Output:** Genome-wide filtered VCF files   (`*.vcf.gz`, `*.tbi`)  
**Tools:** GATK, bcftoolses and filter

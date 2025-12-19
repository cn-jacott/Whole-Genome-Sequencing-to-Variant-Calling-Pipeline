
# Step 2. Read alignment and BAM generation (FASTQ → per sample BAM)

## Description
This step aligns paired-end reads to the reference genome with **BWA-MEM**, handles **multiple sequencing lanes per sample** using lane-specific read groups, merges lanes into a single BAM per sample, marks duplicates, and produces basic alignment/QC metrics. It is implemented as a **SLURM array job**, where each array task processes **one sample folder**.

## Script summary
1. **Align each lane independently** with BWA-MEM and add a lane-specific read group (`@RG`).
2. **Sort** alignments per lane and write a lane BAM.
3. **Merge lanes** to create one BAM per sample (or copy if only one lane).
4. **Mark duplicates** to create a final deduplicated BAM + index.
5. **Generate QC metrics** (WGS metrics, alignment summary, insert size).

## Expected input structure
The script assumes one folder per sample under `../working/data/`, with FASTQ files named as `*_1.fq.gz` and `*_2.fq.gz` (one or more lanes per sample).

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



## Inputs
- Paired-end reads: `../working/data/<SAMPLE>/*_1.fq.gz` and `*_2.fq.gz`
- Reference genome: `../ref/glyma.Wm82.gnm6.S97D.genome_main.fna`  
  (must be indexed for BWA; and have `.fai` + sequence dictionary for GATK tools)


## Outputs (written to `../working/bwa/`)
Per sample:
- Final deduplicated BAM: `<SAMPLE>_sorted_dedup.bam`
- BAM index: `<SAMPLE>_sorted_dedup.bam.bai`
- QC metrics:
  - `<SAMPLE>_wgs_metrics.txt`
  - `<SAMPLE>_alignment_metrics.txt`
  - `<SAMPLE>_insert_size.txt`
  - `<SAMPLE>_insert_size.pdf`

Intermediate files (removed by default):
- Per-lane BAMs: `<SAMPLE>_<RGID>_sorted.bam`
- Merged BAM: `<SAMPLE>_sorted.bam`



## Tools / modules
Loaded within the script:
- **BWA** `0.7.17`
- **SAMtools** `1.17`
- **GATK** `4.5.0.0`
-  **R** – required by GATK to generate insert size histogram PDFs (`CollectInsertSizeMetrics`)
    


## Script summary
1. **Align each lane independently** with BWA-MEM and add a lane-specific read group (`@RG`).
2. **Sort** alignments per lane and write a lane BAM.
3. **Merge lanes** to create one BAM per sample (or copy if only one lane).
4. **Mark duplicates** to create a final deduplicated BAM + index.
5. **Generate QC metrics** (WGS metrics, alignment summary, insert size).




## How to run on the cluster (SLURM)

### 1) Set the SLURM array range

Edit the SLURM header to match the number of sample directories in `../working/data/`:
`#SBATCH --array=1-n` 

To count sample directories:
`find ../working/data -mindepth 1 -maxdepth 1 -type d | wc -l` 


### 2) Submit the job
From the directory containing the script:
`sbatch bwa_lanes.sh` 

Job-specific log files will be generated for each array task:
-   `bwa_<jobID>_<taskID>.log`
-   `bwa_<jobID>_<taskID>.err`

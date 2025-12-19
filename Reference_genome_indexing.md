# Reference genome indexing (required before Step 2)

Before running the alignment step, the reference genome must be indexed for **BWA**, **SAMtools**, and **GATK**. These index files are required for read alignment, BAM processing, and all downstream variant calling steps.

The reference genome used in this pipeline is:

`../ref/glyma.Wm82.gnm6.S97D.genome_main.fna` 

Reference indexing only needs to be performed **once per reference genome**.



## 1) BWA index

Create the BWA index files required for read alignment:

`module load BWA/0.7.17-GCCcore-12.2.0

bwa index glyma.Wm82.gnm6.S97D.genome_main.fna` 

This command generates the following index files:

`.fna.amb
.fna.ann
.fna.bwt
.fna.pac
.fna.sa` 


## 2) SAMtools FASTA index

Create the FASTA index required by SAMtools and GATK:

`module load SAMtools/1.17-GCC-12.2.0

samtools faidx glyma.Wm82.gnm6.S97D.genome_main.fna` 

This generates:

`.fna.fai` 



## 3) GATK sequence dictionary

Create the sequence dictionary required by GATK tools:

`module load GATK/4.5.0.0-GCCcore-12.3.0-Java-17`

`gatk CreateSequenceDictionary -R glyma.Wm82.gnm6.S97D.genome_main.fna`

This generates:

`.dict`

#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH --output=bwa_%A_%a.log
#SBATCH --error=bwa_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1-10               # <-- set this to number of *sample folders*
#SBATCH --mail-type=END,FAIL

start_time=$(date +%s)
echo "Start time: $(date)"

# -------------------------
# Load tools
# -------------------------
module purge
module load R
module load BWA/0.7.17-GCCcore-12.2.0
module load SAMtools/1.17-GCC-12.2.0
module load GATK/4.5.0.0-GCCcore-12.3.0-Java-17

# -------------------------
# Paths
# -------------------------
ref="../ref/glyma.Wm82.gnm6.S97D.genome_main.fna"
DATA_DIR="../working/data"
OUT="../working/bwa"

mkdir -p "$OUT"

# -------------------------
# Sample list = one folder per sample (0, 2, 4, 15, etc.)
# -------------------------
mapfile -t SAMPLE_DIRS < <(find "$DATA_DIR" -mindepth 1 -maxdepth 1 -type d | sort)

INDEX=$((SLURM_ARRAY_TASK_ID - 1))

if [[ $INDEX -ge ${#SAMPLE_DIRS[@]} ]]; then
    echo "No sample folder for index $SLURM_ARRAY_TASK_ID" >&2
    exit 1
fi

SAMPLE_DIR="${SAMPLE_DIRS[$INDEX]}"
SAMPLE=$(basename "$SAMPLE_DIR")

echo "Processing SAMPLE: $SAMPLE"
echo "Sample folder: $SAMPLE_DIR"

# -------------------------
# Find all lanes for this sample
# -------------------------
mapfile -t R1_LANES < <(find "$SAMPLE_DIR" -type f -name "*_1.fq.gz" | sort)

if [[ ${#R1_LANES[@]} -eq 0 ]]; then
    echo "No R1 files found in $SAMPLE_DIR" >&2
    exit 1
fi

lane_bams=()

# -------------------------
# STEP 1: BWA MEM + sort per lane (each with its own RG)
# -------------------------
for R1 in "${R1_LANES[@]}"; do
    R2=${R1/_1.fq.gz/_2.fq.gz}

    if [[ ! -f "$R2" ]]; then
        echo "Missing R2 for $R1 (expected $R2)" >&2
        exit 1
    fi

    RGID=$(basename "$R1" _1.fq.gz)    # e.g. E250090877_L01_1_1
    LANEBAM="${OUT}/${SAMPLE}_${RGID}_sorted.bam"

    echo "  Lane RGID: $RGID"
    echo "    R1: $R1"
    echo "    R2: $R2"
    echo "    BAM: $LANEBAM"

    bwa mem -t 32 \
      -R "@RG\tID:${RGID}\tSM:${SAMPLE}\tPL:DNBSEQ\tLB:lib_${SAMPLE}\tPU:${RGID}" \
      "$ref" "$R1" "$R2" \
      | samtools sort -@ 32 -o "$LANEBAM" -

    lane_bams+=("$LANEBAM")
done

# -------------------------
# STEP 2: Merge lanes -> one BAM per sample
# -------------------------
MERGED_BAM="${OUT}/${SAMPLE}_sorted.bam"

if [[ ${#lane_bams[@]} -eq 1 ]]; then
    echo "Only one lane for $SAMPLE, copying BAM."
    cp "${lane_bams[0]}" "$MERGED_BAM"
else
    echo "Merging ${#lane_bams[@]} lanes for $SAMPLE"
    samtools merge -@ 32 "$MERGED_BAM" "${lane_bams[@]}"
fi

samtools index "$MERGED_BAM"

# (optional) remove per-lane BAMs to save space
rm -f "${lane_bams[@]}"

# -------------------------
# STEP 3: Mark duplicates -> FINAL BAM
# -------------------------
DEDUP_BAM="${OUT}/${SAMPLE}_sorted_dedup.bam"

gatk MarkDuplicatesSpark \
  -I "$MERGED_BAM" \
  -O "$DEDUP_BAM" \
  --create-output-bam-index true

rm -f "$MERGED_BAM" "$MERGED_BAM.bai"

# -------------------------
# STEP 4: QC Metrics
# -------------------------
gatk CollectWgsMetrics \
  -I "$DEDUP_BAM" \
  -O "${OUT}/${SAMPLE}_wgs_metrics.txt" \
  -R "$ref"

gatk CollectAlignmentSummaryMetrics \
  -R "$ref" \
  -I "$DEDUP_BAM" \
  -O "${OUT}/${SAMPLE}_alignment_metrics.txt"

gatk CollectInsertSizeMetrics \
  -I "$DEDUP_BAM" \
  -O "${OUT}/${SAMPLE}_insert_size.txt" \
  -H "${OUT}/${SAMPLE}_insert_size.pdf"

# -------------------------
# Done
# -------------------------
echo "Final BAM: $DEDUP_BAM"
echo "Final BAI: ${DEDUP_BAM}.bai"

end_time=$(date +%s)
echo "End time: $(date)"
echo "Total runtime: $((end_time - start_time)) seconds"

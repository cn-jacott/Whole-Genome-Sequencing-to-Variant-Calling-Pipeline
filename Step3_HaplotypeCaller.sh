#!/bin/bash
#SBATCH --job-name=gatk_haplotype
#SBATCH --output=gatk_haplotype_%A_%a.log
#SBATCH --error=gatk_haplotype_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=96:00:00
#SBATCH --array=1-20   # 20 chromosomes

# ====== PATHS YOU CAN EDIT ======
BAM_DIR="../working/bwa"
GVCF_DIR="../working/gvcf"
REFERENCE="../ref/glyma.Wm82.gnm6.S97D.genome_main.fna"
# ================================

module load GATK/4.5.0.0-GCCcore-12.3.0-Java-17
module load BCFtools/1.18-GCC-12.3.0

# 20 chromosomes — clean, fixed, consistent
declare -a regions=(
"glyma.Wm82.gnm6.Gm01:1-59600650"
"glyma.Wm82.gnm6.Gm02:1-54790513"
"glyma.Wm82.gnm6.Gm03:1-49119424"
"glyma.Wm82.gnm6.Gm04:1-55533332"
"glyma.Wm82.gnm6.Gm05:1-44449196"
"glyma.Wm82.gnm6.Gm06:1-53606008"
"glyma.Wm82.gnm6.Gm07:1-47205696"
"glyma.Wm82.gnm6.Gm08:1-50288313"
"glyma.Wm82.gnm6.Gm09:1-51135784"
"glyma.Wm82.gnm6.Gm10:1-55391814"
"glyma.Wm82.gnm6.Gm11:1-42152160"
"glyma.Wm82.gnm6.Gm12:1-43682070"
"glyma.Wm82.gnm6.Gm13:1-46106413"
"glyma.Wm82.gnm6.Gm14:1-52776670"
"glyma.Wm82.gnm6.Gm15:1-55366112"
"glyma.Wm82.gnm6.Gm16:1-40198974"
"glyma.Wm82.gnm6.Gm17:1-44073542"
"glyma.Wm82.gnm6.Gm18:1-60749638"
"glyma.Wm82.gnm6.Gm19:1-53875614"
"glyma.Wm82.gnm6.Gm20:1-51026854"
)

# Region for this array task
REGION="${regions[$((SLURM_ARRAY_TASK_ID-1))]}"

# Ensure BAMs exist
if ! ls "$BAM_DIR"/*_sorted_dedup*.bam 1>/dev/null 2>&1; then
    echo "No BAMs matching *_sorted_dedup*.bam found in $BAM_DIR"
    exit 1
fi

# Loop over BAMs
for bam in "$BAM_DIR"/*_sorted_dedup*.bam; do
    sample=$(basename "$bam" | sed 's/_sorted_dedup.*\.bam//')

    echo "[$(date)] Sample: $sample | Region: $REGION"

    OUT="$GVCF_DIR/${sample}_${REGION}.g.vcf.gz"

    gatk --java-options "-Xmx24g" HaplotypeCaller \
        -R "$REFERENCE" \
        -I "$bam" \
        -O "$OUT" \
        --emit-ref-confidence GVCF \
        --native-pair-hmm-threads 8 \
        --intervals "$REGION"

    echo "[$(date)] Finished $sample for $REGION → $OUT"
done

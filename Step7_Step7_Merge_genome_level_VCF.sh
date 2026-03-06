#!/bin/bash
#SBATCH --job-name=concat_allchr
#SBATCH --output=04_concat_%j.log
#SBATCH --error=04_concat_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail

module load BCFtools/1.18-GCC-12.3.0 || true
module load HTSlib/1.18-GCC-12.3.0   || true

# =========================
INDIR="../../working/vcfs/"
OUT_VCF="../../working/vcfs/FINAL_genome_level.vcf.gz"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
# =========================

echo "============================================================"
echo "[INFO] Concatenating chromosomes"
echo "[INFO] Input dir : $INDIR"
echo "[INFO] Output    : $OUT_VCF"
echo "============================================================"

vcfs=()
for i in $(seq -w 1 20); do
  f="${INDIR}/Gm${i}_snps_PASS.vcf.gz"
  [[ -f "$f" ]] || { echo "[ERROR] Missing $f"; exit 1; }
  vcfs+=("$f")
done

bcftools concat --threads "$THREADS" -Oz -o "$OUT_VCF" "${vcfs[@]}"
tabix -f -p vcf "$OUT_VCF"

echo "[DONE] $(bcftools index -n "$OUT_VCF") variants"
echo "============================================================"

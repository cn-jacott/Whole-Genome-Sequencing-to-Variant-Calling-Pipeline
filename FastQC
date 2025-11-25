#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --output=fastqc.log
#SBATCH --error=fastqc.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=96:00:00

module load parallel/20230722-GCCcore-12.2.0
module load FastQC/0.12.1-Java-11

DATA_DIR="../working/data"
FASTQC_DIR="../working/fastqc"

mkdir -p "$FASTQC_DIR"

# 1️⃣ Run FastQC on all FASTQs
find "$DATA_DIR" -type f -name "*.fq.gz" | \
    parallel -j 16 fastqc {} -o "$FASTQC_DIR"

# 2️⃣ Rename outputs: prefix with folder name
#    e.g. data/0/E250...fq.gz → 0_E250..._fastqc.html
echo "Renaming FastQC outputs..."
while IFS= read -r fq; do
    sample=$(basename "$fq" .fq.gz)             # E250090877_L01_1_1
    folder=$(basename "$(dirname "$fq")")       # 0, 2, 4, etc.

    for ext in html zip; do
        src="${FASTQC_DIR}/${sample}_fastqc.${ext}"
        dest="${FASTQC_DIR}/${folder}_${sample}_fastqc.${ext}"

        if [[ -f "$src" ]]; then
            mv "$src" "$dest"
        fi
    done
done < <(find "$DATA_DIR" -type f -name "*.fq.gz")

echo "Renaming complete!"

#!/bin/bash
set -euo pipefail

# ============================================================
# CHANGE ONLY THIS
GM="Gm01"
CHUNK=10000000
MAX_CONCURRENT=6
# ============================================================

REFERENCE="../ref/glyma.Wm82.gnm6.S97D.genome_main.fna"
DBDIR="../working/genomicsDB/genomicsdb_${GM}"
INTERVALS_FILE="${DBDIR}/${GM}_${CHUNK}.intervals.list"
OUTDIR="../working/genotype/genotyped_chunks_${GM}_${CHUNK}"
FINAL_VCF="../working/genotype/${GM}_genotyped_chunks.vcf.gz"

mkdir -p "$OUTDIR"

if [[ ! -s "$INTERVALS_FILE" ]]; then
  echo "[ERROR] Missing intervals file: $INTERVALS_FILE"
  echo "        Expected something like: ${DBDIR}/${GM}_${CHUNK}.intervals.list"
  echo "        (This must match the chunks you imported into GenomicsDB.)"
  exit 1
fi

N=$(wc -l < "$INTERVALS_FILE")

echo "============================================================"
echo "Chromosome : $GM"
echo "Chunk size : $CHUNK"
echo "Intervals  : $INTERVALS_FILE  (N=$N)"
echo "DB dir     : $DBDIR"
echo "Outdir     : $OUTDIR"
echo "Final VCF  : $FINAL_VCF"
echo "============================================================"

# Build SLURM array spec safely
ARRAY_SPEC="1-${N}"
if [[ -n "${MAX_CONCURRENT}" ]]; then
  ARRAY_SPEC="1-${N}%${MAX_CONCURRENT}"
fi

# -----------------------
# Submit Genotype array
# -----------------------
ARRAY_JOBID=$(sbatch --array="${ARRAY_SPEC}" \
  --export=ALL,GM="$GM",CHUNK="$CHUNK",REFERENCE="$REFERENCE",DBDIR="$DBDIR",OUTDIR="$OUTDIR",INTERVALS_FILE="$INTERVALS_FILE" \
  <<'EOF' | awk '{print $4}'
#!/bin/bash
#SBATCH --job-name=genotype
#SBATCH --output=genotype_%A_%a.log
#SBATCH --error=genotype_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=96:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cjacott@us.es

set -euo pipefail
module load GATK/4.5.0.0-GCCcore-12.3.0-Java-17

interval="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INTERVALS_FILE")"
if [[ -z "$interval" ]]; then
  echo "[ERROR] No interval for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

# Interval: glyma.Wm82.gnm6.Gm01:1-10000000
# Workspace tag should match how you created them: gendb_Gm01_1_10000000
TAG=$(echo "$interval" | sed 's/.*://; s/-/_/g')
WORKSPACE="${DBDIR}/gendb_${GM}_${TAG}"

# IMPORTANT: Output is relabeled to simple chunk numbers
OUTVCF="${OUTDIR}/${GM}_chunk_${SLURM_ARRAY_TASK_ID}.vcf.gz"

echo "============================================================"
echo "[$(date)] Task     : ${SLURM_ARRAY_TASK_ID}"
echo "Interval           : $interval"
echo "Workspace expected : $WORKSPACE"
echo "Output VCF         : $OUTVCF"
echo "============================================================"

if [[ ! -d "$WORKSPACE" ]]; then
  echo "[ERROR] Missing GenomicsDB workspace: $WORKSPACE"
  echo "[HINT] Check: ls -d ${DBDIR}/gendb_${GM}_* | head"
  exit 1
fi

gatk --java-options "-Xmx28g -XX:+UseParallelGC" GenotypeGVCFs \
  -R "$REFERENCE" \
  -V "gendb://${WORKSPACE}" \
  -L "$interval" \
  -O "$OUTVCF"

echo "[$(date)] Validating gzip integrity: $OUTVCF"
gzip -t "$OUTVCF" 2>/dev/null || { echo "ERROR: gzip test failed for $OUTVCF" >&2; exit 1; }

echo "[$(date)] Indexing: $OUTVCF"
gatk IndexFeatureFile -I "$OUTVCF"

echo "[$(date)] Done chunk ${SLURM_ARRAY_TASK_ID}"
EOF
)

echo "Submitted genotype array job: $ARRAY_JOBID"

# -----------------------
# Submit gather job
# -----------------------
GATHER_JOBID=$(sbatch --dependency=afterok:${ARRAY_JOBID} \
  --export=ALL,GM="$GM",CHUNK="$CHUNK",OUTDIR="$OUTDIR",FINAL_VCF="$FINAL_VCF",INTERVALS_FILE="$INTERVALS_FILE" \
  <<'EOF' | awk '{print $4}'
#!/bin/bash
#SBATCH --job-name=concat
#SBATCH --output=concat_%j.log
#SBATCH --error=concat_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL

set -euo pipefail
module load BCFtools/1.18-GCC-12.3.0

echo "========================================"
echo "[$(date)] CONCAT + CHECK"
echo "Chromosome: $GM"
echo "Chunk size : $CHUNK"
echo "Chunk dir  : $OUTDIR"
echo "Final VCF  : $FINAL_VCF"
echo "========================================"
echo

# Numeric order by chunk ID: chunk_1, chunk_2, ...
mapfile -t files < <(
  ls -1 "${OUTDIR}/${GM}_chunk_"*.vcf.gz 2>/dev/null \
  | sort -t_ -k3,3n
)

echo "Chunks found: ${#files[@]}"
if [[ ${#files[@]} -eq 0 ]]; then
  echo "[ERROR] No chunk VCFs found in $OUTDIR"
  exit 1
fi

echo
echo "Chunks in concat order:"
printf "  %s\n" "${files[@]}"
echo

echo "[$(date)] Concatenating with bcftools concat..."
bcftools concat -a -Oz -o "${FINAL_VCF}" "${files[@]}"

echo "[$(date)] Indexing final VCF..."
bcftools index -t "${FINAL_VCF}"

echo
echo "========================================"
echo "Variant counts per chunk"
echo "========================================"

chunk_sum=0
for f in "${files[@]}"; do
  n=$(bcftools view -H "$f" | wc -l)
  printf "%-35s %10d\n" "$(basename "$f")" "$n"
  chunk_sum=$((chunk_sum + n))
done

final_count=$(bcftools view -H "${FINAL_VCF}" | wc -l)
nsamples=$(bcftools query -l "${FINAL_VCF}" | wc -l)

echo
echo "========================================"
echo "Totals"
echo "========================================"
echo "Sum of chunks : $chunk_sum"
echo "Final VCF     : $final_count"
echo "Samples       : $nsamples"
echo "----------------------------------------"

if [[ "$chunk_sum" -eq "$final_count" ]]; then
  echo "✅ MATCH: chunk sum equals final VCF"
else
  echo "❌ MISMATCH: chunk sum does NOT equal final VCF"
fi

echo
echo "[$(date)] Done: ${FINAL_VCF}"
EOF
)

echo "Submitted concat job: $GATHER_JOBID"
echo "Final VCF will be: $FINAL_VCF"

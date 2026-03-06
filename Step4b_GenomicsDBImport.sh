#!/bin/bash
#SBATCH --job-name=gendb_Gm01
#SBATCH --output=gendb_Gm01_%A_%a.log
#SBATCH --error=gendb_Gm01_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=96:00:00
#SBATCH --array=1-7
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cjacott@us.es

set -euo pipefail
module load GATK/4.5.0.0-GCCcore-12.3.0-Java-17

# ============================================================
# EDIT THESE
GM="Gm01"
CHUNK=10000000
THREADS=4
JAVA_MEM_GB=28
# ============================================================

# Paths
IN_DIR="../working/genomicsDB${GM}"
OUT_DIR="../working/genomicsDB/genomicsdb_${GM}"
REFERENCE="../ref/glyma.Wm82.gnm6.S97D.genome_main.fna"
CONTIG="glyma.Wm82.gnm6.${GM}"

SAMPLE_MAP="${OUT_DIR}/sample_map_${GM}.txt"
INTERVALS="${OUT_DIR}/${GM}_${CHUNK}.intervals.list"

TMP="${SLURM_TMPDIR:-/tmp}"

echo "============================================================"
echo "[INFO] GM        : $GM"
echo "[INFO] CONTIG    : $CONTIG"
echo "[INFO] IN_DIR    : $IN_DIR"
echo "[INFO] OUT_DIR   : $OUT_DIR"
echo "[INFO] CHUNK     : $CHUNK"
echo "[INFO] THREADS   : $THREADS"
echo "[INFO] SLURM_JOB : ${SLURM_JOB_ID:-NA}  TASK: ${SLURM_ARRAY_TASK_ID:-NA}"
echo "============================================================"

# Require prep outputs
[[ -s "$SAMPLE_MAP" ]] || { echo "[ERROR] Missing/empty sample map: $SAMPLE_MAP (run prep script)"; exit 1; }
[[ -s "$INTERVALS" ]]  || { echo "[ERROR] Missing/empty intervals list: $INTERVALS (run prep script)"; exit 1; }

echo "[INFO] sample_map samples: $(wc -l < "$SAMPLE_MAP")"
echo "[INFO] intervals n      : $(wc -l < "$INTERVALS")"

# Pick this task's interval
INTERVAL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INTERVALS" || true)
if [[ -z "${INTERVAL}" ]]; then
  echo "[INFO] No interval for task ${SLURM_ARRAY_TASK_ID}; exiting."
  exit 0
fi

TAG=$(echo "$INTERVAL" | sed 's/.*://; s/-/_/g')
WORKSPACE="${OUT_DIR}/gendb_${GM}_${TAG}"

echo "------------------------------------------------------------"
echo "[INFO] Interval  : $INTERVAL"
echo "[INFO] Workspace : $WORKSPACE"
echo "------------------------------------------------------------"

rm -rf "$WORKSPACE"

gatk --java-options "-Xmx${JAVA_MEM_GB}g -Djava.io.tmpdir=${TMP} -XX:+UseParallelGC" GenomicsDBImport \
  --genomicsdb-workspace-path "$WORKSPACE" \
  --sample-name-map "$SAMPLE_MAP" \
  -L "$INTERVAL" \
  --reader-threads "$THREADS" \
  --tmp-dir "$TMP"

# Truth: report imported callsets from workspace output
if command -v jq >/dev/null 2>&1 && [[ -f "${WORKSPACE}/callset.json" ]]; then
  echo "[INFO] Imported callsets (workspace): $(jq '.callsets | length' "${WORKSPACE}/callset.json")"
fi

echo "[DONE] GenomicsDBImport for $INTERVAL"

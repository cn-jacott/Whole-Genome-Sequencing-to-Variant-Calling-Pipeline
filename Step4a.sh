#!/bin/bash
#SBATCH --job-name=prep_gendb_Gm18
#SBATCH --output=prep_gendb_Gm18_%j.log
#SBATCH --error=prep_gendb_Gm18_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cjacott@us.es

set -euo pipefail

# ============================================================
# EDIT THESE
GM="Gm01"
CHUNK=10000000
# ============================================================

# Paths
IN_DIR="../working/genomicsDB/${GM}"
OUT_DIR="../working/genomicsDB/genomicsdb_${GM}"
REFERENCE="../ref/glyma.Wm82.gnm6.S97D.genome_main.fna"
CONTIG="glyma.Wm82.gnm6.${GM}"

SAMPLE_MAP="${OUT_DIR}/sample_map_${GM}.txt"
INTERVALS="${OUT_DIR}/${GM}_${CHUNK}.intervals.list"

mkdir -p "$OUT_DIR"

echo "============================================================"
echo "[INFO] PREP INPUTS"
echo "[INFO] GM        : $GM"
echo "[INFO] CONTIG    : $CONTIG"
echo "[INFO] IN_DIR    : $IN_DIR"
echo "[INFO] OUT_DIR   : $OUT_DIR"
echo "[INFO] CHUNK     : $CHUNK"
echo "============================================================"

# -----------------------
# Build sample map (atomic write)
# -----------------------
echo "[INFO] Building sample map: $SAMPLE_MAP"
[[ -d "$IN_DIR" ]] || { echo "[ERROR] Input directory not found: $IN_DIR"; exit 1; }

tmp_map="${SAMPLE_MAP}.tmp.$$"
: > "$tmp_map"

shopt -s nullglob
files=("${IN_DIR}"/*.g.vcf.gz)
shopt -u nullglob
[[ ${#files[@]} -gt 0 ]] || { echo "[ERROR] No .g.vcf.gz files in $IN_DIR"; exit 1; }

# Optional: sort for reproducibility
IFS=$'\n' files=($(printf "%s\n" "${files[@]}" | sort))
unset IFS

for f in "${files[@]}"; do
  [[ -f "${f}.tbi" ]] || { echo "[ERROR] Missing index for $f (expected ${f}.tbi)"; exit 1; }
  bn=$(basename "$f")
  sample="${bn%.g.vcf.gz}"
  printf "%s\t%s\n" "$sample" "$f" >> "$tmp_map"
done

mv -f "$tmp_map" "$SAMPLE_MAP"
echo "[INFO] Samples in sample_map: $(wc -l < "$SAMPLE_MAP")"

# Quick sanity: duplicates?
dups=$(cut -f1 "$SAMPLE_MAP" | sort | uniq -d | head -n 5 || true)
if [[ -n "$dups" ]]; then
  echo "[ERROR] Duplicate sample names found in sample_map (showing up to 5):"
  echo "$dups"
  exit 1
fi

# -----------------------
# Build intervals list (atomic write)
# -----------------------
echo "[INFO] Building intervals list: $INTERVALS"
[[ -f "${REFERENCE}.fai" ]] || { echo "[ERROR] Missing ${REFERENCE}.fai"; exit 1; }

LEN=$(awk -v c="$CONTIG" '$1==c{print $2}' "${REFERENCE}.fai" || true)
[[ -n "$LEN" ]] || { echo "[ERROR] Contig not found in .fai: $CONTIG"; exit 1; }

tmp_int="${INTERVALS}.tmp.$$"
python - <<PY > "$tmp_int"
chr="${CONTIG}"
L=int("${LEN}")
step=int("${CHUNK}")
for s in range(1, L+1, step):
    e=min(L, s+step-1)
    print(f"{chr}:{s}-{e}")
PY

mv -f "$tmp_int" "$INTERVALS"
echo "[INFO] Intervals in list: $(wc -l < "$INTERVALS")"

echo "[DONE] Prep complete."
echo "[DONE] sample_map: $SAMPLE_MAP"
echo "[DONE] intervals : $INTERVALS"

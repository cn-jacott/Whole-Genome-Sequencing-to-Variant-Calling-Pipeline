#!/bin/bash
#SBATCH --job-name=filter_variants_Gm01
#SBATCH --output=00_filter_variants_Gm01_%j.log
#SBATCH --error=00_filter_variants_Gm01_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cjacott@us.es

set -euo pipefail

module load GATK/4.5.0.0-GCCcore-12.3.0-Java-17
module load BCFtools/1.18-GCC-12.3.0
module load HTSlib/1.18-GCC-12.3.0 2>/dev/null || true

# ============================================================
# EDIT ONLY THIS
# ============================================================

GM="Gm01"
VARIANT_TYPE="snps"   # options: snps or indels

REFERENCE="../../ref/glyma.Wm82.gnm6.S97D.genome_main.fna"
IN_VCF="../../working/vcfs/gmax/${GM}_genotyped_chunks.vcf.gz"

# Variant filtering thresholds
SNP_FILTER_EXPR="QD < 26.0 || FS > 60.0 || MQ < 40.0"
INDEL_FILTER_EXPR="QD < 26.0 || FS > 200.0 || MQ < 40.0"

OUTDIR="../../working/vcfs/basic_gatk_filtered"

# ============================================================

mkdir -p "$OUTDIR"

# Select filtering settings
if [[ "$VARIANT_TYPE" == "snps" ]]; then
    FILTER_EXPR="$SNP_FILTER_EXPR"
    FILTER_NAME="SNP_filter"
    FINAL_PASS="${OUTDIR}/${GM}_GATK_snps_PASS.vcf.gz"
    LABEL="SNPs"
elif [[ "$VARIANT_TYPE" == "indels" ]]; then
    FILTER_EXPR="$INDEL_FILTER_EXPR"
    FILTER_NAME="INDEL_filter"
    FINAL_PASS="${OUTDIR}/${GM}_GATK_indels_PASS.vcf.gz"
    LABEL="INDELs"
else
    echo "[ERROR] VARIANT_TYPE must be 'snps' or 'indels'"
    exit 1
fi

TMPDIR=$(mktemp -d "${OUTDIR}/GW_tmp.XXXXXX")
trap 'rm -rf "$TMPDIR"' EXIT

TMP_VARIANTS="${TMPDIR}/GW_${VARIANT_TYPE}_tmp.vcf.gz"
VARIANTS_LABELED="${TMPDIR}/GW_${VARIANT_TYPE}_labeled.vcf.gz"

echo "============================================================"
echo "[INFO] Chromosome-level ${LABEL} filtering"
echo "[INFO] Variant type : $VARIANT_TYPE"
echo "[INFO] Filter expr  : $FILTER_EXPR"
echo "[INFO] Input VCF    : $IN_VCF"
echo "============================================================"

[[ -f "$IN_VCF" ]] || { echo "[ERROR] Missing input VCF: $IN_VCF"; exit 1; }

# Index if needed
if [[ ! -f "${IN_VCF}.tbi" && ! -f "${IN_VCF}.csi" ]]; then
    tabix -p vcf "$IN_VCF"
fi

count_records () { bcftools view -H "$1" | wc -l; }

echo "=== COUNTS BEFORE ==="
TOTAL_TYPE=$(bcftools view -H -v "$VARIANT_TYPE" "$IN_VCF" | wc -l)
echo "Total ${LABEL}: $TOTAL_TYPE"
echo "------------------------------------------------------------"

echo "=== EXTRACT ${LABEL} ==="
bcftools view -v "$VARIANT_TYPE" -Oz -o "$TMP_VARIANTS" "$IN_VCF"
gatk IndexFeatureFile -I "$TMP_VARIANTS"

echo "=== GATK VariantFiltration ==="
gatk VariantFiltration \
    -R "$REFERENCE" \
    -V "$TMP_VARIANTS" \
    --filter-name "$FILTER_NAME" \
    --filter-expression "$FILTER_EXPR" \
    -O "$VARIANTS_LABELED"

gatk IndexFeatureFile -I "$VARIANTS_LABELED"

echo "=== KEEP PASS ONLY ==="
bcftools view -f PASS -Oz -o "$FINAL_PASS" "$VARIANTS_LABELED"
tabix -f -p vcf "$FINAL_PASS"

PASS_COUNT=$(count_records "$FINAL_PASS")

echo "============================================================"
echo "${LABEL} PASS : $PASS_COUNT"
echo "Output VCF : $FINAL_PASS"
echo "============================================================"

#!/bin/bash
# Downsample FASTQ files to match read counts for fair comparison
# This ensures both samples have the same number of reads before mapping

set -euo pipefail

# Configuration
THREADS=16
DATA_DIR="/mnt/ssd-8tb/atrx_china/data"
OUTPUT_DIR="/mnt/ssd-8tb/atrx_china/tair12_indel_comparison/data/downsampled"

# Input FASTQ files
COL_FQ="${DATA_DIR}/Col-9day-cotyledon-ul-long.fq.gz"
ATXR56_FQ="${DATA_DIR}/atxr56-9day-cotyledon-ul-long.fq.gz"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "=========================================="
echo "FASTQ Downsampling for Indel Comparison"
echo "=========================================="

# Get read counts (count @ lines / 4 for FASTQ)
echo ""
echo "Counting reads in FASTQ files..."
echo "(This may take a few minutes for large files)"

COL_READS=$(zcat "${COL_FQ}" | awk 'END{print NR/4}')
ATXR56_READS=$(zcat "${ATXR56_FQ}" | awk 'END{print NR/4}')

echo "  Col_9day:    ${COL_READS} reads"
echo "  atxr56_9day: ${ATXR56_READS} reads"

# Determine target (use the smaller count)
if [ "${COL_READS}" -gt "${ATXR56_READS}" ]; then
    TARGET_READS=${ATXR56_READS}
    DOWNSAMPLE_FQ="${COL_FQ}"
    DOWNSAMPLE_NAME="Col_9day"
    KEEP_FQ="${ATXR56_FQ}"
    KEEP_NAME="atxr56_9day"
else
    TARGET_READS=${COL_READS}
    DOWNSAMPLE_FQ="${ATXR56_FQ}"
    DOWNSAMPLE_NAME="atxr56_9day"
    KEEP_FQ="${COL_FQ}"
    KEEP_NAME="Col_9day"
fi

echo ""
echo "Downsampling ${DOWNSAMPLE_NAME} to ${TARGET_READS} reads"
echo "to match ${KEEP_NAME}"

# Use seqtk for downsampling (fast and reliable)
echo ""
echo "Running seqtk sample..."
SEED=42

seqtk sample -s${SEED} "${DOWNSAMPLE_FQ}" ${TARGET_READS} | pigz -p ${THREADS} > "${OUTPUT_DIR}/${DOWNSAMPLE_NAME}.downsampled.fq.gz"

echo "Downsampling complete: ${OUTPUT_DIR}/${DOWNSAMPLE_NAME}.downsampled.fq.gz"

# Create symlink for the file that doesn't need downsampling
echo "Creating symlink for ${KEEP_NAME} (no downsampling needed)..."
ln -sf "${KEEP_FQ}" "${OUTPUT_DIR}/${KEEP_NAME}.downsampled.fq.gz"

# Verify counts
echo ""
echo "Verifying read counts..."
for SAMPLE in Col_9day atxr56_9day; do
    FQ="${OUTPUT_DIR}/${SAMPLE}.downsampled.fq.gz"
    COUNT=$(zcat "${FQ}" | awk 'END{print NR/4}')
    echo "  ${SAMPLE}: ${COUNT} reads"
done

echo ""
echo "=========================================="
echo "Downsampling Complete!"
echo "=========================================="
echo ""
echo "Output files in: ${OUTPUT_DIR}"
echo ""
echo "Next step: Run mapping on the downsampled files"
echo "  bash scripts/01-map_reads_winnowmap.sh ${OUTPUT_DIR}/Col_9day.downsampled.fq.gz results/mapping_downsampled/Col_9day"
echo "  bash scripts/01-map_reads_winnowmap.sh ${OUTPUT_DIR}/atxr56_9day.downsampled.fq.gz results/mapping_downsampled/atxr56_9day"

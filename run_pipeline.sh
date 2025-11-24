#!/bin/bash

# SV Satellite Analyzer - Pipeline Runner
# This script coordinates the full analysis pipeline

set -euo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================

# Default config file
CONFIG="config/config_template.yaml"

# Help message
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Run the complete SV Satellite Analyzer pipeline

OPTIONS:
    -c, --config FILE       Configuration file (default: config/config_template.yaml)
    -b, --bam FILE          Input BAM file (required)
    -s, --sample NAME       Sample name (required)
    -o, --output DIR        Output directory (default: results)
    -m, --min-sv INT        Minimum SV size in bp (default: 50)
    --skip-detection        Skip SV detection (use existing catalogs)
    --skip-visualization    Skip visualization step
    -h, --help              Show this help message

EXAMPLES:
    # Run full pipeline
    $0 -b sample.bam -s MySample -o results -c config/myconfig.yaml

    # Just detect SVs
    $0 -b sample.bam -s MySample --skip-visualization

    # Just visualize existing results
    $0 -s MySample -o results --skip-detection

EOF
    exit 1
}

# Parse arguments
BAM=""
SAMPLE=""
OUTPUT="results"
MIN_SV=50
SKIP_DETECTION=false
SKIP_VIZ=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--config)
            CONFIG="$2"
            shift 2
            ;;
        -b|--bam)
            BAM="$2"
            shift 2
            ;;
        -s|--sample)
            SAMPLE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT="$2"
            shift 2
            ;;
        -m|--min-sv)
            MIN_SV="$2"
            shift 2
            ;;
        --skip-detection)
            SKIP_DETECTION=true
            shift
            ;;
        --skip-visualization)
            SKIP_VIZ=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Validate required arguments
if [[ -z "$SAMPLE" ]]; then
    echo "Error: Sample name (-s) is required"
    usage
fi

if [[ "$SKIP_DETECTION" = false && -z "$BAM" ]]; then
    echo "Error: BAM file (-b) is required unless --skip-detection is set"
    usage
fi

if [[ ! -f "$CONFIG" ]]; then
    echo "Error: Config file not found: $CONFIG"
    exit 1
fi

# ============================================================================
# SETUP
# ============================================================================

echo "============================================================"
echo "SV Satellite Analyzer Pipeline"
echo "============================================================"
echo "Sample: $SAMPLE"
echo "Config: $CONFIG"
echo "Output: $OUTPUT"
echo "============================================================"

# Create output directories
mkdir -p "$OUTPUT/sv_molecules"
mkdir -p "$OUTPUT/statistics"
mkdir -p "$OUTPUT/plots"

# ============================================================================
# STEP 1: SV DETECTION
# ============================================================================

if [[ "$SKIP_DETECTION" = false ]]; then
    echo ""
    echo "[Step 1/2] Detecting SVs in $SAMPLE..."
    echo "--------------------------------------------------------"

    if [[ ! -f "$BAM" ]]; then
        echo "Error: BAM file not found: $BAM"
        exit 1
    fi

    if [[ ! -f "${BAM}.bai" ]]; then
        echo "BAM index not found. Indexing..."
        samtools index "$BAM"
    fi

    python scripts/sv_detection/detect_sv_molecules.py \
        "$BAM" \
        "$OUTPUT/sv_molecules" \
        "$SAMPLE" \
        "$MIN_SV"

    echo "✓ SV detection complete"
    echo ""
else
    echo "[Step 1/2] Skipping SV detection (using existing catalogs)"
fi

# ============================================================================
# STEP 2: VISUALIZATION
# ============================================================================

if [[ "$SKIP_VIZ" = false ]]; then
    echo ""
    echo "[Step 2/2] Generating visualizations..."
    echo "--------------------------------------------------------"

    CATALOG="$OUTPUT/sv_molecules/${SAMPLE}_sv_catalog.tsv"

    if [[ ! -f "$CATALOG" ]]; then
        echo "Error: SV catalog not found: $CATALOG"
        echo "Run without --skip-detection first"
        exit 1
    fi

    # Generate genome-wide ideogram
    echo "Creating genome-wide ideogram..."
    python scripts/visualization/plot_genome_wide_ideogram.py \
        --sv-catalogs "$CATALOG" \
        --sample-names "$SAMPLE" \
        --config "$CONFIG" \
        --output "$OUTPUT/plots/${SAMPLE}_genome_wide_ideogram.png"

    echo "✓ Visualization complete"
    echo ""
else
    echo "[Step 2/2] Skipping visualization"
fi

# ============================================================================
# SUMMARY
# ============================================================================

echo ""
echo "============================================================"
echo "Pipeline Complete!"
echo "============================================================"
echo ""
echo "Output files:"
echo "  SV Catalog:  $OUTPUT/sv_molecules/${SAMPLE}_sv_catalog.tsv"
echo "  Molecules:   $OUTPUT/sv_molecules/${SAMPLE}_sv_molecules.tsv"
echo "  Statistics:  $OUTPUT/sv_molecules/${SAMPLE}_sv_molecule_stats.tsv"

if [[ "$SKIP_VIZ" = false ]]; then
    echo "  Ideogram:    $OUTPUT/plots/${SAMPLE}_genome_wide_ideogram.png"
fi

echo ""
echo "Next steps:"
echo "  - Review SV catalog: less $OUTPUT/sv_molecules/${SAMPLE}_sv_catalog.tsv"
echo "  - Check statistics: less $OUTPUT/sv_molecules/${SAMPLE}_sv_molecule_stats.tsv"

if [[ "$SKIP_VIZ" = false ]]; then
    echo "  - View plots: open $OUTPUT/plots/${SAMPLE}_genome_wide_ideogram.png"
fi

echo ""
echo "For multiple samples:"
echo "  - Run this script for each sample"
echo "  - Use plot_genome_wide_ideogram.py with multiple --sv-catalogs for comparison"
echo ""
echo "============================================================"

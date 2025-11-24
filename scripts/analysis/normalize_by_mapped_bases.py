#!/usr/bin/env python3
"""
Normalize indel rates by MAPPED BASES instead of read count.
This accounts for differences in read length (N50) between samples.

Key insight: atxr5/6 N50 ~61.5kb, Col-0 N50 ~44.7kb
So normalizing by read count underestimates rates in Col-0.

Usage: python 09-normalize_by_mapped_bases.py <sample1_bam> <sample1_indels_dir> <sample1_name> <sample2_bam> <sample2_indels_dir> <sample2_name> <output_dir>
"""

import os
import sys
import pysam
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import defaultdict

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300

# Centromere coordinates for TAIR12 (Chr1-5)
CENTROMERE_REGIONS = {
    'Chr1': (14853761, 17129892),
    'Chr2': (9388656, 11655818),
    'Chr3': (13593457, 15740665),
    'Chr4': (8384369, 11037620),
    'Chr5': (12371485, 15175824),
}

PERICENTROMERE_EXTENSION = 500000

# rDNA regions (will be loaded from BED files)
RDNA_5S_REGIONS = []
RDNA_45S_REGIONS = []


def load_bed_regions(bed_file):
    """Load regions from BED file."""
    regions = []
    if not os.path.exists(bed_file):
        return regions

    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                regions.append((chrom, start, end))

    return regions


def load_rdna_regions():
    """Load rDNA regions from BED files."""
    global RDNA_5S_REGIONS, RDNA_45S_REGIONS

    anno_dir = '/mnt/ssd-8tb/atrx_china/TAIR12/curated_anno'
    bed_5s = os.path.join(anno_dir, '5s_rdna_regions.bed')
    bed_45s = os.path.join(anno_dir, '45s_rdna_regions.bed')

    RDNA_5S_REGIONS = load_bed_regions(bed_5s)
    RDNA_45S_REGIONS = load_bed_regions(bed_45s)

    print(f"  Loaded {len(RDNA_5S_REGIONS)} 5S rDNA regions")
    print(f"  Loaded {len(RDNA_45S_REGIONS)} 45S rDNA regions")


def classify_position(chrom, pos):
    """Classify a genomic position."""
    # Check 5S rDNA regions first
    for region_chrom, start, end in RDNA_5S_REGIONS:
        if chrom == region_chrom and start <= pos < end:
            return '5s_rdna'

    # Check 45S rDNA regions
    for region_chrom, start, end in RDNA_45S_REGIONS:
        if chrom == region_chrom and start <= pos < end:
            return '45s_rdna'

    # Check centromere/pericentromere/arms
    if chrom not in CENTROMERE_REGIONS:
        return 'other'

    cen_start, cen_end = CENTROMERE_REGIONS[chrom]

    if cen_start <= pos <= cen_end:
        return 'centromere'
    elif (cen_start - PERICENTROMERE_EXTENSION <= pos < cen_start or
          cen_end < pos <= cen_end + PERICENTROMERE_EXTENSION):
        return 'pericentromere'
    else:
        return 'arms'


def count_mapped_bases_per_region(bam_file):
    """
    Count MAPPED BASES (not reads) per genomic region.
    Uses the aligned length of each read.
    """
    print(f"  Counting mapped bases in {os.path.basename(bam_file)}...")

    region_bases = defaultdict(int)
    region_reads = defaultdict(int)
    chrom_region_bases = defaultdict(lambda: defaultdict(int))

    total_bases = 0
    total_reads = 0

    try:
        bamfile = pysam.AlignmentFile(bam_file, "rb")
    except Exception as e:
        print(f"Error opening BAM file: {e}")
        return None, None, None, None

    nuclear_chroms = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5']

    for chrom in nuclear_chroms:
        print(f"    Processing {chrom}...")
        for read in bamfile.fetch(chrom):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            # Get aligned length (reference length consumed)
            aligned_length = read.reference_length
            if aligned_length is None or aligned_length == 0:
                continue

            region = classify_position(chrom, read.reference_start)

            region_bases[region] += aligned_length
            region_reads[region] += 1
            chrom_region_bases[chrom][region] += aligned_length

            total_bases += aligned_length
            total_reads += 1

    bamfile.close()

    return dict(region_bases), dict(region_reads), dict(chrom_region_bases), total_bases, total_reads


def load_indel_data(indels_dir, sample_name):
    """Load indel catalog."""
    catalog_file = os.path.join(indels_dir, f"{sample_name}_large_indels_catalog.tsv")

    if not os.path.exists(catalog_file):
        print(f"Error: Catalog file not found: {catalog_file}")
        return None

    df = pd.read_csv(catalog_file, sep='\t')
    nuclear_chroms = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5']
    df = df[df['chromosome'].isin(nuclear_chroms)]

    return df


def calculate_normalized_rates(indels_df, region_bases, per_mb=True):
    """Calculate normalized indel rates per region (per Mb of mapped bases)."""

    rates = {}
    regions = ['5s_rdna', '45s_rdna', 'centromere', 'pericentromere', 'arms']

    for region in regions:
        indel_count = len(indels_df[indels_df['region'] == region])
        mapped_bases = region_bases.get(region, 0)

        if mapped_bases > 0:
            if per_mb:
                rate = (indel_count / mapped_bases) * 1_000_000  # per Mb
            else:
                rate = indel_count / mapped_bases
        else:
            rate = 0

        rates[region] = {
            'indel_count': indel_count,
            'mapped_bases': mapped_bases,
            'mapped_mb': mapped_bases / 1_000_000,
            'rate_per_mb': rate
        }

    return rates


def create_comparison_plots(rates1, rates2, name1, name2, total_bases1, total_bases2, output_dir):
    """Create comparison plots with mapped bases normalization."""

    regions = ['5s_rdna', '45s_rdna', 'centromere', 'pericentromere', 'arms']
    region_labels = ['5S rDNA', '45S rDNA', 'Centromere', 'Pericentromere', 'Arms']

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f'Indel Rates Normalized by Mapped Bases: {name1} vs {name2}',
                fontsize=16, fontweight='bold')

    x = np.arange(len(regions))
    width = 0.35

    # Plot 1: Raw indel counts
    ax = axes[0, 0]
    counts1 = [rates1[r]['indel_count'] for r in regions]
    counts2 = [rates2[r]['indel_count'] for r in regions]

    ax.bar(x - width/2, counts1, width, label=name1, color='#3498db')
    ax.bar(x + width/2, counts2, width, label=name2, color='#e74c3c')

    ax.set_ylabel('Indel Count (raw)', fontsize=12)
    ax.set_title('Raw Indel Counts', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(region_labels, rotation=45, ha='right')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    # Plot 2: Mapped bases per region
    ax = axes[0, 1]
    bases1 = [rates1[r]['mapped_mb'] for r in regions]
    bases2 = [rates2[r]['mapped_mb'] for r in regions]

    ax.bar(x - width/2, bases1, width, label=name1, color='#3498db')
    ax.bar(x + width/2, bases2, width, label=name2, color='#e74c3c')

    ax.set_ylabel('Mapped Bases (Mb)', fontsize=12)
    ax.set_title('Mapped Bases per Region', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(region_labels, rotation=45, ha='right')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    # Plot 3: Normalized rates (per Mb)
    ax = axes[1, 0]
    rates1_norm = [rates1[r]['rate_per_mb'] for r in regions]
    rates2_norm = [rates2[r]['rate_per_mb'] for r in regions]

    ax.bar(x - width/2, rates1_norm, width, label=name1, color='#3498db')
    ax.bar(x + width/2, rates2_norm, width, label=name2, color='#e74c3c')

    ax.set_ylabel('Indels per Mb of mapped bases', fontsize=12)
    ax.set_title('NORMALIZED Indel Rate (per Mb)', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(region_labels, rotation=45, ha='right')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    # Plot 4: Fold change
    ax = axes[1, 1]
    fold_changes = []
    for r in regions:
        if rates1[r]['rate_per_mb'] > 0:
            fc = rates2[r]['rate_per_mb'] / rates1[r]['rate_per_mb']
        else:
            fc = 0
        fold_changes.append(fc)

    colors = ['#27ae60' if fc > 1 else '#e74c3c' for fc in fold_changes]
    ax.bar(x, fold_changes, width=0.6, color=colors, edgecolor='black')
    ax.axhline(y=1, color='black', linestyle='--', linewidth=1)

    ax.set_ylabel(f'Fold Change ({name2}/{name1})', fontsize=12)
    ax.set_title('Fold Change in Indel Rate', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(region_labels, rotation=45, ha='right')
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()

    output_file = os.path.join(output_dir, 'indel_rates_normalized_by_mapped_bases.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved comparison plot: {output_file}")
    plt.close()


def print_summary(rates1, rates2, name1, name2, total_bases1, total_bases2, total_reads1, total_reads2):
    """Print summary statistics."""

    print("\n" + "=" * 80)
    print("INDEL RATES NORMALIZED BY MAPPED BASES")
    print("=" * 80)

    print(f"\nTotal mapped bases:")
    print(f"  {name1}: {total_bases1:,} bp ({total_bases1/1e9:.2f} Gb) from {total_reads1:,} reads")
    print(f"  {name2}: {total_bases2:,} bp ({total_bases2/1e9:.2f} Gb) from {total_reads2:,} reads")
    print(f"  Average read length {name1}: {total_bases1/total_reads1:.0f} bp")
    print(f"  Average read length {name2}: {total_bases2/total_reads2:.0f} bp")

    regions = ['5s_rdna', '45s_rdna', 'centromere', 'pericentromere', 'arms']

    print(f"\n{'Region':<18} {'Sample':<12} {'Indels':<10} {'Mapped Mb':<12} {'Rate/Mb':<12}")
    print("-" * 80)

    for region in regions:
        r1 = rates1[region]
        r2 = rates2[region]

        region_label = region.replace('_', ' ').upper() if 'rdna' in region else region.capitalize()

        print(f"{region_label:<18} {name1:<12} {r1['indel_count']:<10} {r1['mapped_mb']:<12.2f} {r1['rate_per_mb']:<12.2f}")
        print(f"{'':<18} {name2:<12} {r2['indel_count']:<10} {r2['mapped_mb']:<12.2f} {r2['rate_per_mb']:<12.2f}")

        if r1['rate_per_mb'] > 0:
            fold_change = r2['rate_per_mb'] / r1['rate_per_mb']
            print(f"{'':<18} {'Fold change:':<12} {'':<10} {'':<12} {fold_change:<12.2f}x")
        print()

    print("=" * 80)


def main():
    if len(sys.argv) != 8:
        print("Usage: python 09-normalize_by_mapped_bases.py <sample1_bam> <sample1_indels_dir> <sample1_name> <sample2_bam> <sample2_indels_dir> <sample2_name> <output_dir>")
        sys.exit(1)

    sample1_bam = sys.argv[1]
    sample1_indels_dir = sys.argv[2]
    sample1_name = sys.argv[3]
    sample2_bam = sys.argv[4]
    sample2_indels_dir = sys.argv[5]
    sample2_name = sys.argv[6]
    output_dir = sys.argv[7]

    os.makedirs(output_dir, exist_ok=True)

    print("=" * 80)
    print("Indel Rate Analysis - Normalized by MAPPED BASES")
    print("=" * 80)
    print(f"Sample 1: {sample1_name}")
    print(f"Sample 2: {sample2_name}")
    print(f"Output: {output_dir}")

    # Load rDNA regions
    print("\n[Step 0] Loading rDNA regions...")
    load_rdna_regions()

    # Count mapped bases per region
    print("\n[Step 1] Counting mapped bases per region...")
    region_bases1, region_reads1, chrom_bases1, total_bases1, total_reads1 = count_mapped_bases_per_region(sample1_bam)
    region_bases2, region_reads2, chrom_bases2, total_bases2, total_reads2 = count_mapped_bases_per_region(sample2_bam)

    if region_bases1 is None or region_bases2 is None:
        print("Error counting mapped bases!")
        sys.exit(1)

    # Load indel data
    print("\n[Step 2] Loading indel data...")
    indels1 = load_indel_data(sample1_indels_dir, sample1_name)
    indels2 = load_indel_data(sample2_indels_dir, sample2_name)

    if indels1 is None or indels2 is None:
        print("Error loading indels!")
        sys.exit(1)

    print(f"  {sample1_name}: {len(indels1)} indels")
    print(f"  {sample2_name}: {len(indels2)} indels")

    # Calculate normalized rates
    print("\n[Step 3] Calculating normalized rates (per Mb of mapped bases)...")
    rates1 = calculate_normalized_rates(indels1, region_bases1)
    rates2 = calculate_normalized_rates(indels2, region_bases2)

    # Save rates to file
    print("\n[Step 4] Saving results...")

    rates_df = []
    for sample_name, rates, total_bases in [(sample1_name, rates1, total_bases1), (sample2_name, rates2, total_bases2)]:
        for region in ['5s_rdna', '45s_rdna', 'centromere', 'pericentromere', 'arms']:
            rates_df.append({
                'sample': sample_name,
                'region': region,
                'indel_count': rates[region]['indel_count'],
                'mapped_bases': rates[region]['mapped_bases'],
                'mapped_mb': rates[region]['mapped_mb'],
                'indels_per_mb': rates[region]['rate_per_mb'],
                'total_mapped_bases': total_bases
            })

    rates_file = os.path.join(output_dir, 'indel_rates_by_mapped_bases.tsv')
    pd.DataFrame(rates_df).to_csv(rates_file, sep='\t', index=False)
    print(f"  Saved rates: {rates_file}")

    # Create plots
    print("\n[Step 5] Creating plots...")
    create_comparison_plots(rates1, rates2, sample1_name, sample2_name, total_bases1, total_bases2, output_dir)

    # Print summary
    print_summary(rates1, rates2, sample1_name, sample2_name, total_bases1, total_bases2, total_reads1, total_reads2)


if __name__ == "__main__":
    main()

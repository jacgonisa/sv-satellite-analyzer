#!/usr/bin/env python3
"""
Parse 5S and 45S rDNA GFF files and create BED files.
Merges adjacent regions and excludes overlaps with centromeres/pericentromeres.

Usage: python 08-create_rdna_beds.py
"""

import os
import sys

# Chromosome accession to name mapping
CHROM_MAP = {
    'CP116280.1': 'Chr1',
    'CP116281.2': 'Chr2',
    'CP116282.1': 'Chr3',
    'CP116283.2': 'Chr4',
    'CP116284.1': 'Chr5',
}

# Centromere coordinates (1-based, inclusive)
CENTROMERE_REGIONS = {
    'Chr1': (14853761, 17129892),
    'Chr2': (9388656, 11655818),
    'Chr3': (13593457, 15740665),
    'Chr4': (8384369, 11037620),
    'Chr5': (12371485, 15175824),
}

PERICENTROMERE_EXTENSION = 500000


def parse_gff_file(gff_file):
    """
    Parse GFF file that has all entries on one line (no line breaks).
    Returns list of (chrom, start, end) tuples.
    """
    regions = []

    with open(gff_file, 'r') as f:
        content = f.read()

    # Split by the pattern that indicates a new entry
    # Each entry ends with "Composition:XXX" followed by the next chromosome ID
    entries = []
    current_entry = []

    for char in content:
        current_entry.append(char)
        # Check if we hit a tab followed by chromosome ID pattern
        entry_str = ''.join(current_entry)
        if '\t' in entry_str and any(acc in entry_str for acc in CHROM_MAP.keys()):
            # Try to extract a complete entry
            parts = entry_str.split('\t')
            if len(parts) >= 9:
                # We have a complete entry
                chrom_acc = parts[0]
                try:
                    start = int(parts[3])
                    end = int(parts[4])

                    # Convert to standard chromosome name
                    if chrom_acc in CHROM_MAP:
                        chrom = CHROM_MAP[chrom_acc]
                        regions.append((chrom, start, end))
                except:
                    pass

                # Start new entry with remaining characters
                # Find where the next entry starts
                last_tab_idx = entry_str.rfind('\t')
                if last_tab_idx > 0:
                    # Check if there's a chromosome ID after the last tab
                    remaining = entry_str[last_tab_idx+1:]
                    if any(acc in remaining for acc in CHROM_MAP.keys()):
                        current_entry = list(remaining)
                    else:
                        current_entry = []

    print(f"  Parsed {len(regions)} regions from {os.path.basename(gff_file)}")
    return regions


def parse_gff_simple(gff_file):
    """
    Simpler GFF parser - split by accession IDs.
    """
    regions = []

    with open(gff_file, 'r') as f:
        content = f.read()

    # Replace each accession ID with newline + accession to create line breaks
    for acc in CHROM_MAP.keys():
        content = content.replace(acc, '\n' + acc)

    # Now parse line by line
    for line in content.split('\n'):
        line = line.strip()
        if not line:
            continue

        fields = line.split('\t')
        if len(fields) >= 5:
            chrom_acc = fields[0]
            try:
                start = int(fields[3])
                end = int(fields[4])

                if chrom_acc in CHROM_MAP:
                    chrom = CHROM_MAP[chrom_acc]
                    regions.append((chrom, start, end))
            except:
                continue

    print(f"  Parsed {len(regions)} regions from {os.path.basename(gff_file)}")
    return regions


def merge_regions(regions):
    """
    Merge overlapping or adjacent regions.
    Input: list of (chrom, start, end) tuples (1-based, inclusive)
    Output: merged list
    """
    if not regions:
        return []

    # Sort by chromosome and start position
    sorted_regions = sorted(regions, key=lambda x: (x[0], x[1]))

    merged = []
    current_chrom, current_start, current_end = sorted_regions[0]

    for chrom, start, end in sorted_regions[1:]:
        if chrom == current_chrom and start <= current_end + 1:
            # Overlapping or adjacent - merge
            current_end = max(current_end, end)
        else:
            # Not overlapping - save current and start new
            merged.append((current_chrom, current_start, current_end))
            current_chrom, current_start, current_end = chrom, start, end

    # Don't forget the last region
    merged.append((current_chrom, current_start, current_end))

    print(f"  Merged {len(regions)} regions into {len(merged)} blocks")
    return merged


def overlaps_with_centromere_or_peri(chrom, start, end):
    """
    Check if region overlaps with centromere or pericentromere.
    Coordinates are 1-based, inclusive.
    """
    if chrom not in CENTROMERE_REGIONS:
        return False

    cen_start, cen_end = CENTROMERE_REGIONS[chrom]
    peri_start = cen_start - PERICENTROMERE_EXTENSION
    peri_end = cen_end + PERICENTROMERE_EXTENSION

    # Check for overlap
    return not (end < peri_start or start > peri_end)


def exclude_overlapping_regions(regions):
    """
    Exclude regions that overlap with centromeres or pericentromeres.
    """
    filtered = []
    excluded_count = 0

    for chrom, start, end in regions:
        if overlaps_with_centromere_or_peri(chrom, start, end):
            excluded_count += 1
        else:
            filtered.append((chrom, start, end))

    print(f"  Excluded {excluded_count} regions overlapping centromere/pericentromere")
    print(f"  Retained {len(filtered)} non-overlapping regions")
    return filtered


def write_bed_file(regions, output_file, feature_name):
    """
    Write regions to BED file.
    BED format is 0-based, half-open [start, end)
    GFF is 1-based, inclusive [start, end]
    """
    with open(output_file, 'w') as f:
        f.write(f"# {feature_name} rDNA regions (TAIR12)\n")
        f.write(f"# BED format: 0-based, half-open [start, end)\n")
        f.write(f"# chrom\tstart\tend\tname\n")

        for i, (chrom, start, end) in enumerate(regions, 1):
            # Convert from 1-based inclusive to 0-based half-open
            bed_start = start - 1
            bed_end = end
            name = f"{feature_name}_rDNA_{i}"
            f.write(f"{chrom}\t{bed_start}\t{bed_end}\t{name}\n")

    print(f"  Wrote {len(regions)} regions to {output_file}")


def calculate_stats(regions):
    """Calculate statistics about regions."""
    by_chrom = {}
    total_bp = 0

    for chrom, start, end in regions:
        length = end - start + 1  # 1-based inclusive
        total_bp += length

        if chrom not in by_chrom:
            by_chrom[chrom] = {'count': 0, 'bp': 0}
        by_chrom[chrom]['count'] += 1
        by_chrom[chrom]['bp'] += length

    return by_chrom, total_bp


def main():
    # Input GFF files
    gff_5s = "/mnt/ssd-8tb/atrx_china/TAIR12/curated_anno/CC_Col_v2_5s_02022024.gff"
    gff_45s = "/mnt/ssd-8tb/atrx_china/TAIR12/curated_anno/CC_Col_v2_45s_02022024.gff"

    # Output directory
    output_dir = "/mnt/ssd-8tb/atrx_china/TAIR12/curated_anno"

    print("=" * 70)
    print("Processing rDNA GFF files")
    print("=" * 70)

    # Process 5S rDNA
    print("\n[1] Processing 5S rDNA...")
    regions_5s = parse_gff_simple(gff_5s)
    regions_5s = merge_regions(regions_5s)
    regions_5s = exclude_overlapping_regions(regions_5s)

    bed_5s = os.path.join(output_dir, "5s_rdna_regions.bed")
    write_bed_file(regions_5s, bed_5s, "5S")

    # Process 45S rDNA
    print("\n[2] Processing 45S rDNA...")
    regions_45s = parse_gff_simple(gff_45s)
    regions_45s = merge_regions(regions_45s)
    regions_45s = exclude_overlapping_regions(regions_45s)

    bed_45s = os.path.join(output_dir, "45s_rdna_regions.bed")
    write_bed_file(regions_45s, bed_45s, "45S")

    # Generate summary report
    print("\n" + "=" * 70)
    print("SUMMARY REPORT")
    print("=" * 70)

    print("\n5S rDNA regions:")
    by_chrom_5s, total_bp_5s = calculate_stats(regions_5s)
    for chrom in sorted(by_chrom_5s.keys()):
        print(f"  {chrom}: {by_chrom_5s[chrom]['count']} regions, {by_chrom_5s[chrom]['bp']:,} bp")
    print(f"  Total: {len(regions_5s)} regions, {total_bp_5s:,} bp")

    print("\n45S rDNA regions:")
    by_chrom_45s, total_bp_45s = calculate_stats(regions_45s)
    for chrom in sorted(by_chrom_45s.keys()):
        print(f"  {chrom}: {by_chrom_45s[chrom]['count']} regions, {by_chrom_45s[chrom]['bp']:,} bp")
    print(f"  Total: {len(regions_45s)} regions, {total_bp_45s:,} bp")

    print("\n" + "=" * 70)
    print("BED files created:")
    print(f"  {bed_5s}")
    print(f"  {bed_45s}")
    print("=" * 70)


if __name__ == "__main__":
    main()

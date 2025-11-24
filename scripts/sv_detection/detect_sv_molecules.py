#!/usr/bin/env python3
"""
Identify single molecules (reads) with structural variations.
Extracts detailed information about SV position on both genome and read.

Output:
- Per-read summary with all SVs found in each molecule
- Detailed catalog with genomic and read-relative positions

Usage: python 10-identify_sv_molecules.py <bam_file> <output_dir> <sample_name> [min_sv_size]
"""

import os
import sys
import pysam
import pandas as pd
import csv
from collections import defaultdict

# Centromere coordinates for TAIR12 (Chr1-5)
CENTROMERE_REGIONS = {
    'Chr1': (14853761, 17129892),
    'Chr2': (9388656, 11655818),
    'Chr3': (13593457, 15740665),
    'Chr4': (8384369, 11037620),
    'Chr5': (12371485, 15175824),
}

PERICENTROMERE_EXTENSION = 500000

# rDNA regions (loaded from BED files)
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


def extract_sv_molecules(bam_file, sv_catalog_file, molecule_summary_file, min_sv_size=50):
    """
    Extract molecules (reads) with SVs and detailed positional information.

    For each SV, records:
    - Genomic position (ref_pos)
    - Read position (read_pos) - absolute position on the read
    - Relative position (sv_pos_pct) - percentage along the read
    """
    try:
        bamfile = pysam.AlignmentFile(bam_file, "rb")
    except Exception as e:
        print(f"Error opening BAM file: {e}")
        return None

    nuclear_chroms = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5']

    # Per-SV catalog
    sv_fieldnames = [
        'read_id', 'read_length', 'chromosome',
        'read_start', 'read_end',  # Read alignment span on genome
        'sv_ref_pos', 'sv_read_pos', 'sv_pos_pct',  # SV position info
        'sv_type', 'sv_size', 'region',
        'mapping_quality', 'is_178_multiple', 'strand'
    ]

    # Per-molecule summary
    mol_fieldnames = [
        'read_id', 'read_length', 'chromosome',
        'read_start', 'read_end', 'region',
        'mapping_quality', 'strand',
        'n_insertions', 'n_deletions', 'n_total_svs',
        'total_ins_bp', 'total_del_bp',
        'sv_positions_on_read',  # Comma-separated list
        'sv_positions_on_genome',  # Comma-separated list
        'sv_sizes',  # Comma-separated list
        'sv_types'  # Comma-separated list
    ]

    total_reads = 0
    reads_with_sv = 0
    total_svs = 0

    with open(sv_catalog_file, 'w', newline='') as f_sv, \
         open(molecule_summary_file, 'w', newline='') as f_mol:

        sv_writer = csv.DictWriter(f_sv, fieldnames=sv_fieldnames, delimiter='\t')
        sv_writer.writeheader()

        mol_writer = csv.DictWriter(f_mol, fieldnames=mol_fieldnames, delimiter='\t')
        mol_writer.writeheader()

        sv_buffer = []
        mol_buffer = []
        chunk_size = 1000

        for chrom in nuclear_chroms:
            print(f"  Processing {chrom}...")
            chrom_reads = 0
            chrom_sv_reads = 0

            for read in bamfile.fetch(chrom):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue

                total_reads += 1
                chrom_reads += 1

                read_length = read.query_length if read.query_length else 0
                if read_length == 0:
                    continue

                # Track positions on read and genome
                read_pos = 0
                ref_pos = read.reference_start

                # Read alignment span
                read_start = read.reference_start
                read_end = read.reference_end if read.reference_end else read.reference_start

                # Strand
                strand = '-' if read.is_reverse else '+'

                # Region classification at read start
                region = classify_position(chrom, read_start)

                # Collect SVs for this read
                read_svs = []

                if read.cigartuples:
                    for op, length in read.cigartuples:
                        if op == 0:  # M (match/mismatch)
                            read_pos += length
                            ref_pos += length
                        elif op == 1:  # I (insertion)
                            if length >= min_sv_size:
                                sv_pos_pct = (read_pos / read_length) * 100
                                is_178_mult = (length % 178 == 0)
                                sv_region = classify_position(chrom, ref_pos)

                                sv_record = {
                                    'read_id': read.query_name,
                                    'read_length': read_length,
                                    'chromosome': chrom,
                                    'read_start': read_start,
                                    'read_end': read_end,
                                    'sv_ref_pos': ref_pos,
                                    'sv_read_pos': read_pos,
                                    'sv_pos_pct': round(sv_pos_pct, 2),
                                    'sv_type': 'INS',
                                    'sv_size': length,
                                    'region': sv_region,
                                    'mapping_quality': read.mapping_quality,
                                    'is_178_multiple': is_178_mult,
                                    'strand': strand
                                }
                                read_svs.append(sv_record)
                                sv_buffer.append(sv_record)
                                total_svs += 1

                            read_pos += length
                        elif op == 2:  # D (deletion)
                            if length >= min_sv_size:
                                sv_pos_pct = (read_pos / read_length) * 100
                                is_178_mult = (length % 178 == 0)
                                sv_region = classify_position(chrom, ref_pos)

                                sv_record = {
                                    'read_id': read.query_name,
                                    'read_length': read_length,
                                    'chromosome': chrom,
                                    'read_start': read_start,
                                    'read_end': read_end,
                                    'sv_ref_pos': ref_pos,
                                    'sv_read_pos': read_pos,
                                    'sv_pos_pct': round(sv_pos_pct, 2),
                                    'sv_type': 'DEL',
                                    'sv_size': length,
                                    'region': sv_region,
                                    'mapping_quality': read.mapping_quality,
                                    'is_178_multiple': is_178_mult,
                                    'strand': strand
                                }
                                read_svs.append(sv_record)
                                sv_buffer.append(sv_record)
                                total_svs += 1

                            ref_pos += length
                        elif op == 4:  # S (soft clip)
                            read_pos += length
                        elif op == 3:  # N (skip)
                            ref_pos += length

                # Create molecule summary if SVs found
                if read_svs:
                    reads_with_sv += 1
                    chrom_sv_reads += 1

                    n_ins = sum(1 for sv in read_svs if sv['sv_type'] == 'INS')
                    n_del = sum(1 for sv in read_svs if sv['sv_type'] == 'DEL')
                    total_ins_bp = sum(sv['sv_size'] for sv in read_svs if sv['sv_type'] == 'INS')
                    total_del_bp = sum(sv['sv_size'] for sv in read_svs if sv['sv_type'] == 'DEL')

                    mol_record = {
                        'read_id': read.query_name,
                        'read_length': read_length,
                        'chromosome': chrom,
                        'read_start': read_start,
                        'read_end': read_end,
                        'region': region,
                        'mapping_quality': read.mapping_quality,
                        'strand': strand,
                        'n_insertions': n_ins,
                        'n_deletions': n_del,
                        'n_total_svs': len(read_svs),
                        'total_ins_bp': total_ins_bp,
                        'total_del_bp': total_del_bp,
                        'sv_positions_on_read': ','.join(str(sv['sv_read_pos']) for sv in read_svs),
                        'sv_positions_on_genome': ','.join(str(sv['sv_ref_pos']) for sv in read_svs),
                        'sv_sizes': ','.join(str(sv['sv_size']) for sv in read_svs),
                        'sv_types': ','.join(sv['sv_type'] for sv in read_svs)
                    }
                    mol_buffer.append(mol_record)

                # Write buffers periodically
                if len(sv_buffer) >= chunk_size:
                    sv_writer.writerows(sv_buffer)
                    sv_buffer = []
                if len(mol_buffer) >= chunk_size:
                    mol_writer.writerows(mol_buffer)
                    mol_buffer = []

            print(f"    {chrom}: {chrom_reads:,} reads, {chrom_sv_reads:,} with SVs")

        # Write remaining buffers
        if sv_buffer:
            sv_writer.writerows(sv_buffer)
        if mol_buffer:
            mol_writer.writerows(mol_buffer)

    bamfile.close()

    return {
        'total_reads': total_reads,
        'reads_with_sv': reads_with_sv,
        'total_svs': total_svs
    }


def generate_summary_stats(molecule_summary_file, output_dir, sample_name):
    """Generate summary statistics from molecule data."""

    df = pd.read_csv(molecule_summary_file, sep='\t')

    stats = {
        'sample': sample_name,
        'total_molecules_with_sv': len(df),
        'unique_reads': df['read_id'].nunique(),
        'mean_read_length': df['read_length'].mean(),
        'median_read_length': df['read_length'].median(),
        'mean_svs_per_read': df['n_total_svs'].mean(),
        'max_svs_per_read': df['n_total_svs'].max(),
    }

    # By region
    for region in ['5s_rdna', '45s_rdna', 'centromere', 'pericentromere', 'arms']:
        region_df = df[df['region'] == region]
        stats[f'{region}_molecules'] = len(region_df)
        if len(region_df) > 0:
            stats[f'{region}_mean_svs'] = region_df['n_total_svs'].mean()
            stats[f'{region}_mean_read_len'] = region_df['read_length'].mean()
        else:
            stats[f'{region}_mean_svs'] = 0
            stats[f'{region}_mean_read_len'] = 0

    # Save stats
    stats_file = os.path.join(output_dir, f"{sample_name}_sv_molecule_stats.tsv")
    pd.DataFrame([stats]).to_csv(stats_file, sep='\t', index=False)
    print(f"  Saved statistics: {stats_file}")

    return stats


def main():
    if len(sys.argv) < 4:
        print("Usage: python 10-identify_sv_molecules.py <bam_file> <output_dir> <sample_name> [min_sv_size]")
        print("  min_sv_size: minimum SV size to extract (default: 50)")
        sys.exit(1)

    bam_file = sys.argv[1]
    output_dir = sys.argv[2]
    sample_name = sys.argv[3]
    min_sv_size = int(sys.argv[4]) if len(sys.argv) > 4 else 50

    os.makedirs(output_dir, exist_ok=True)

    print("=" * 70)
    print(f"Identifying SV-containing molecules for: {sample_name}")
    print("=" * 70)
    print(f"BAM file: {bam_file}")
    print(f"Output directory: {output_dir}")
    print(f"Minimum SV size: {min_sv_size}bp")

    # Load rDNA regions
    print(f"\n[Step 0] Loading rDNA regions...")
    load_rdna_regions()

    # Output files
    sv_catalog_file = os.path.join(output_dir, f"{sample_name}_sv_catalog.tsv")
    molecule_summary_file = os.path.join(output_dir, f"{sample_name}_sv_molecules.tsv")

    # Extract SVs
    print(f"\n[Step 1] Extracting SV-containing molecules...")
    results = extract_sv_molecules(bam_file, sv_catalog_file, molecule_summary_file, min_sv_size)

    if results is None:
        print("Error extracting SVs!")
        sys.exit(1)

    print(f"\n  Total reads processed: {results['total_reads']:,}")
    print(f"  Reads with SVs (>={min_sv_size}bp): {results['reads_with_sv']:,}")
    print(f"  Total SVs found: {results['total_svs']:,}")
    print(f"\n  Saved SV catalog: {sv_catalog_file}")
    print(f"  Saved molecule summary: {molecule_summary_file}")

    # Generate summary stats
    print(f"\n[Step 2] Generating summary statistics...")
    stats = generate_summary_stats(molecule_summary_file, output_dir, sample_name)

    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Molecules with SVs: {stats['total_molecules_with_sv']:,}")
    print(f"Mean read length: {stats['mean_read_length']:,.0f} bp")
    print(f"Mean SVs per read: {stats['mean_svs_per_read']:.2f}")
    print(f"Max SVs in single read: {stats['max_svs_per_read']}")
    print(f"\nBy region:")
    for region in ['5s_rdna', '45s_rdna', 'centromere', 'pericentromere', 'arms']:
        region_name = region.replace('_', ' ').upper() if 'rdna' in region else region.capitalize()
        print(f"  {region_name}: {stats[f'{region}_molecules']:,} molecules")
    print("=" * 70)


if __name__ == "__main__":
    main()

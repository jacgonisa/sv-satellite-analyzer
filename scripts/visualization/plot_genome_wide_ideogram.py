#!/usr/bin/env python3
"""
Generate genome-wide ideogram visualization of 178bp-multiple SVs.

Creates a comprehensive summary plot showing:
- All chromosomes with genomic features (centromeres, pericentromeres, rDNA)
- SV positions (insertions above, deletions below)
- Count summaries per chromosome

Usage:
    python plot_genome_wide_ideogram.py \\
        --sv-catalogs sample1.tsv sample2.tsv \\
        --sample-names Sample1 Sample2 \\
        --config config.yaml \\
        --output results/plots
"""

import os
import sys
import argparse
import yaml
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np

sns.set_style("white")


def load_config(config_file):
    """Load configuration from YAML file."""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)


def load_bed_regions(bed_file):
    """Load regions from BED file."""
    regions = []
    if not os.path.exists(bed_file):
        print(f"Warning: BED file not found: {bed_file}")
        return regions

    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                regions.append((fields[0], int(fields[1]), int(fields[2])))
    return regions


def load_pericentromeres(bed_file):
    """Load pericentromere coordinates from BED file, grouped by chromosome."""
    peri_dict = {}
    if not os.path.exists(bed_file):
        print(f"Warning: Pericentromere BED file not found: {bed_file}")
        return peri_dict

    # Mapping RefSeq accession IDs to common chromosome names
    chr_map = {
        'CP116280.1': 'Chr1',
        'CP116281.2': 'Chr2',
        'CP116282.1': 'Chr3',
        'CP116283.2': 'Chr4',
        'CP116284.1': 'Chr5',
    }

    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                accession = fields[0]
                # Convert RefSeq ID to common name if needed
                chrom = chr_map.get(accession, accession)
                start = int(fields[1])
                end = int(fields[2])

                if chrom not in peri_dict:
                    peri_dict[chrom] = []
                peri_dict[chrom].append((start, end))

    return peri_dict


def plot_genome_wide_ideogram(sv_catalogs, sample_names, config, output_file):
    """
    Generate genome-wide ideogram plot.

    Parameters:
    -----------
    sv_catalogs : list of str
        Paths to SV catalog TSV files
    sample_names : list of str
        Display names for each sample
    config : dict
        Configuration dictionary
    output_file : str
        Output file path
    """
    # Load SV data
    sv_data = []
    for catalog_file, sample_name in zip(sv_catalogs, sample_names):
        df = pd.read_csv(catalog_file, sep='\t')
        df['sample'] = sample_name
        sv_data.append(df)

    # Extract config
    chroms = config['genome']['chromosomes']
    chrom_lengths = config['genome']['lengths']
    centromeres = config['genome']['centromeres']

    # Load pericentromere regions from BED file if provided
    if 'pericentromere_bed' in config and config['pericentromere_bed']:
        pericentromeres = load_pericentromeres(config['pericentromere_bed'])
    else:
        # Fallback to extension-based calculation if BED not provided
        peri_extension = config.get('pericentromere_extension', 500000)
        pericentromeres = {}
        for chrom in chroms:
            cen_start, cen_end = centromeres[chrom]
            chrom_len = chrom_lengths[chrom]
            pericentromeres[chrom] = [
                (max(0, cen_start - peri_extension), cen_start),
                (cen_end, min(chrom_len, cen_end + peri_extension))
            ]

    # Load rDNA regions
    rdna_5s = load_bed_regions(config['rdna_regions']['5s'])
    rdna_45s = load_bed_regions(config['rdna_regions']['45s'])

    # Create figure
    n_samples = len(sv_data)
    fig = plt.figure(figsize=(20, 5 + n_samples * 2))
    gs = fig.add_gridspec(n_samples + 2, 1,
                          height_ratios=[1] * n_samples + [0.2, 1],
                          hspace=0.4)

    fig.suptitle('Genome-wide Summary: 178bp-Multiple SVs',
                 fontsize=20, fontweight='bold', y=0.98)

    chrom_spacing = 1.5
    max_len = max(chrom_lengths.values())

    # Plot each sample
    for sample_idx, (df, sample_name) in enumerate(zip(sv_data, sample_names)):
        ax = fig.add_subplot(gs[sample_idx, 0])
        ax.set_title(sample_name, fontsize=16, fontweight='bold', loc='left')
        ax.set_xlim(-1, max_len/1e6 + 1)
        ax.set_ylim(-0.5, len(chroms) * chrom_spacing)
        ax.set_ylabel('Chromosome', fontsize=12)
        ax.set_yticks([i * chrom_spacing for i in range(len(chroms))])
        ax.set_yticklabels(chroms)

        if sample_idx == n_samples - 1:
            ax.set_xlabel('Position (Mb)', fontsize=12)
        else:
            ax.set_xlabel('')

        for i, chrom in enumerate(chroms):
            y_pos = i * chrom_spacing
            chrom_len = chrom_lengths[chrom]
            cen_start, cen_end = centromeres[chrom]

            # Draw chromosome backbone
            ax.plot([0, chrom_len/1e6], [y_pos, y_pos], 'k-',
                   linewidth=10, solid_capstyle='round', alpha=0.2)

            # Draw pericentromeres from loaded regions
            if chrom in pericentromeres:
                for peri_start, peri_end in pericentromeres[chrom]:
                    ax.plot([peri_start/1e6, peri_end/1e6], [y_pos, y_pos],
                           color=config['visualization']['colors']['pericentromere'],
                           linewidth=10, solid_capstyle='round', alpha=0.6)

            # Draw centromere
            ax.plot([cen_start/1e6, cen_end/1e6], [y_pos, y_pos],
                   color=config['visualization']['colors']['centromere'],
                   linewidth=10, solid_capstyle='round', alpha=0.8)

            # Draw 5S rDNA
            for rdna_chrom, rdna_start, rdna_end in rdna_5s:
                if rdna_chrom == chrom:
                    ax.plot([rdna_start/1e6, rdna_end/1e6], [y_pos, y_pos],
                           color=config['visualization']['colors']['rdna_5s'],
                           linewidth=10, solid_capstyle='round', alpha=0.9)
                    mid = (rdna_start + rdna_end) / 2 / 1e6
                    ax.text(mid, y_pos + 0.55, '5S', fontsize=8, ha='center',
                           fontweight='bold',
                           color=config['visualization']['colors']['rdna_5s'],
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                   edgecolor=config['visualization']['colors']['rdna_5s'],
                                   linewidth=1.5))

            # Draw 45S rDNA
            for rdna_chrom, rdna_start, rdna_end in rdna_45s:
                if rdna_chrom == chrom:
                    ax.plot([rdna_start/1e6, rdna_end/1e6], [y_pos, y_pos],
                           color=config['visualization']['colors']['rdna_45s'],
                           linewidth=10, solid_capstyle='round', alpha=0.9)
                    mid = (rdna_start + rdna_end) / 2 / 1e6
                    ax.text(mid, y_pos + 0.55, '45S', fontsize=8, ha='center',
                           fontweight='bold',
                           color=config['visualization']['colors']['rdna_45s'],
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                   edgecolor=config['visualization']['colors']['rdna_45s'],
                                   linewidth=1.5))

            # Add centromere label
            cen_mid = (cen_start + cen_end) / 2 / 1e6
            ax.text(cen_mid, y_pos + 0.55, 'CEN', fontsize=8, ha='center',
                   fontweight='bold',
                   color=config['visualization']['colors']['centromere'],
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                           edgecolor=config['visualization']['colors']['centromere'],
                           linewidth=1.5))

            # Add pericentromere labels at actual region midpoints
            if chrom in pericentromeres:
                for peri_start, peri_end in pericentromeres[chrom]:
                    peri_mid = (peri_start + peri_end) / 2 / 1e6
                    ax.text(peri_mid, y_pos - 0.55, 'PERI', fontsize=7, ha='center',
                           fontweight='normal',
                           color=config['visualization']['colors']['pericentromere'],
                           alpha=0.8,
                           bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                                   edgecolor=config['visualization']['colors']['pericentromere'],
                                   linewidth=1))

            # Plot SVs
            sample_chrom = df[(df['chromosome'] == chrom) & (df['is_178_multiple'] == True)]
            sample_ins = sample_chrom[sample_chrom['sv_type'] == 'INS']
            sample_del = sample_chrom[sample_chrom['sv_type'] == 'DEL']

            ax.scatter(sample_ins['sv_ref_pos']/1e6, [y_pos + 0.25] * len(sample_ins),
                      s=8, c=config['visualization']['colors']['insertions'],
                      alpha=0.6, marker='^', zorder=10)
            ax.scatter(sample_del['sv_ref_pos']/1e6, [y_pos - 0.25] * len(sample_del),
                      s=8, c=config['visualization']['colors']['deletions'],
                      alpha=0.6, marker='v', zorder=10)

        # Add legend (only once)
        if sample_idx == n_samples - 1:
            from matplotlib.lines import Line2D
            legend_elements = [
                Line2D([0], [0], marker='^', color='w',
                       markerfacecolor=config['visualization']['colors']['insertions'],
                       markersize=10, label='Insertions'),
                Line2D([0], [0], marker='v', color='w',
                       markerfacecolor=config['visualization']['colors']['deletions'],
                       markersize=10, label='Deletions'),
                Line2D([0], [0], color=config['visualization']['colors']['centromere'],
                       linewidth=8, label='Centromere', alpha=0.8),
                Line2D([0], [0], color=config['visualization']['colors']['pericentromere'],
                       linewidth=8, label='Pericentromere', alpha=0.6),
                Line2D([0], [0], color=config['visualization']['colors']['rdna_5s'],
                       linewidth=8, label='5S rDNA', alpha=0.9),
                Line2D([0], [0], color=config['visualization']['colors']['rdna_45s'],
                       linewidth=8, label='45S rDNA', alpha=0.9),
            ]
            ax.legend(handles=legend_elements, loc='upper right', fontsize=11, ncol=3)

    # Summary bar plots (if multiple samples)
    if n_samples > 1:
        ax_bars = fig.add_subplot(gs[n_samples + 1, 0])

        # Collect counts
        data = []
        for chrom in chroms:
            row = {'Chromosome': chrom}
            for df, sample_name in zip(sv_data, sample_names):
                chrom_data = df[(df['chromosome'] == chrom) & (df['is_178_multiple'] == True)]
                row[f'{sample_name}_INS'] = len(chrom_data[chrom_data['sv_type'] == 'INS'])
                row[f'{sample_name}_DEL'] = len(chrom_data[chrom_data['sv_type'] == 'DEL'])
            data.append(row)

        plot_df = pd.DataFrame(data)

        x = np.arange(len(chroms))
        width = 0.8 / (n_samples * 2)

        for idx, sample_name in enumerate(sample_names):
            offset = (idx - n_samples/2 + 0.5) * 2 * width
            ax_bars.bar(x + offset - width/2, plot_df[f'{sample_name}_INS'], width,
                       label=f'{sample_name} INS',
                       color=config['visualization']['colors']['insertions'],
                       alpha=0.7)
            ax_bars.bar(x + offset + width/2, plot_df[f'{sample_name}_DEL'], width,
                       label=f'{sample_name} DEL',
                       color=config['visualization']['colors']['deletions'],
                       alpha=0.7)

        ax_bars.set_ylabel('Count', fontsize=12)
        ax_bars.set_title('178bp-Multiple SV Counts by Chromosome',
                         fontsize=14, fontweight='bold')
        ax_bars.set_xticks(x)
        ax_bars.set_xticklabels(chroms)
        ax_bars.set_xlabel('Chromosome', fontsize=12)
        ax_bars.legend(fontsize=9, ncol=n_samples*2)
        ax_bars.grid(axis='y', alpha=0.3)

    plt.savefig(output_file, dpi=config['output']['dpi'], bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Generate genome-wide ideogram visualization',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)

    parser.add_argument('--sv-catalogs', nargs='+', required=True,
                       help='SV catalog TSV files')
    parser.add_argument('--sample-names', nargs='+', required=True,
                       help='Sample names (same order as SV catalogs)')
    parser.add_argument('--config', required=True,
                       help='Configuration YAML file')
    parser.add_argument('--output', required=True,
                       help='Output file path')

    args = parser.parse_args()

    if len(args.sv_catalogs) != len(args.sample_names):
        print("Error: Number of SV catalogs must match number of sample names")
        sys.exit(1)

    # Load config
    config = load_config(args.config)

    # Create output directory
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    # Generate plot
    plot_genome_wide_ideogram(args.sv_catalogs, args.sample_names, config, args.output)


if __name__ == '__main__':
    main()

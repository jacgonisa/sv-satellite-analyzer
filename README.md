# SV Satellite Analyzer

**A comprehensive pipeline for analyzing structural variations in satellite repeats from long-read sequencing data**

Developed for analyzing centromeric satellite repeat dynamics in *Arabidopsis thaliana*, but generalizable to any organism with tandem repeats.

## Overview

This pipeline identifies, quantifies, and visualizes structural variations (insertions and deletions) in satellite repeat arrays, with special focus on:
- Exact satellite unit multiples (e.g., 178bp centromeric repeats)
- Regional enrichment (centromeres, pericentromeres, rDNA loci)
- Single-molecule resolution SV tracking
- Genome-wide distribution patterns
- Mb-normalized comparisons

## Features

- **Memory-efficient**: Processes large BAM files without overwhelming system resources
- **Modular**: Each script performs a specific task and can be run independently
- **Flexible**: Configurable repeat sizes, genomic regions, and filtering thresholds
- **Comprehensive visualization**: Generates publication-ready plots including:
  - Genome-wide ideograms with annotated features
  - Density heatmaps
  - Size distribution histograms
  - Fold-change comparisons
  - Single-molecule SV positions

## Requirements

### Software
- Python 3.7+
- pysam
- pandas
- numpy
- matplotlib
- seaborn
- samtools (for preprocessing)
- seqtk (optional, for downsampling)

### Input Data
- Mapped long-read sequencing data (BAM files, indexed)
- Reference genome coordinates (chromosome lengths)
- Genomic feature annotations (centromeres, pericentromeres, rDNA regions)

## Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/yourusername/sv_satellite_analyzer.git
cd sv_satellite_analyzer

# 2. Install dependencies
conda env create -f environment.yml
conda activate sv_satellite

# 3. Configure your analysis
cp config/config_template.yaml config/my_analysis.yaml
# Edit config/my_analysis.yaml with your paths and parameters

# 4. Run SV detection
python scripts/sv_detection/detect_sv_molecules.py \
    --bam /path/to/sample.bam \
    --output results/sv_molecules \
    --sample sample_name \
    --config config/my_analysis.yaml

# 5. Generate visualizations
python scripts/visualization/plot_genome_wide_ideogram.py \
    --sv-catalog results/sv_molecules/sample_name_sv_catalog.tsv \
    --config config/my_analysis.yaml \
    --output results/plots
```

## Pipeline Structure

### 1. Preprocessing (Optional)
- **downsample_fastq.sh**: Downsample FASTQ files to equalize coverage between samples

### 2. SV Detection
- **detect_sv_molecules.py**: Extract SVs from BAM files with genomic and read-level positions
- **create_rdna_beds.py**: Generate BED files for rDNA regions from GFF annotations

### 3. Analysis
- **compare_samples.py**: Compare SV rates between samples with statistical tests
- **normalize_by_mapped_bases.py**: Calculate region-specific normalization factors
- **analyze_repeat_multiples.py**: Identify exact satellite unit multiples

### 4. Visualization
- **plot_genome_wide_ideogram.py**: Create comprehensive genome-wide summary plots
- **plot_density_heatmap.py**: Generate density tracks across chromosomes
- **plot_size_distributions.py**: Visualize SV size distributions
- **plot_fold_changes.py**: Compare enrichment between samples
- **plot_insertions_vs_deletions.py**: Side-by-side INS/DEL comparison

## Configuration

All analyses are controlled via YAML configuration files in `config/`:

```yaml
# Genome configuration
genome:
  name: "TAIR12"
  chromosomes: ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
  lengths:
    Chr1: 30427671
    Chr2: 19698289
    # ...

# Centromere coordinates
centromeres:
  Chr1: [14853761, 17129892]
  Chr2: [9388656, 11655818]
  # ...

# Analysis parameters
analysis:
  min_sv_size: 50  # Minimum SV size in bp
  satellite_repeat_size: 178  # Centromeric repeat unit size
  pericentromere_extension: 500000  # bp flanking centromere
  window_size: 500000  # For density calculations

# Paths
paths:
  rdna_5s: "/path/to/5s_rdna_regions.bed"
  rdna_45s: "/path/to/45s_rdna_regions.bed"
```

## Output Files

### SV Catalogs
- `{sample}_sv_catalog.tsv`: Per-SV details with genomic and read positions
- `{sample}_sv_molecules.tsv`: Per-molecule summary with all SVs
- `{sample}_insertions_178bp_multiples.tsv`: Satellite-unit multiple insertions
- `{sample}_deletions_178bp_multiples.tsv`: Satellite-unit multiple deletions

### Statistics
- `{sample}_sv_molecule_stats.tsv`: Summary statistics per region
- `178bp_multiples_comparison.tsv`: Cross-sample comparison table
- `mapped_bases_per_region.tsv`: Normalization factors

### Visualizations
- `genome_wide_ideogram_summary.png`: Comprehensive genome overview
- `genome_wide_density_heatmap.png`: Regional density tracks
- `genome_wide_fold_change.png`: Sample comparison in sliding windows
- `178bp_multiples_comparison_normalized.png`: Mb-normalized comparisons
- `insertions_vs_deletions_comparison.png`: INS vs DEL analysis

## Example Analysis

See `docs/example_workflow.md` for a complete walkthrough analyzing centromeric repeat dynamics in *Arabidopsis* ATXR5/6 mutants.

## Citation

If you use this pipeline in your research, please cite:

```
[Your paper citation here]
```

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request with clear description

## License

MIT License - see LICENSE file for details

## Contact

- Issues: https://github.com/yourusername/sv_satellite_analyzer/issues
- Email: your.email@institution.edu

## Acknowledgments

Developed for analyzing centromeric satellite repeat dynamics in *Arabidopsis thaliana* ATXR5/6 mutants using PacBio/ONT long-read sequencing data.

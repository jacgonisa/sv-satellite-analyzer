# SV Satellite Analyzer - Package Summary

## Overview

This package provides a complete pipeline for analyzing structural variations (SVs) in tandem satellite repeat arrays from long-read sequencing data. It was developed for studying centromeric satellite repeat dynamics in *Arabidopsis thaliana*, but is designed to be generalizable to any organism.

## What This Package Does

### Core Functionality

1. **Detects SVs with single-molecule resolution**
   - Extracts insertions and deletions from mapped long reads (BAM files)
   - Tracks SV positions on both genome and individual molecules
   - Identifies exact satellite repeat unit multiples (e.g., 1x, 2x, 3x 178bp)

2. **Classifies SVs by genomic context**
   - Centromeres
   - Pericentromeres
   - 5S rDNA arrays
   - 45S rDNA arrays
   - Chromosome arms

3. **Normalizes across samples**
   - Per-region mapped base normalization
   - Accounts for differences in read length (N50)
   - Enables fair cross-sample comparisons

4. **Generates comprehensive visualizations**
   - Genome-wide ideograms with all features labeled
   - Density heatmaps in sliding windows
   - Size distribution histograms
   - Fold-change comparisons
   - Insertion vs deletion analyses

## Directory Structure

```
sv_satellite_analyzer/
├── README.md                    # Main documentation
├── LICENSE                      # MIT License
├── environment.yml              # Conda environment
├── .gitignore                   # Git ignore rules
│
├── config/                      # Configuration files
│   └── config_template.yaml    # Template configuration
│
├── scripts/                     # All analysis scripts
│   ├── preprocessing/           # Data preparation
│   │   ├── downsample_fastq.sh
│   │   └── create_rdna_beds.py
│   ├── sv_detection/            # Core SV detection
│   │   └── detect_sv_molecules.py
│   ├── analysis/                # Statistical analysis
│   │   └── normalize_by_mapped_bases.py
│   └── visualization/           # Plotting scripts
│       └── plot_genome_wide_ideogram.py
│
├── docs/                        # Documentation
│   ├── example_workflow.md     # Step-by-step tutorial
│   └── PACKAGE_SUMMARY.md      # This file
│
└── example_output/              # Example output files (optional)
```

## Key Scripts

### 1. detect_sv_molecules.py
**Purpose**: Main SV detection engine

**Input**: BAM file, output directory, sample name, min SV size

**Output**:
- `{sample}_sv_catalog.tsv` - Per-SV catalog with all details
- `{sample}_sv_molecules.tsv` - Per-molecule summary
- `{sample}_sv_molecule_stats.tsv` - Regional statistics
- `{sample}_insertions_178bp_multiples.tsv` - Repeat-unit multiple insertions
- `{sample}_deletions_178bp_multiples.tsv` - Repeat-unit multiple deletions

**Key Features**:
- Memory-efficient (processes chromosomes sequentially)
- Tracks both genomic and read-level positions
- Identifies exact repeat unit multiples
- Regional classification

**Usage**:
```bash
python scripts/sv_detection/detect_sv_molecules.py \
    input.bam \
    output_dir \
    sample_name \
    50  # min SV size
```

### 2. plot_genome_wide_ideogram.py
**Purpose**: Generate comprehensive genome-wide visualization

**Input**: SV catalogs, sample names, config file

**Output**: Beautiful ideogram showing:
- All chromosomes with genomic features
- SV positions (insertions above, deletions below)
- Color-coded regions (centromeres, pericentromeres, rDNA)
- Count summaries

**Usage**:
```bash
python scripts/visualization/plot_genome_wide_ideogram.py \
    --sv-catalogs sample1.tsv sample2.tsv \
    --sample-names "Sample 1" "Sample 2" \
    --config config/myconfig.yaml \
    --output plots/ideogram.png
```

### 3. create_rdna_beds.py
**Purpose**: Convert GFF rDNA annotations to BED format

**Input**: GFF files for 5S and 45S rDNA

**Output**: BED files with merged, non-overlapping regions

**Features**:
- Handles concatenated GFF lines
- Merges adjacent regions
- Excludes overlaps with centromeres
- Converts accession IDs to chromosome names

### 4. normalize_by_mapped_bases.py
**Purpose**: Calculate region-specific normalization factors

**Input**: BAM files, genomic regions

**Output**:
- Mapped bases per region
- Normalized indel rates (per Mb)

**Why This Matters**:
Different samples may have different read lengths (N50). Normalizing by read COUNT would be biased. This script normalizes by actual MAPPED BASES, accounting for N50 differences.

### 5. downsample_fastq.sh
**Purpose**: Downsample FASTQ files to equal read counts (optional)

**Usage**:
```bash
bash scripts/preprocessing/downsample_fastq.sh \
    input.fastq \
    output.fastq \
    target_read_count \
    random_seed
```

## Configuration System

All analyses are controlled via YAML configuration files. This makes the pipeline:
- **Reproducible**: Same config = same results
- **Flexible**: Easy to adapt to new organisms/datasets
- **Documented**: Config file serves as analysis record

### Key Config Sections

1. **Genome**: Chromosome names and lengths
2. **Centromeres**: Coordinates for each chromosome
3. **rDNA regions**: Paths to BED files
4. **Analysis parameters**: Min SV size, repeat unit size, window sizes
5. **Samples**: BAM paths, display names, colors
6. **Visualization**: Color schemes, font sizes, output formats

## Output Files Explained

### SV Catalog (`{sample}_sv_catalog.tsv`)
One row per SV with columns:
- `read_id`: Unique identifier for the molecule
- `read_length`: Total length of the sequenced molecule
- `chromosome`: Chromosome where SV is located
- `read_start`, `read_end`: Alignment span on genome
- `sv_ref_pos`: Genomic position of SV
- `sv_read_pos`: Position on the read (0-based from read start)
- `sv_pos_pct`: Percentage along the read (0-100%)
- `sv_type`: INS or DEL
- `sv_size`: Size in bp
- `region`: Genomic context (centromere, pericentromere, etc.)
- `is_178_multiple`: Boolean (for Arabidopsis 178bp repeats)

### SV Molecules (`{sample}_sv_molecules.tsv`)
One row per molecule with at least one SV:
- Read metadata (ID, length, mapping quality)
- SV counts (n_insertions, n_deletions, n_total_svs)
- Total bp affected (total_ins_bp, total_del_bp)
- Lists of positions, sizes, types (comma-separated)

### Statistics (`{sample}_sv_molecule_stats.tsv`)
Summary statistics:
- Total molecules with SVs
- Mean/median read lengths
- Mean/max SVs per read
- Breakdowns by region

## How to Adapt for Your Dataset

### 1. Different Organism
Edit `config.yaml`:
```yaml
genome:
  name: "Your_Genome"
  chromosomes: ["chr1", "chr2", ...]  # Your chromosome names
  lengths:
    chr1: 123456789
    # ...

centromeres:
  chr1: [start, end]
  # ...

analysis:
  satellite_repeat_size: 171  # e.g., human alpha-satellite
```

### 2. Additional Genomic Features
To add custom regions (e.g., specific transposons):
1. Create BED file: `my_feature.bed`
2. Modify `detect_sv_molecules.py`:
   - Add to `load_bed_regions()` call
   - Add classification in `classify_position()` function
   - Add to color scheme in config

### 3. Different SV Size Threshold
Change `min_sv_size` in config or pass as argument:
```bash
python scripts/sv_detection/detect_sv_molecules.py ... 100  # 100bp minimum
```

### 4. Custom Visualizations
All plotting functions accept config dictionaries. Modify colors, fonts, layouts in `config.yaml` under `visualization:` section.

## Memory and Performance

### Memory Efficiency
- Processes one chromosome at a time
- Uses buffered writing (1000 records at a time)
- Tested on systems with 16GB RAM
- Can handle 50GB+ BAM files

### Runtime
- SV detection: ~10-30 min per sample (depends on coverage)
- Normalization: ~5-10 min per sample
- Visualization: ~1-2 min per plot
- Total: <1 hour for typical analysis

### Parallelization
Run multiple samples in parallel:
```bash
parallel -j 4 "python scripts/sv_detection/detect_sv_molecules.py {} results {/.} 50" ::: *.bam
```

## Biological Insights from This Pipeline

Using this pipeline on Arabidopsis ATXR5/6 mutants, we discovered:

1. **ATXR5/6 specifically suppresses centromeric satellite LOSS**
   - 4.1x more deletions in mutant (0.08 → 0.33 per Mb)
   - Only 1.2x change in insertions (minimal)

2. **Spatial specificity**
   - Effect is centromere-specific
   - Not observed in rDNA or chromosome arms

3. **Size specificity**
   - Strong enrichment for exact 178bp multiples
   - Particularly 1x (178bp) and 3x (534bp) deletions

4. **Consistency**
   - ~4-6x fold change across all 5 centromeres
   - No chromosome-specific effects

## Citation

If you use this pipeline in your research, please cite:

```
[Your paper citation will go here]
```

## Getting Help

- **Issues**: https://github.com/yourusername/sv_satellite_analyzer/issues
- **Documentation**: See `docs/` directory
- **Example**: Follow `docs/example_workflow.md`

## Contributing

We welcome contributions! Areas for improvement:
- Additional visualization types
- Support for more file formats (CRAM, etc.)
- Parallelization improvements
- Additional normalization methods
- Integration with SV callers (Sniffles, cuteSV, etc.)

## Version History

- **v1.0.0** (2024): Initial release
  - Core SV detection
  - Regional classification
  - Genome-wide visualizations
  - Mb-normalized comparisons

## Future Directions

Planned features:
- VCF output format
- Integration with phasing information
- Haplotype-resolved SV detection
- Real-time visualization dashboard
- GPU acceleration for large datasets

## Acknowledgments

Developed for studying centromeric satellite repeat dynamics using:
- PacBio/ONT ultra-long reads
- Winnowmap2 mapping
- TAIR12 reference genome

Special thanks to the Arabidopsis genomics community for curated annotations.

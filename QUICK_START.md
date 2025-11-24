# Quick Start Guide

## 🚀 Get Started in 5 Minutes

### 1. Setup Environment

```bash
cd /mnt/ssd-8tb/atrx_china/sv_satellite_analyzer

# Create conda environment
conda env create -f environment.yml

# Activate environment
conda activate sv_satellite
```

### 2. Prepare Your Configuration

```bash
# Copy template
cp config/config_template.yaml config/my_analysis.yaml

# Edit with your settings
nano config/my_analysis.yaml
```

**Minimum required edits:**
- Update `genome:` section with your chromosome names/lengths
- Update `centromeres:` section with your coordinates
- Set `satellite_repeat_size:` to your repeat unit (e.g., 178 for Arabidopsis, 171 for human)

### 3. Run Analysis on One Sample

**Option A: Use the pipeline runner (easiest)**
```bash
bash run_pipeline.sh \
    -b /path/to/sample.bam \
    -s MySample \
    -o results \
    -c config/my_analysis.yaml
```

**Option B: Run scripts individually**
```bash
# Detect SVs
python scripts/sv_detection/detect_sv_molecules.py \
    /path/to/sample.bam \
    results/sv_molecules \
    MySample \
    50  # min SV size

# Visualize
python scripts/visualization/plot_genome_wide_ideogram.py \
    --sv-catalogs results/sv_molecules/MySample_sv_catalog.tsv \
    --sample-names "MySample" \
    --config config/my_analysis.yaml \
    --output results/plots/MySample_ideogram.png
```

### 4. Compare Multiple Samples

```bash
# Run detection for each sample
bash run_pipeline.sh -b sample1.bam -s Sample1 -o results -c config/my_analysis.yaml
bash run_pipeline.sh -b sample2.bam -s Sample2 -o results -c config/my_analysis.yaml

# Compare with visualization
python scripts/visualization/plot_genome_wide_ideogram.py \
    --sv-catalogs results/sv_molecules/Sample1_sv_catalog.tsv \
                  results/sv_molecules/Sample2_sv_catalog.tsv \
    --sample-names "Sample 1" "Sample 2" \
    --config config/my_analysis.yaml \
    --output results/plots/comparison_ideogram.png
```

### 5. Explore Results

```bash
# View SV catalog
less results/sv_molecules/MySample_sv_catalog.tsv

# Check statistics
less results/sv_molecules/MySample_sv_molecule_stats.tsv

# Open plots
open results/plots/MySample_ideogram.png
```

## 📊 What You Get

For each sample:
- **SV Catalog** (`*_sv_catalog.tsv`): Every SV with genomic and read positions
- **Molecule Summary** (`*_sv_molecules.tsv`): Per-read SV summary
- **Statistics** (`*_sv_molecule_stats.tsv`): Regional breakdown
- **178bp Multiples** (`*_insertions_178bp_multiples.tsv`, `*_deletions_178bp_multiples.tsv`): Satellite repeat units
- **Ideogram Plot** (PNG): Beautiful genome-wide visualization

## 🔍 Example Analysis

See complete walkthrough in `docs/example_workflow.md`

## 💡 Tips

**For large datasets:**
```bash
# Run samples in parallel
parallel -j 4 "bash run_pipeline.sh -b {} -s {/.} -o results -c config/my_analysis.yaml" ::: *.bam
```

**Different repeat sizes:**
```yaml
# In config file
analysis:
  satellite_repeat_size: 171  # human alpha-satellite
  # or 120 for mouse minor satellite
  # or 359 for Drosophila
```

**Custom regions:**
Edit `scripts/sv_detection/detect_sv_molecules.py` to add your own genomic features.

## 📚 Documentation

- **Full documentation**: See `README.md`
- **Package details**: See `PACKAGE_SUMMARY.md`
- **Example workflow**: See `docs/example_workflow.md`
- **GitHub setup**: See `GITHUB_SETUP.md`

## 🐛 Troubleshooting

**"BAM index not found"**
```bash
samtools index your_file.bam
```

**"Config file not found"**
```bash
# Use absolute path
bash run_pipeline.sh -c /full/path/to/config.yaml ...
```

**"Memory error"**
The pipeline is designed to be memory-efficient. If you still have issues:
- Process fewer chromosomes at once
- Increase min SV size (`-m 100`)
- Close other applications

## ✨ Next Steps

1. ✅ Run on your first sample
2. 📊 Generate visualizations
3. 🔬 Analyze results (look for enrichment patterns)
4. 📈 Compare multiple samples
5. 📝 Publish your findings!

## 🤝 Need Help?

- Read `PACKAGE_SUMMARY.md` for detailed documentation
- Check `docs/example_workflow.md` for step-by-step tutorial
- Open an issue on GitHub
- Email your questions

**Happy analyzing! 🧬**

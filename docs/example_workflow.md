# Example Workflow: Analyzing ATXR5/6 Centromeric Repeat Dynamics

This example demonstrates the complete workflow for analyzing 178bp centromeric satellite repeat dynamics in *Arabidopsis thaliana* ATXR5/6 mutants using PacBio long-read sequencing data.

## Dataset

- **Organism**: *Arabidopsis thaliana* (TAIR12 reference)
- **Samples**: Col-0 (wildtype) and atxr5/6 (double mutant)
- **Sequencing**: PacBio ultra-long reads (N50 ~45-62 kb)
- **Mapping**: Winnowmap2 (optimized for repetitive regions)
- **Focus**: 178bp centromeric satellite repeats

## Step 1: Prepare Configuration

```bash
cp config/config_template.yaml config/atxr_analysis.yaml
```

Edit `config/atxr_analysis.yaml`:
```yaml
genome:
  name: "TAIR12"
  # ... (see template)

samples:
  wildtype:
    bam: "results/mapping/Col_9day.bam"
    label: "Col-0"
    color: "#2ecc71"

  mutant:
    bam: "results/mapping/atxr56_9day.bam"
    label: "atxr56"
    color: "#9b59b6"

analysis:
  min_sv_size: 50
  satellite_repeat_size: 178  # Arabidopsis centromeric repeat
```

## Step 2: Create rDNA Annotation BEDs

```bash
python scripts/preprocessing/create_rdna_beds.py \
    --gff-5s TAIR12/5s_rdna.gff \
    --gff-45s TAIR12/45s_rdna.gff \
    --output-dir TAIR12/curated_anno
```

This generates:
- `5s_rdna_regions.bed`
- `45s_rdna_regions.bed`

## Step 3: Detect SVs in Each Sample

```bash
# Col-0
python scripts/sv_detection/detect_sv_molecules.py \
    results/mapping/Col_9day.bam \
    results/sv_molecules \
    Col0_full \
    50

# atxr56
python scripts/sv_detection/detect_sv_molecules.py \
    results/mapping/atxr56_9day.bam \
    results/sv_molecules \
    atxr56_full \
    50
```

**Outputs per sample:**
- `{sample}_sv_catalog.tsv` - All SVs with positions
- `{sample}_sv_molecules.tsv` - Per-molecule summaries
- `{sample}_sv_molecule_stats.tsv` - Summary statistics

## Step 4: Generate Visualizations

### Genome-wide Ideogram
```bash
python scripts/visualization/plot_genome_wide_ideogram.py \
    --sv-catalogs results/sv_molecules/Col0_full_sv_catalog.tsv \
                  results/sv_molecules/atxr56_full_sv_catalog.tsv \
    --sample-names "Col-0" "atxr56" \
    --config config/atxr_analysis.yaml \
    --output results/plots/genome_wide_ideogram.png
```

### Additional Plots
For more comprehensive visualizations, see the scripts in `scripts/visualization/`:
- `plot_density_heatmap.py` - Sliding window density tracks
- `plot_fold_change.py` - Sample comparison ratios
- `plot_size_distributions.py` - SV size histograms
- `plot_insertions_vs_deletions.py` - INS vs DEL comparison

## Step 5: Analyze Results

### Key Findings

**178bp-multiple deletions (per Mb, centromere):**
- Col-0: 0.08 deletions/Mb
- atxr56: 0.33 deletions/Mb
- **Fold change: 4.14x**

**178bp-multiple insertions (per Mb, centromere):**
- Col-0: 0.08 insertions/Mb
- atxr56: 0.09 insertions/Mb
- **Fold change: 1.21x**

### Biological Interpretation

1. **ATXR5/6 specifically suppresses centromeric satellite LOSS**
   - Deletions are 4x more frequent in mutant
   - Insertions remain relatively constant

2. **Spatial specificity**
   - Effect is centromere-specific (not seen in rDNA or arms)
   - Consistent across all 5 chromosomes

3. **Repeat unit specificity**
   - Enrichment for exact 178bp multiples (1x, 2x, 3x, etc.)
   - Particularly strong for 1x (178bp) and 3x (534bp) deletions

4. **Mechanism**
   - H3.1 deposition by ATXR5/6 may stabilize satellite arrays
   - Loss of this modification increases recombination/contraction

## Custom Analyses

### Extract Specific Repeat Sizes

```python
import pandas as pd

# Load SV catalog
df = pd.read_csv('results/sv_molecules/atxr56_full_sv_catalog.tsv', sep='\t')

# Filter for 3x 178bp deletions in centromeres
target = df[
    (df['region'] == 'centromere') &
    (df['sv_type'] == 'DEL') &
    (df['sv_size'] == 534) &  # 3x 178bp
    (df['is_178_multiple'] == True)
]

# Get read IDs for extraction
read_ids = target['read_id'].unique()
print(f"Found {len(read_ids)} molecules with 534bp deletions")
```

### Normalize by Custom Regions

If you have additional genomic features (e.g., specific transposons), you can:
1. Create a BED file of your regions
2. Modify `detect_sv_molecules.py` to add your classification
3. Re-run the analysis

## Tips

1. **Memory efficiency**: The pipeline processes chromosomes sequentially to avoid memory issues with large BAMs

2. **Parallelization**: Run multiple samples in parallel:
   ```bash
   parallel -j 4 "python scripts/sv_detection/detect_sv_molecules.py {} results/sv_molecules {/.} 50" ::: *.bam
   ```

3. **Different repeat sizes**: Change `satellite_repeat_size` in config for other organisms:
   - Human alpha-satellite: 171bp
   - Mouse minor satellite: 120bp
   - Drosophila: 359bp

4. **Quality filtering**: Adjust `min_mapping_quality` in config if needed

## Expected Runtime

- SV detection: ~10-30 minutes per sample (depends on BAM size and coverage)
- Visualization: ~1-2 minutes per plot
- Total: <1 hour for 2 samples with full visualization suite

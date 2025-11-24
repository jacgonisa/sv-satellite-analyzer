# Mb-Normalized Analysis Guide

## When to Run Normalization

**After** you've run SV detection on all samples, use this to calculate SVs per Mb for fair comparisons.

## Why Normalize by Mapped Bases?

Different samples may have different read length distributions (N50):
- Sample A: N50 = 45 kb
- Sample B: N50 = 62 kb

If you normalize by **read count**, you underestimate rates in Sample A (shorter reads).
If you normalize by **mapped bases**, you get fair per-Mb rates regardless of read length.

## Quick Example

```bash
# 1. Run SV detection first (for both samples)
python scripts/sv_detection/detect_sv_molecules.py \
    sample1.bam results/sv_molecules Sample1 50

python scripts/sv_detection/detect_sv_molecules.py \
    sample2.bam results/sv_molecules Sample2 50

# 2. Then run normalization
python scripts/analysis/normalize_by_mapped_bases.py \
    sample1.bam \
    results/sv_molecules \
    Sample1 \
    sample2.bam \
    results/sv_molecules \
    Sample2 \
    results/normalized
```

## What It Does

### Step 1: Count Mapped Bases per Region
For each sample, it calculates how many bases are mapped to each genomic region:

```
Region          Sample1 (Mb)    Sample2 (Mb)
Centromere      4,414           5,403
Pericentromere  2,247           2,595
5S rDNA         284             402
45S rDNA        1,840           2,669
Arms            45,953          53,602
```

### Step 2: Count SVs (>=50bp) per Region
From the SV catalog files:

```
Region          Sample1 SVs     Sample2 SVs
Centromere      7,605           15,909
Pericentromere  6,156           7,480
...
```

### Step 3: Calculate Rate per Mb

```
Rate = (SV count / mapped Mb)
```

Example for centromere:
- Sample1: 7,605 SVs / 4,414 Mb = **1.72 SVs per Mb**
- Sample2: 15,909 SVs / 5,403 Mb = **2.94 SVs per Mb**

**Fold change: 2.94 / 1.72 = 1.71x**

## Output Files

Running normalization creates:

```
results/normalized/
├── indel_rates_by_mapped_bases.tsv       # Main results table
├── normalized_comparison.png              # Bar plot comparison
└── normalized_comparison_by_region.png    # Detailed by region
```

## Results Table Format

`indel_rates_by_mapped_bases.tsv`:

```
sample      region          indel_count  mapped_bases_in_region  indels_per_mb  total_mapped_bases
Sample1     centromere      7605         4414400000              1.72           181234000000
Sample1     pericentromere  6156         2247000000              2.74           181234000000
Sample2     centromere      15909        5403000000              2.94           234567000000
...
```

## Visualizations

### Plot 1: Overall Comparison
Shows total indel rate per Mb for each sample across all regions.

### Plot 2: By Region
Bar plots showing indel rates per Mb for each genomic region, with error bars.

## Using with Your Own Data

### Option 1: Two Samples
```bash
python scripts/analysis/normalize_by_mapped_bases.py \
    wildtype.bam results/sv_molecules WT \
    mutant.bam results/sv_molecules MUT \
    results/normalized
```

### Option 2: Multiple Samples
Run pairwise comparisons or modify the script to handle more samples.

### Option 3: Just Calculate Rates
If you only want the rates without comparison plots, modify the script to output just the TSV.

## Advanced: Filter by SV Type

The current script calculates rates for **all SVs (>=50bp)**. To get rates for:

### Only 178bp-multiple deletions:
```python
# In your analysis script
import pandas as pd

# Load SV catalog
df = pd.read_csv('results/sv_molecules/Sample1_sv_catalog.tsv', sep='\t')

# Filter for 178bp-multiple deletions
filtered = df[(df['sv_type'] == 'DEL') & (df['is_178_multiple'] == True)]

# Load mapped bases
bases_df = pd.read_csv('results/normalized/mapped_bases_per_region.tsv', sep='\t')

# Calculate rate
for region in ['centromere', 'pericentromere', 'arms']:
    sv_count = len(filtered[filtered['region'] == region])
    mapped_mb = bases_df[(bases_df['sample'] == 'Sample1') &
                         (bases_df['region'] == region)]['mapped_Mb'].values[0]
    rate = sv_count / mapped_mb
    print(f"{region}: {rate:.2f} 178bp-DEL per Mb")
```

## Integration with Visualization Pipeline

After running normalization, you can:

1. **Update config.yaml** with normalization results
2. **Annotate ideogram plots** with per-Mb rates
3. **Generate publication tables** with normalized values

Example workflow:
```bash
# Full analysis pipeline
bash run_pipeline.sh -b sample1.bam -s Sample1 -o results -c config.yaml
bash run_pipeline.sh -b sample2.bam -s Sample2 -o results -c config.yaml

# Normalization
python scripts/analysis/normalize_by_mapped_bases.py \
    sample1.bam results/sv_molecules Sample1 \
    sample2.bam results/sv_molecules Sample2 \
    results/normalized

# Comparison visualization
python scripts/visualization/plot_genome_wide_ideogram.py \
    --sv-catalogs results/sv_molecules/Sample1_sv_catalog.tsv \
                  results/sv_molecules/Sample2_sv_catalog.tsv \
    --sample-names "Sample 1 (1.72 SVs/Mb)" "Sample 2 (2.94 SVs/Mb)" \
    --config config.yaml \
    --output results/plots/comparison_normalized.png
```

## Tips

1. **Run normalization AFTER** all SV detection is complete
2. **Use the same min_sv_size** threshold for all samples
3. **Check N50 values** - large differences require Mb normalization
4. **Save the mapped_bases TSV** - reuse for different SV filters
5. **Document rates in Methods** - report "per Mb" not "per read"

## Troubleshooting

**"BAM index not found"**
```bash
samtools index sample.bam
```

**"Different numbers of chromosomes detected"**
Make sure all BAMs use the same reference genome.

**"Memory error"**
The script processes chromosomes sequentially, so it should be fine. If issues persist, close other programs.

## Example Results (Arabidopsis ATXR5/6)

From the paper analysis:

### Centromeric SVs (>=50bp):
- Col-0: **1.72 SVs per Mb**
- atxr56: **2.94 SVs per Mb**
- Fold change: **1.71x**

### Centromeric 178bp-multiple deletions:
- Col-0: **0.08 per Mb**
- atxr56: **0.33 per Mb**
- Fold change: **4.14x** ← Key finding!

This normalization revealed that ATXR5/6 specifically suppresses satellite repeat loss.

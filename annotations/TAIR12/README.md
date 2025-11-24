# TAIR12 Genome Annotations

This directory contains genomic feature annotations for *Arabidopsis thaliana* TAIR12 reference genome.

## Files

### Core Genomic Features

- **centromeres.bed** - Centromeric regions with 178bp satellite repeats
- **pericentromeres.bed** - Pericentromeric regions (flanking centromeres)
- **5s_rdna_regions.bed** - 5S ribosomal DNA arrays
- **45s_rdna_regions.bed** - 45S ribosomal DNA arrays

### Format

All files are in BED format (0-based, half-open coordinates: [start, end))

```
chromosome    start       end         [name]
```

## Chromosome Naming

TAIR12 uses RefSeq accession IDs:

| Accession ID | Common Name | Length (bp) |
|--------------|-------------|-------------|
| CP116280.1   | Chr1        | 30,427,671  |
| CP116281.2   | Chr2        | 19,698,289  |
| CP116282.1   | Chr3        | 23,459,830  |
| CP116283.2   | Chr4        | 18,585,056  |
| CP116284.1   | Chr5        | 26,975,502  |

**Note:** rDNA BED files use `Chr1-5` naming for compatibility. Scripts automatically convert between formats.

## Feature Details

### Centromeres

Centromeres are defined by dense arrays of 178bp satellite repeats.

| Chromosome | Start      | End        | Size (Mb) |
|------------|------------|------------|-----------|
| Chr1       | 14,853,761 | 17,129,892 | 2.28      |
| Chr2       | 9,388,656  | 11,655,818 | 2.27      |
| Chr3       | 13,593,457 | 15,740,665 | 2.15      |
| Chr4       | 8,384,369  | 11,037,620 | 2.65      |
| Chr5       | 12,371,485 | 15,175,824 | 2.80      |

### Pericentromeres

Pericentromeric regions are defined as flanking regions on both sides of each centromere.

- Variable size depending on chromosome arm length
- Contain transposons, heterochromatin marks
- Total: 10 regions (2 per chromosome)

### 5S rDNA

5S ribosomal DNA arrays (~120bp repeat units):

| Location | Start     | End       | Size (kb) |
|----------|-----------|-----------|-----------|
| Chr3     | 16,497,982| 16,885,266| 387       |
| Chr4     | 6,980,847 | 7,497,506 | 517       |

### 45S rDNA

45S ribosomal DNA arrays (~9.1kb repeat units containing 18S, 5.8S, and 25S rRNA genes):

| Location | Start | End       | Size (Mb) |
|----------|-------|-----------|-----------|
| Chr2     | 311   | 5,412,003 | 5.41      |
| Chr4     | 1,735 | 3,848,707 | 3.85      |

**Note:** These are the **Nucleolar Organizing Regions (NORs)** in Arabidopsis.

## Coordinates System

- **BED format**: 0-based, half-open [start, end)
- **GFF format**: 1-based, closed [start, end]

When converting:
- BED → GFF: start = start + 1
- GFF → BED: start = start - 1

## Source

- **Centromeres & Pericentromeres**: Curated from TAIR12 reference genome annotation
- **5S rDNA**: Generated from `CC_Col_v2_5s_02022024.gff` (Feb 2024)
- **45S rDNA**: Generated from `CC_Col_v2_45s_02022024.gff` (Feb 2024)

## Usage in Pipeline

These annotations are automatically loaded by the SV detection pipeline when specified in the configuration file:

```yaml
# config.yaml
rdna_regions:
  5s: "annotations/TAIR12/5s_rdna_regions.bed"
  45s: "annotations/TAIR12/45s_rdna_regions.bed"
```

Centromere and pericentromere coordinates are typically defined directly in the config file.

## Citation

If you use these annotations, please cite:

- **TAIR12 Reference**: [Add citation]
- **rDNA annotations**: [Add citation if from publication]

## Notes

- Organellar chromosomes (ChrM, ChrC) are excluded from nuclear analyses
- Some chromosome arms may have overlapping feature annotations (resolved by priority in scripts)
- Coordinates are based on the TAIR12 assembly (RefSeq: GCF_000001735.4)

## Version History

- **2024-11**: Initial curated annotation set
- **2024-02**: rDNA annotations updated

# GitHub Setup Guide

## Quick Setup

Follow these steps to push this package to GitHub:

### 1. Create a new repository on GitHub
1. Go to https://github.com/new
2. Repository name: `sv-satellite-analyzer`
3. Description: "Pipeline for analyzing structural variations in satellite repeats from long-read sequencing"
4. Choose Public or Private
5. **DO NOT** initialize with README (we already have one)
6. Click "Create repository"

### 2. Push your local repository

```bash
cd /mnt/ssd-8tb/atrx_china/sv_satellite_analyzer

# Configure git (if not already done)
git config user.name "Your Name"
git config user.email "your.email@example.com"

# Add all files
git add .

# Make initial commit
git commit -m "Initial commit: SV Satellite Analyzer v1.0

- Core SV detection with single-molecule resolution
- Regional classification (centromeres, pericentromeres, rDNA)
- Comprehensive genome-wide visualizations
- Mb-normalized comparisons
- Configurable for any organism
- Example workflow for Arabidopsis ATXR5/6"

# Add your GitHub repository as remote (replace with your username!)
git remote add origin https://github.com/YOUR_USERNAME/sv-satellite-analyzer.git

# Push to GitHub
git branch -M main
git push -u origin main
```

### 3. Verify on GitHub
Go to your repository URL: `https://github.com/YOUR_USERNAME/sv-satellite-analyzer`

You should see:
- README.md displayed automatically
- All directory structure
- MIT License badge

### 4. Add topics/tags (optional)
On GitHub repository page:
- Click "Settings" → "About" (gear icon)
- Add topics: `bioinformatics`, `long-read-sequencing`, `structural-variation`, `satellite-repeats`, `centromere`, `genomics`, `arabidopsis`, `pacbio`, `nanopore`

### 5. Enable GitHub Pages for documentation (optional)
1. Settings → Pages
2. Source: Deploy from branch
3. Branch: main, folder: /docs
4. Save

## Repository Structure on GitHub

```
yourusername/sv-satellite-analyzer/
├── README.md                       ← Displayed on repo homepage
├── LICENSE                         ← MIT license
├── PACKAGE_SUMMARY.md              ← Comprehensive documentation
├── GITHUB_SETUP.md                 ← This file
├── run_pipeline.sh                 ← Pipeline runner script
├── environment.yml                 ← Conda environment
├── .gitignore                      ← Excludes data files
│
├── config/
│   └── config_template.yaml        ← Example configuration
│
├── scripts/
│   ├── preprocessing/
│   ├── sv_detection/
│   ├── analysis/
│   └── visualization/
│
├── docs/
│   └── example_workflow.md         ← Tutorial
│
└── example_output/
    ├── genome_wide_ideogram_summary.png
    └── 178bp_multiples_comparison_normalized.png
```

## Best Practices

### Commit Messages
Use descriptive commit messages:
```bash
git commit -m "Add: New visualization for density heatmaps"
git commit -m "Fix: Memory issue in large BAM processing"
git commit -m "Docs: Update example workflow with new parameters"
```

### Versioning
When making releases:
```bash
git tag -a v1.0.0 -m "Version 1.0.0: Initial release"
git push origin v1.0.0
```

### Branching Strategy
```bash
# Create feature branch
git checkout -b feature/new-visualization

# Make changes, commit
git add scripts/visualization/new_plot.py
git commit -m "Add: New circular genome plot"

# Push feature branch
git push origin feature/new-visualization

# Create Pull Request on GitHub
# After review, merge to main
```

## Sharing with Collaborators

### Adding Collaborators
1. Repository → Settings → Collaborators
2. Add GitHub usernames
3. They can now push to the repository

### For Users (Not Collaborators)
Users can fork and clone:
```bash
# Fork on GitHub (click "Fork" button)

# Clone their fork
git clone https://github.com/THEIR_USERNAME/sv-satellite-analyzer.git

# They can submit pull requests with improvements
```

## Making a Release

### Create a Release on GitHub
1. Go to "Releases" on repository page
2. Click "Draft a new release"
3. Tag version: `v1.0.0`
4. Release title: "SV Satellite Analyzer v1.0.0"
5. Description:
   ```markdown
   ## Features
   - Single-molecule SV detection
   - Regional classification
   - Genome-wide visualizations
   - Mb-normalized comparisons

   ## Installation
   See README.md

   ## Example Data
   See docs/example_workflow.md
   ```
6. Attach files (optional): Precomputed example outputs
7. Publish release

### Zenodo DOI (for Citation)
1. Link your GitHub to Zenodo: https://zenodo.org/account/settings/github/
2. Enable repository
3. Create a release on GitHub
4. Zenodo automatically creates a DOI
5. Add DOI badge to README

```markdown
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
```

## Maintenance

### Keep Dependencies Updated
```bash
# Update conda environment
conda env update -f environment.yml

# Test pipeline
bash run_pipeline.sh -b test.bam -s test_sample
```

### Issue Template
Create `.github/ISSUE_TEMPLATE.md`:
```markdown
### Description
[Describe the issue]

### Steps to Reproduce
1.
2.
3.

### Expected Behavior
[What should happen]

### Actual Behavior
[What actually happens]

### Environment
- OS:
- Python version:
- Package version:

### Additional Context
[Any other information]
```

## Documentation Website (Optional)

Use GitHub Pages or Read the Docs:

1. **GitHub Pages** (Simple):
   - Put docs in `docs/` folder
   - Enable GitHub Pages in Settings
   - Access at: `https://yourusername.github.io/sv-satellite-analyzer/`

2. **Read the Docs** (Advanced):
   - Connect repository to readthedocs.org
   - Add `docs/conf.py` for Sphinx
   - Access at: `https://sv-satellite-analyzer.readthedocs.io/`

## Continuous Integration (Advanced)

Add `.github/workflows/test.yml`:
```yaml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Set up Python
        uses: actions/setup-python@v2
        with:
          python-version: 3.8
      - name: Install dependencies
        run: |
          pip install pysam pandas numpy matplotlib seaborn pyyaml
      - name: Run tests
        run: |
          python -m pytest tests/
```

## Questions?

- Open an issue on GitHub
- Email: your.email@institution.edu
- Documentation: See PACKAGE_SUMMARY.md

## License

This project is licensed under the MIT License - see LICENSE file for details.

Pan-cancer ZPR1 survival analysis

Files:
- `pan_cancer_zpr1.R` — main R script. Downloads TCGA data via TCGAbiolinks, normalizes with DESeq2 VST (when counts), computes Kaplan–Meier and Cox models per project, and writes `ZPR1_pan_cancer_survival_results.csv`.

Usage:
1. Make sure R (>= 4.x) is installed and you have internet access.
2. From the terminal run:

```bash
Rscript /home/ganesh/pan_cancer_zpr1.R
```

Notes and tips:
- The script will attempt to install missing packages (CRAN and Bioconductor) non-interactively.
- Downloading TCGA data is large. Consider running on a machine with sufficient disk and network bandwidth.
- To publish to GitHub, see instructions below.

Publish to GitHub (example)

If you have the GitHub CLI (`gh`) configured locally and want to create a repo and push:

```bash
cd /home/ganesh
git init
git add pan_cancer_zpr1.R README_pan_cancer_zpr1.md
git commit -m "Add pan-cancer ZPR1 survival analysis script"
# create repo (replace <repo-name> if desired)
gh repo create <your-username>/pan-cancer-zpr1 --public --source=. --remote=origin
git push -u origin main
```

If you prefer to create a repo on github.com manually, create it and then run:

```bash
git remote add origin https://github.com/<your-username>/<repo>.git
git branch -M main
git push -u origin main
```

If you want me to push the code for you, provide the repository URL and confirm you want me to run git commands in this environment (I will not store or request your credentials — the environment must already be authenticated).
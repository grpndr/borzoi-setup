# Borzoi Setup and Inference

A setup script for running Borzoi model inference from SNP input. It automates repository cloning, conda environment creation, package installation, and does a test inference run.

## Contents

- `setup_borzoi.sh` — Main setup script.
- `environment.yml` — Conda environment.
- `export_vars.sh` — Helper script.
- `borzoi_inference.py` — Python script for loading a Borzoi model and scoring SNPs from a VCF-like input file. 
- `snps_dummy.vcf` — Example SNP input file used by the test run.

## Prerequisites

- `conda` installed and available on `PATH`
- `git` installed and available on `PATH`
- Network access on your machine
- macOS or Linux (Or WSL - good luck)

## Setup

From the repository root, run:

```bash
bash setup_borzoi.sh
```

This will:

1. Create a local `borzoi` working directory
2. Clone the `baskerville`, `borzoi`, and `westminster` repositories
3. Create a conda environment named `borzoi` from `environment.yml`
4. Install the cloned repositories into the active environment
5. Configure environment variables via `export_vars.sh`
6. Download example Borzoi pre-trained models
7. Run a sample inference using `borzoi_inference.py`

## Notes

- `borzoi_inference.py` expects a VCF-like table as shown in `snps_dummy.vcf`.
- The example test run in `setup_borzoi.sh` uses the example models bundled with the cloned `borzoi` repository.

## Example output

The output file is a tab-separated table containing the original SNP data plus one or more `L2NORM_<target>_<suffix>` columns for each scored target.

## Legal stuff?

- The three repos cloned by the setup script are https://github.com/calico/baskerville, https://github.com/calico/borzoi and https://github.com/calico/westminster. Please cite them instead of this repo if you use this code. 

- `Borzoi_inference.py` is a modified version of these scripts: https://github.com/calico/borzoi/blob/main/src/scripts/borzoi_sed.py; https://github.com/calico/borzoi/blob/main/src/scripts/borzoi_sad.py. 

- `export_vars.sh` is a collated (and modified) from the scripts of the same name in the three repos.

- Lastly, I work at Earlham Institute and funded by NextGen so all this code belongs to them I guess. I just hope this is useful code, anyone is free to claim whatever from it. 
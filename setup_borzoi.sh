#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Borzoi inference setup script
# Just a wrapper to automate the setup steps in the Borzoi README file.
# In addition there is a test to verify that the model inference runs as intended.
# Note: this script assumes you have conda installed and available on PATH.
# -----------------------------

ENV_NAME="borzoi"
BASE_DIR="borzoi"
# conda environment YAML file
ENV_YML="environment.yml"

echo "Creating base directory: ${BASE_DIR}"
mkdir -p "${BASE_DIR}"
cd "${BASE_DIR}"

echo "Cloning repositories..."

if [ ! -d "baskerville" ]; then
    git clone https://github.com/calico/baskerville.git
else
    echo "baskerville already exists; skipping."
fi

if [ ! -d "borzoi" ]; then
    git clone https://github.com/calico/borzoi.git
else
    echo "borzoi already exists; skipping."
fi

if [ ! -d "westminster" ]; then
    git clone https://github.com/calico/westminster.git
else
    echo "westminster already exists; skipping."
fi

cd ..

echo "Creating conda environment from ${ENV_YML}" #TODO

if ! command -v conda >/dev/null 2>&1; then
    echo "Error: conda is not available on PATH."
    exit 1
fi

if [ ! -f "${ENV_YML}" ]; then
    echo "Error: ${ENV_YML} not found in $(pwd)."
    exit 1
fi

# Enable conda activation inside sh script
CONDA_BASE="$(conda info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"

if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
    echo "Conda environment '${ENV_NAME}' already exists; skipping."
else
    conda env create -n "${ENV_NAME}" -f "${ENV_YML}"
fi

echo "Activating conda environment: ${ENV_NAME}"
conda activate "${ENV_NAME}"

echo "Installing repositories mode..."
cd "${BASE_DIR}"
for repo in baskerville borzoi westminster; do
    # check if they are already installed in the current environment
    if python -c "import ${repo}" &> /dev/null; then
        echo "${repo} is already installed in the current conda environment; skipping."
        ALREADY_INSTALLED=true
    else
        echo "Installing ${repo}..."
        cd "${repo}"
        pip install -e .
        cd ..
    fi
done

cd ..

if [ "${ALREADY_INSTALLED:-false}" = true ]; then
    echo "Baskerville, Borzoi, and Westminster were already installed in the current conda environment; skipping running env_vars.sh"
else
    echo "Exporting repository paths with export_vars.sh" 
    if [ -f "export_vars.sh" ]; then
        source ./export_vars.sh ${BASE_DIR}
    else
        echo "Error: export_vars.sh not found in $(pwd)."
        exit 1
    fi
fi


echo "Downloading Borzoi example pre-trained models..." 

cd ${BASE_DIR}/borzoi

if [ -f "download_models.sh" ]; then
    bash ./download_models.sh
else
    echo "Error: download_models.sh not found in $(pwd)."
    exit 1
fi

cd ..
cd ..

echo "Running Borzoi inference test..." 

if [ -f "borzoi_inference.py" ]; then
    python borzoi_inference.py --snp-vcf-file="snps_dummy.vcf" \
                            --output-file="snp_dummy_predictions.tsv" \
                            --model-file="${BASE_DIR}/borzoi/examples/saved_models/f3c0/train/model0_best.h5" \
                            --model-parameters="${BASE_DIR}/borzoi/examples/params_pred.json" \
                            --targets-file="${BASE_DIR}/borzoi/examples/targets_gtex.txt" \
                            --hg38-fasta="${BASE_DIR}/borzoi/examples/hg38/assembly/ucsc/hg38.fa" \
                            --column-suffix="f3c0_model0"

else
    echo "Error: borzoi_inference.py not found in $(pwd)."
    exit 1
fi

echo "Deactivating conda environment"
conda deactivate

echo "Borzoi inference setup completed successfully. Have fun."

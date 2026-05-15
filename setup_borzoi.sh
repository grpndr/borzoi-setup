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
    echo "Conda environment '${ENV_NAME}' already exists; skipping creation."
else
    conda env create -n "${ENV_NAME}" -f "${ENV_YML}"
fi

echo "Activating conda environment: ${ENV_NAME}"
conda activate "${ENV_NAME}"

echo "Installing repositories mode..."

for repo in baskerville borzoi westminster; do
    echo "Installing ${repo}..."
    cd "${repo}"
    pip install -e .
    cd ..
done

echo "Exporting repository paths with export_vars.sh" #TODO

if [ -f "export_vars.sh" ]; then
    source ./export_vars.sh
else
    echo "Error: export_vars.sh not found in $(pwd)."
    exit 1
fi

echo "Downloading Borzoi example pre-trained models..." #TODO

cd borzoi

if [ -f "download_models.sh" ]; then
    bash ./download_models.sh
else
    echo "Error: download_models.sh not found in $(pwd)."
    exit 1
fi

cd ..

echo "Running Borzoi inference test..." #TODO

if [ -f "test_borzoi_inference.py" ]; then
    python test_borzoi_inference.py

else
    echo "Error: test_borzoi_inference.py not found in $(pwd)."
    exit 1
fi

echo "Deactivating conda environment"
conda deactivate

echo "Borzoi inference setup completed successfully. Have fun."

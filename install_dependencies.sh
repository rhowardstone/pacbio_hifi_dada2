#!/bin/bash
# Installation script for PacBio DADA2 pipeline dependencies

echo "Installing dependencies for PacBio DADA2 pipeline..."
echo "=============================================="

# Check if we're in a conda environment
if [ -z "$CONDA_DEFAULT_ENV" ]; then
    echo "Warning: Not in a conda environment. Consider activating one first."
    echo "Example: conda activate myenv"
    echo ""
fi

# Install cutadapt via conda (if conda is available)
if command -v conda &> /dev/null; then
    echo "Installing cutadapt via conda..."
    conda install -y -c bioconda cutadapt
else
    echo "Conda not found. Please install cutadapt manually:"
    echo "pip install cutadapt"
fi

# Install R packages
echo ""
echo "Installing R packages..."
echo "This may take several minutes..."

Rscript -e '
# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos="https://cloud.r-project.org")
}

# Set Bioconductor version (use latest stable)
BiocManager::install(version = "3.18", ask = FALSE, update = TRUE)

# Install Bioconductor packages
cat("Installing Bioconductor packages...\n")
BiocManager::install(c("dada2", "ShortRead", "Biostrings"), 
                     ask = FALSE, 
                     update = TRUE,
                     force = TRUE)

# Install CRAN packages
cat("Installing CRAN packages...\n")
install.packages(c("optparse", "data.table"), 
                 repos="https://cloud.r-project.org",
                 dependencies = TRUE)

# Check installations
cat("\nChecking installations...\n")
packages <- c("dada2", "ShortRead", "Biostrings", "optparse", "data.table")
for (pkg in packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf("✓ %s installed successfully (version %s)\n", 
                    pkg, packageVersion(pkg)))
    } else {
        cat(sprintf("✗ %s failed to install\n", pkg))
    }
}
'

echo ""
echo "Checking cutadapt installation..."
if command -v cutadapt &> /dev/null; then
    echo "✓ cutadapt installed (version: $(cutadapt --version))"
else
    echo "✗ cutadapt not found in PATH"
fi

echo ""
echo "=============================================="
echo "Installation complete!"
echo ""
echo "To test the pipeline:"
echo "Rscript dada2_pacbio_pipeline.R --help"

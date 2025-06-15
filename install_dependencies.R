#!/usr/bin/env Rscript
# install_dependencies.R
# Installation script for PacBio DADA2 pipeline R dependencies

cat("Installing R dependencies for PacBio DADA2 pipeline...\n")
cat("=====================================================\n\n")

# Function to install packages
install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    if (bioc) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org", dependencies = TRUE)
    }
  } else {
    cat(sprintf("%s is already installed\n", pkg))
  }
}

# Install BiocManager first
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Install CRAN packages
cat("\nInstalling CRAN packages...\n")
install_if_missing("optparse", bioc = FALSE)
install_if_missing("data.table", bioc = FALSE)

# Install Bioconductor packages
cat("\nInstalling Bioconductor packages...\n")
install_if_missing("dada2", bioc = TRUE)
install_if_missing("ShortRead", bioc = TRUE) 
install_if_missing("Biostrings", bioc = TRUE)

# Verify installations
cat("\n\nVerifying installations...\n")
cat("==========================\n")

required_packages <- list(
  list(name = "dada2", bioc = TRUE),
  list(name = "ShortRead", bioc = TRUE),
  list(name = "Biostrings", bioc = TRUE),
  list(name = "optparse", bioc = FALSE),
  list(name = "data.table", bioc = FALSE)
)

all_installed <- TRUE

for (pkg_info in required_packages) {
  pkg <- pkg_info$name
  if (requireNamespace(pkg, quietly = TRUE)) {
    version <- packageVersion(pkg)
    cat(sprintf("✓ %-15s version %s\n", pkg, version))
  } else {
    cat(sprintf("✗ %-15s FAILED TO INSTALL\n", pkg))
    all_installed <- FALSE
  }
}

# Final message
cat("\n")
if (all_installed) {
  cat("All R packages installed successfully!\n")
  cat("\nDon't forget to also install cutadapt:\n")
  cat("  conda install -c bioconda cutadapt\n")
  cat("  OR\n")
  cat("  pip install cutadapt\n")
} else {
  cat("Some packages failed to install. Please check the errors above.\n")
  quit(status = 1)
}

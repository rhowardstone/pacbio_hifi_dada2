# PacBio HiFi DADA2 Pipeline with Read Tracking

A command-line DADA2 pipeline optimized for PacBio HiFi amplicon sequencing data with optional read ID tracking functionality.

## Features

- **Multiple amplicon support**: Pre-configured for V3V4, full-length 16S, Titan, and full operon amplicons
- **Read ID tracking**: Maintains mapping from original read IDs to final ASVs
- **Command-line interface**: Easy to integrate into automated workflows
- **Customizable parameters**: Override default settings for non-standard amplicons
- **Parallel processing**: Multi-threaded support for faster processing

## Installation

### Prerequisites

1. **R (≥4.2.0)** with the following packages:
```R
install.packages("optparse")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("dada2", "ShortRead", "Biostrings"))
install.packages("data.table")
```

2. **cutadapt** (for primer removal):
```bash
# Using conda/mamba
conda install -c bioconda cutadapt

# Or using pip
pip install cutadapt
```

### Download the pipeline
```bash
wget https://raw.githubusercontent.com/[your-repo]/dada2_pacbio_pipeline.R
chmod +x dada2_pacbio_pipeline.R
```

## Usage

### Basic usage
```bash
Rscript dada2_pacbio_pipeline.R -i input_directory -o output_directory -a amplicon_type
```

### Command-line options

| Option | Description | Default |
|--------|-------------|---------|
| `-i, --input` | Input FASTQ files (directory or comma-separated list) | Required |
| `-o, --output` | Output directory | dada2_output |
| `-a, --amplicon` | Amplicon type: V3V4, FL-16S, Titan, full-operon | FL-16S |
| `-t, --threads` | Number of threads | 4 |
| `-p, --pool` | Pooling method: pseudo, independent, pool | pseudo |
| `-r, --track-reads` | Track reads from input to ASV | TRUE |
| `--no-track-reads` | Disable read tracking | - |
| `--minLen` | Override minimum sequence length | Amplicon-specific |
| `--maxLen` | Override maximum sequence length | Amplicon-specific |
| `--maxEE` | Override maximum expected errors | Amplicon-specific |
| `--fwd-primer` | Override forward primer sequence | Amplicon-specific |
| `--rev-primer` | Override reverse primer sequence | Amplicon-specific |
| `--taxonomy` | Path to taxonomy reference database | None |
| `--skip-primers` | Skip primer removal (if already removed) | FALSE |

## Amplicon configurations

| Amplicon | Length range | MaxEE | Forward primer | Reverse primer |
|----------|--------------|-------|----------------|----------------|
| V3V4 | 400-600 bp | 2 | CCTACGGGNGGCWGCAG | GACTACHVGGGTATCTAATCC |
| FL-16S | 1400-1600 bp | 3 | AGAGTTTGATCMTGGCTCAG | AAGGAGGTGATCCAGCCGCA |
| Titan | 2000-2500 bp | 4 | AGAGTTTGATCMTGGCTCAG | Custom 23S primer |
| full-operon | 4000-5000 bp | 5 | AGAGTTTGATCMTGGCTCAG | GCGTGTGTACAAGGCCCGGGAACG |

## Output files

The pipeline generates the following outputs in the specified directory:

- `asv_table.csv` - ASV abundance table (samples × ASVs)
- `asv_sequences.fasta` - Representative sequences for each ASV
- `read_to_asv_mapping.csv` - Read ID to ASV mapping (if tracking enabled)
- `read_tracking.csv` - Read counts through each pipeline step
- `taxonomy_assignments.csv` - Taxonomic assignments (if database provided)
- `error_rates.pdf` - Visualization of learned error rates
- `summary_stats.rds` - Summary statistics in R format

## Examples

### Full-length 16S with taxonomy
```bash
Rscript dada2_pacbio_pipeline.R \
  -i /data/pacbio/fl16s/ \
  -o fl16s_results \
  -a FL-16S \
  -t 16 \
  --taxonomy silva_nr99_v138.1_train_set.fa.gz
```

### V3V4 with custom filtering
```bash
Rscript dada2_pacbio_pipeline.R \
  -i sample1.fastq.gz,sample2.fastq.gz \
  -o v3v4_results \
  -a V3V4 \
  --maxEE 1 \
  --minLen 420
```

### Pre-trimmed reads
```bash
Rscript dada2_pacbio_pipeline.R \
  -i trimmed_reads/ \
  -o results \
  -a Titan \
  --skip-primers \
  --no-track-reads
```

### Custom amplicon
```bash
Rscript dada2_pacbio_pipeline.R \
  -i reads/ \
  -o custom_results \
  -a FL-16S \
  --minLen 800 \
  --maxLen 1200 \
  --fwd-primer GTGCCAGCMGCCGCGGTAA \
  --rev-primer GGACTACHVGGGTWTCTAAT
```

## Troubleshooting

### Common issues

1. **"cutadapt not found"**: Install cutadapt or use `--skip-primers` if primers are already removed

2. **High chimera rates (>30%)**: Usually indicates primers were not properly removed

3. **Memory errors**: Reduce pooling by using `--pool independent` or process samples in batches

4. **No ASVs found**: Check quality filtering parameters, especially `--maxEE`

### Performance tips

- Use at least 8-16 threads for large datasets
- For very large datasets (>100 samples), consider using `--pool independent`
- Pre-filter very long reads (>maxLen) before running to save memory

## Citation

If you use this pipeline, please cite:

- DADA2: Callahan et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods, 13(7), 581-583.
- For PacBio HiFi: Callahan et al. (2019). High-throughput amplicon sequencing of the full-length 16S rRNA gene with single-nucleotide resolution. Nucleic Acids Research, 47(18), e103.

## License

This pipeline is provided as-is under the MIT license.

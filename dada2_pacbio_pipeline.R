#!/usr/bin/env Rscript
# PacBio HiFi DADA2 Pipeline with Read ID Tracking
# Supports multiple amplicon types: V3V4, FL-16S, Titan, full operons

# Parse command line arguments
library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input FASTQ files (comma-separated list or directory path)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="dada2_output",
              help="Output directory [default: %default]", metavar="character"),
  make_option(c("-a", "--amplicon"), type="character", default="FL-16S",
              help="Amplicon type: V3V4, V1V9, FL-16S, Titan, full-operon [default: %default]", metavar="character"),
  make_option(c("-t", "--threads"), type="integer", default=4,
              help="Number of threads [default: %default]", metavar="integer"),
  make_option(c("-p", "--pool"), type="character", default="pseudo",
              help="Pooling method: pseudo, independent, or pool [default: %default]", metavar="character"),
  make_option(c("-r", "--track-reads"), type="logical", default=TRUE, action="store_true",
              help="Track reads from input to ASV [default: %default]"),
  make_option(c("--no-track-reads"), dest="track_reads", action="store_false",
              help="Disable read tracking"),
  make_option(c("--minLen"), type="integer", default=NULL,
              help="Override minimum sequence length", metavar="integer"),
  make_option(c("--maxLen"), type="integer", default=NULL,
              help="Override maximum sequence length", metavar="integer"),
  make_option(c("--maxEE"), type="numeric", default=NULL,
              help="Override maximum expected errors", metavar="numeric"),
  make_option(c("--fwd-primer"), type="character", default=NULL,
              help="Override forward primer sequence", metavar="character"),
  make_option(c("--rev-primer"), type="character", default=NULL,
              help="Override reverse primer sequence", metavar="character"),
  make_option(c("--taxonomy"), type="character", default=NULL,
              help="Path to taxonomy reference database (optional)", metavar="character"),
  make_option(c("--skip-primers"), type="logical", default=FALSE, action="store_true",
              help="Skip primer removal step (if already removed)")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="PacBio HiFi DADA2 Pipeline for ASV calling with read tracking")
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input files must be specified with -i/--input", call.=FALSE)
}

# Load required libraries
suppressPackageStartupMessages({
  library(dada2)
  library(ShortRead)
  library(Biostrings)
  library(data.table)
})

# Configuration for different amplicon types
AMPLICON_CONFIGS <- list(
  "V3V4" = list(
    minLen = 400,
    maxLen = 600,
    maxEE = 2,
    primers = list(
      forward = "CCTACGGGNGGCNGCAG",
      reverse = "GACTACNNGGGTATCTAATCC"
    )
  ),
  "V1V9" = list(
    minLen = 1400,
    maxLen = 1600,
    maxEE = 3,
    primers = list(
      forward = "AGRGTTYGATYMTGGCTCAG",
      reverse = "RGYTACCTTGTTACGACTT"
    )
  ),
  "FL-16S" = list(
    minLen = 1400,
    maxLen = 1600,
    maxEE = 3,
    primers = list(
      forward = "AGRGTTYGATYMTGGCTCAG",
      reverse = "RGYTACCTTGTTACGACTT"
    )
  ),
  "Titan" = list(
    minLen = 2000,
    maxLen = 2500,
    maxEE = 4,
    primers = list(
      forward = "AGRRTTYGATYHTDGYTYAG",
      reverse = "YCNTTCCYTYDYRGTACT"
    )
  ),
  "full-operon" = list(
    minLen = 4000,
    maxLen = 5000,
    maxEE = 5,
    primers = list(
      forward = "AGRGTTTGATYHTGGCTCAG",
      reverse = "CCRAMCTGTCTCACGACG"
    )
  )
)

# Check for required external tools
if (!opt$`skip-primers`) {
  cutadapt_check <- system("which cutadapt", intern = FALSE)
  if (cutadapt_check != 0) {
    stop("cutadapt is required for primer removal but not found in PATH.\n",
         "Install with: conda install -c bioconda cutadapt\n",
         "Or use --skip-primers if primers are already removed.", call.=FALSE)
  }
}

# Process input files
if (dir.exists(opt$input)) {
  # If input is a directory, find all FASTQ files
  input_files <- list.files(opt$input, pattern="\\.(fastq|fq)(\\.gz)?$", 
                           full.names=TRUE, recursive=FALSE)
  if (length(input_files) == 0) {
    stop("No FASTQ files found in input directory: ", opt$input, call.=FALSE)
  }
} else {
  # If input is comma-separated list of files
  input_files <- unlist(strsplit(opt$input, ","))
  # Check if files exist
  missing_files <- input_files[!file.exists(input_files)]
  if (length(missing_files) > 0) {
    stop("Input files not found: ", paste(missing_files, collapse=", "), call.=FALSE)
  }
}

# Override amplicon config if custom parameters provided
if (!is.null(opt$minLen) || !is.null(opt$maxLen) || !is.null(opt$maxEE) || 
    !is.null(opt$`fwd-primer`) || !is.null(opt$`rev-primer`)) {
  
  # Get base config
  config <- AMPLICON_CONFIGS[[opt$amplicon]]
  
  # Override with custom values
  if (!is.null(opt$minLen)) config$minLen <- opt$minLen
  if (!is.null(opt$maxLen)) config$maxLen <- opt$maxLen
  if (!is.null(opt$maxEE)) config$maxEE <- opt$maxEE
  if (!is.null(opt$`fwd-primer`)) config$primers$forward <- opt$`fwd-primer`
  if (!is.null(opt$`rev-primer`)) config$primers$reverse <- opt$`rev-primer`
  
  # Store custom config
  AMPLICON_CONFIGS[["custom"]] <- config
  amplicon_type <- "custom"
} else {
  amplicon_type <- opt$amplicon
}

# Print run information
cat("\n=== PacBio HiFi DADA2 Pipeline ===\n")
cat("Input files:", length(input_files), "FASTQ files\n")
cat("Output directory:", opt$output, "\n")
cat("Amplicon type:", amplicon_type, "\n")
cat("Threads:", opt$threads, "\n")
cat("Pooling method:", opt$pool, "\n")
cat("Read tracking:", opt$`track-reads`, "\n")
if (!is.null(opt$taxonomy)) {
  cat("Taxonomy database:", opt$taxonomy, "\n")
}
cat("\n")

# Main pipeline function
run_pacbio_dada2_pipeline <- function(
  input_files,
  output_dir,
  amplicon_type = "FL-16S",
  threads = 4,
  pool = "pseudo",
  track_reads = TRUE,
  skip_primers = FALSE,
  taxonomy_db = NULL
) {
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get amplicon configuration
  config <- AMPLICON_CONFIGS[[amplicon_type]]
  if (is.null(config)) {
    stop("Unknown amplicon type. Choose from: ", paste(names(AMPLICON_CONFIGS), collapse = ", "))
  }
  
  # Step 1: Read tracking setup
  if (track_reads) {
    read_mapping <- list()
    original_read_ids <- list()
  }
  
  # Step 2: Primer removal with cutadapt (external)
  if (skip_primers) {
    message("Step 1: Skipping primer removal (--skip-primers specified)...")
    trimmed_files <- input_files
  } else {
    message("Step 1: Remove primers using cutadapt...")
    trimmed_files <- file.path(output_dir, paste0("trimmed_", basename(input_files)))
    
    for (i in seq_along(input_files)) {
      cutadapt_cmd <- sprintf(
        "cutadapt -g '^%s' -a '%s$' -o %s %s --discard-untrimmed",
        config$primers$forward,
        as.character(reverseComplement(DNAString(config$primers$reverse))),
        trimmed_files[i],
        input_files[i]
      )
      system(cutadapt_cmd)
    }
  }
  
  # Step 3: Quality filtering with read tracking
  message("Step 2: Quality filtering...")
  filtered_files <- file.path(output_dir, paste0("filtered_", basename(trimmed_files)))
  
  if (track_reads) {
    # Custom filtering to preserve read IDs
    for (i in seq_along(trimmed_files)) {
      # Read sequences with IDs
      reads <- readFastq(trimmed_files[i])
      read_ids <- as.character(id(reads))
      seqs <- as.character(sread(reads))
      quals <- as(quality(reads), "matrix")
      
      # Calculate expected errors
      ee <- rowSums(10^(-quals/10))
      
      # Apply filters
      keep <- ee <= config$maxEE & 
              width(reads) >= config$minLen & 
              width(reads) <= config$maxLen &
              !grepl("N", seqs)
      
      # Store mapping of sequences to original IDs
      filtered_seqs <- seqs[keep]
      filtered_ids <- read_ids[keep]
      
      # Create mapping
      seq_to_ids <- split(filtered_ids, filtered_seqs)
      read_mapping[[basename(input_files[i])]] <- seq_to_ids
      
      # Write filtered file
      writeFastq(reads[keep], filtered_files[i])
    }
  } else {
    # Standard DADA2 filtering
    out <- filterAndTrim(
      trimmed_files, filtered_files,
      maxN = 0,
      maxEE = config$maxEE,
      truncQ = 2,
      minLen = config$minLen,
      maxLen = config$maxLen,
      rm.phix = TRUE,
      compress = TRUE,
      multithread = threads
    )
  }
  
  # Step 4: Learn error rates
  message("Step 3: Learning error rates...")
  setDADAOpt(PACBIO = TRUE)  # Set PacBio-specific parameters
  
  err <- learnErrors(
    filtered_files,
    nbases = 1e8,
    randomize = TRUE,
    multithread = threads
  )
  
  # Plot error rates
  pdf(file.path(output_dir, "error_rates.pdf"))
  plotErrors(err, nominalQ = TRUE)
  dev.off()
  
  # Step 5: Dereplicate and run DADA
  message("Step 4: Dereplicating and denoising...")
  derep <- derepFastq(filtered_files, verbose = TRUE)
  names(derep) <- sapply(strsplit(basename(filtered_files), "_"), `[`, 1)
  
  dada_res <- dada(
    derep, 
    err = err,
    pool = pool,
    multithread = threads
  )
  
  # Step 6: Construct sequence table
  message("Step 5: Constructing sequence table...")
  seqtab <- makeSequenceTable(dada_res)
  
  # Step 7: Remove chimeras
  message("Step 6: Removing chimeras...")
  seqtab_nochim <- removeBimeraDenovo(
    seqtab,
    method = "consensus",
    minFoldParentOverAbundance = 3,
    multithread = threads
  )
  
  # Step 8: Create read-to-ASV mapping if tracking
  if (track_reads) {
    message("Step 7: Creating read-to-ASV mapping...")
    asv_sequences <- colnames(seqtab_nochim)
    read_to_asv_mapping <- data.table()
    
    for (sample in names(read_mapping)) {
      sample_name <- gsub("\\.fastq.*$", "", sample)
      
      # Get ASVs present in this sample
      sample_asvs <- asv_sequences[seqtab_nochim[sample_name, ] > 0]
      
      for (asv_idx in seq_along(sample_asvs)) {
        asv_seq <- sample_asvs[asv_idx]
        
        # Find original read IDs for this sequence
        if (asv_seq %in% names(read_mapping[[sample]])) {
          original_ids <- read_mapping[[sample]][[asv_seq]]
          
          mapping_rows <- data.table(
            sample = sample_name,
            read_id = original_ids,
            asv_id = paste0("ASV", asv_idx),
            asv_sequence = asv_seq,
            abundance_in_sample = seqtab_nochim[sample_name, asv_seq]
          )
          
          read_to_asv_mapping <- rbind(read_to_asv_mapping, mapping_rows)
        }
      }
    }
    
    # Save mapping
    fwrite(read_to_asv_mapping, file.path(output_dir, "read_to_asv_mapping.csv"))
  }
  
  # Step 9: Assign taxonomy (optional)
  if (!is.null(taxonomy_db)) {
    message("Step 8: Assigning taxonomy...")
    tax <- assignTaxonomy(seqtab_nochim, taxonomy_db, multithread = threads, tryRC = TRUE)
    write.csv(tax, file.path(output_dir, "taxonomy_assignments.csv"))
  } else {
    message("Step 8: Skipping taxonomy assignment (no database provided)")
    tax <- NULL
  }
  
  # Step 10: Save outputs
  message("Step 9: Saving outputs...")
  
  # ASV table
  write.csv(seqtab_nochim, file.path(output_dir, "asv_table.csv"))
  
  # ASV sequences
  asv_seqs <- DNAStringSet(colnames(seqtab_nochim))
  names(asv_seqs) <- paste0("ASV", seq_along(asv_seqs))
  writeXStringSet(asv_seqs, file.path(output_dir, "asv_sequences.fasta"))
  
  # Track reads through pipeline
  track <- cbind(
    out,
    sapply(dada_res, getN),
    rowSums(seqtab_nochim)
  )
  colnames(track) <- c("input", "filtered", "denoised", "non-chimeric")
  write.csv(track, file.path(output_dir, "read_tracking.csv"))
  
  # Summary statistics
  summary_stats <- list(
    total_asvs = ncol(seqtab_nochim),
    total_reads_retained = sum(seqtab_nochim),
    reads_per_sample = rowSums(seqtab_nochim),
    chimera_percent = (1 - sum(seqtab_nochim) / sum(seqtab)) * 100
  )
  
  saveRDS(summary_stats, file.path(output_dir, "summary_stats.rds"))
  
  return(list(
    seqtab_final = seqtab_nochim,
    read_mapping = if (track_reads) read_to_asv_mapping else NULL,
    summary = summary_stats
  ))
}

# Run the pipeline
tryCatch({
  results <- run_pacbio_dada2_pipeline(
    input_files = input_files,
    output_dir = opt$output,
    amplicon_type = amplicon_type,
    threads = opt$threads,
    pool = opt$pool,
    track_reads = opt$`track-reads`,
    skip_primers = opt$`skip-primers`,
    taxonomy_db = opt$taxonomy
  )
  
  # Print summary
  cat("\n=== Pipeline Complete ===\n")
  cat("Total ASVs identified:", results$summary$total_asvs, "\n")
  cat("Total reads retained:", results$summary$total_reads_retained, "\n")
  cat("Chimera percentage:", round(results$summary$chimera_percent, 2), "%\n")
  cat("\nOutputs saved to:", opt$output, "\n")
  
}, error = function(e) {
  cat("\nError: Pipeline failed with error:\n")
  cat(conditionMessage(e), "\n")
  quit(status = 1)
})

# Example usage from command line:
# Rscript dada2_pacbio_pipeline.R -i sample1.fastq.gz,sample2.fastq.gz -o results -a FL-16S -t 8
# Rscript dada2_pacbio_pipeline.R -i /path/to/fastq/directory -o results -a Titan --fwd-primer AGAGTTTGATCMTGGCTCAG
# Rscript dada2_pacbio_pipeline.R -i samples/ -o results -a V3V4 --no-track-reads --taxonomy silva_db.fa.gz

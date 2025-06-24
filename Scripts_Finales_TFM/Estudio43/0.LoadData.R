##########################################################################
# Datasets: GSE199243, GSE198323 y GSE185553.
# Article: A Single-Cell Transcriptome Atlas of Glia Diversity in the
# Human Hippocampus across the Lifespan and in Alzheimer’s Disease
# PMID: 36332572
##########################################################################

# Author: Patricia del Carmen Vaca Rubio
# Date: 03/03/25

setwd("/clinicfs/projects/i63/tfm_hipocampo/Study43")
system("mkdir -p scripts")

pacman::p_load(GEOquery, SingleCellExperiment, DropletUtils, stringr, data.table)

###############################################################################
# 1. Download data --> I've downloaded everything by bash scripting (wget,
# ftp links), but here is the equivalent form with R
###############################################################################

studies = c("GSE199243", "GSE198323", "GSE185553")
options(timeout = max(5000, getOption("timeout")))
lapply(studies, getGEOSuppFiles)
system("mkdir -p GSE_files")
system("mv GSE199243 GSE_files/")
system("mv GSE198323 GSE_files/")
system("mv GSE185553 GSE_files/")

system("tar -xvf GSE_files/GSE199243/GSE199243_RAW.tar -C GSE_files/GSE199243/")
system("tar -xvf GSE_files/GSE185553/GSE185553_RAW.tar -C GSE_files/GSE185553/")
system("tar -xvf GSE_files/GSE185553/GSE185553.barcodes.tar.gz -C GSE_files/GSE185553/")

###############################################################################
# 2. Generate the individual SingleCellExperiment objects (one per sample)
###############################################################################

###### Function to read and format data (common for GSE198323 and GSE185553) ######
load_counts <- function(expr_file, barcode_file) {
  counts <- fread(expr_file)
  colnames(counts)[1] <- "Gene"
  
  barcodes <- fread(barcode_file, header = FALSE)$V1
  barcodes[1] <- "Gene"
  
  colnames(counts) <- c(barcodes)
  
  counts_matrix <- as.matrix(counts[, -1, with = FALSE])
  rownames(counts_matrix) <- counts$Gene
  
  return(counts_matrix)
}

###### GSE198323 ######
base_dir <- "GSE_files/GSE198323"
output_dir <- getwd()

# Define the samples that will be processed
samples_one_lib <- c("26", "28", "38", "27", "35", "36")
samples_two_libs <- c("24", "29", "30", "31", "34", "23", "25", "32", "33", "37")
all_samples <- c(samples_one_lib, samples_two_libs)

for (sample in all_samples) {
  sample_prefix <- paste0("GSE198323_Sample", sample)
  
  if (sample %in% samples_two_libs) {
    # Files with two libraries
    expr_file1 <- file.path(base_dir, paste0(sample_prefix, "_lib1.txt.gz"))
    barcode_file1 <- file.path(base_dir, paste0(sample_prefix, "_lib1.barcode.txt.gz"))
    
    expr_file2 <- file.path(base_dir, paste0(sample_prefix, "_lib2.txt.gz"))
    barcode_file2 <- file.path(base_dir, paste0(sample_prefix, "_lib2.barcode.txt.gz"))
    
    if (file.exists(expr_file1) && file.exists(barcode_file1) &&
        file.exists(expr_file2) && file.exists(barcode_file2)) {
      counts_lib1 <- load_counts(expr_file1, barcode_file1)
      counts_lib2 <- load_counts(expr_file2, barcode_file2)
      
      colnames(counts_lib1) <- paste0("GSE198323_Sample", sample, "_",
                                      colnames(counts_lib1), "_lib1")
      colnames(counts_lib2) <- paste0("GSE198323_Sample", sample, "_",
                                      colnames(counts_lib2), "_lib2")
      
      counts_merged <- merge(counts_lib1, counts_lib2, by = "row.names", all = TRUE)
      rownames(counts_merged) <- counts_merged[, 1]
      counts_merged <- counts_merged[, -1]
      counts_merged[is.na(counts_merged)] <- 0
    }
  } else {
    # Files with one library
    expr_file <- file.path(base_dir, paste0(sample_prefix, ".txt.gz"))
    barcode_file <- file.path(base_dir, paste0(sample_prefix, ".barcode.txt.gz"))
    
    if (file.exists(expr_file) && file.exists(barcode_file)) {
      counts_merged <- load_counts(expr_file, barcode_file)
      colnames(counts_merged) <- paste0("GSE198323_Sample", sample, "_",
                                        colnames(counts_merged))
      counts_merged[is.na(counts_merged)] <- 0
    }
  }
  
  # Create and save the SCE object
  if (exists("counts_merged")) {
    sce <- SingleCellExperiment(assays = list(counts = counts_merged))
    colData(sce)$'Sample' <- paste0("Sample", sample)
    colData(sce)$'GSE' <- "GSE198323"
    saveRDS(sce, file.path(output_dir, paste0(sample_prefix, "_SCE.rds")))
  }
}

###### GSE185553 #######
base_dir <- "GSE_files/GSE185553"
output_dir <- getwd()

# Define the samples that will be processed (samples with just one library)
all_samples <- c("18", "19")

for (sample in all_samples) {
  sample_prefix <- paste0("GSE185553_Sample", sample)
  
  # Find the files that include "Sample18" or "Sample19" in their filename 
  sample_files <- list.files(path = base_dir, pattern = paste0("Sample", sample),
                             full.names = TRUE)
  
  # Filter the matrix and corresponding barcode files
  expr_file <- sample_files[grepl("\\.txt\\.gz$", sample_files)]
  barcode_file <- sample_files[grepl("\\.barcode\\.txt$", sample_files)]
  
  # Verify the existence of files and load them
  if (length(expr_file) > 0 && length(barcode_file) > 0) {
    counts_merged <- load_counts(expr_file, barcode_file)
    
    colnames(counts_merged) <- paste0("GSE185553_Sample", sample, "_",
                                      colnames(counts_merged))
    
    counts_merged[is.na(counts_merged)] <- 0
    
    # Create and save the SCE object
    sce <- SingleCellExperiment(assays = list(counts = counts_merged))
    colData(sce)$'Sample' <- paste0("Sample", sample)
    colData(sce)$'GSE' <- "GSE185553"
    saveRDS(sce, file.path(output_dir, paste0(sample_prefix, "_SCE.rds")))
  } else {
    message(paste("Archivos no encontrados para", sample_prefix))
  }
}

###### GSE199243 ###### 
# Create a new function adapted to GSE199243 because the barcode files of this
# dataset have less cells (lower number of barcodes; less number of rows in the
# barcode files) than the number of cells the count matrix has (higher number
# of columns). This may be due to the authors uploading the already quality
# filtered barcodes but the raw count matrix.
base_dir <- "GSE_files/GSE199243"
output_dir <- getwd()

# Function to load and filter the expression data according to filtered barcodes
load_counts43 <- function(expr_file, barcode_file) {
  # Read expression data
  counts <- data.table::fread(expr_file)
  colnames(counts)[1] <- "Gene"
  # Read barcodes
  barcodes <- data.table::fread(barcode_file, header = FALSE)$V1
  barcodes[1] <- "Gene"
  
  # Keep only the count matrix cells that have their corresponding barcode
  counts <- counts[, ..barcodes] 
  
  colnames(counts) <- c(barcodes)
  
  # Convert to matrix (excluding the genes column)
  counts_matrix <- as.matrix(counts[, -1, with = FALSE])
  rownames(counts_matrix) <- counts$Gene
  
  return(counts_matrix)
}

# Define the samples that will be processed (samples with just one library)
all_samples <- c("50", "51", "52")

for (sample in all_samples) {
  sample_prefix <- paste0("GSE199243_Sample", sample)
  # Search for the corresponding expression and barcode files for each sample
  expr_file <- list.files(path = base_dir,
                          pattern = paste0("Sample", sample, ".txt.gz"), full.names = TRUE)
  barcode_file <- list.files(path = base_dir,
                             pattern = paste0("Sample", sample, ".barcode.txt.gz"), full.names = TRUE)
  
  # Verify the existence of files and load them
  if (length(expr_file) > 0 && length(barcode_file) > 0) {
    counts_merged <- load_counts43(expr_file, barcode_file)
    
    colnames(counts_merged) <- paste0("GSE199243_Sample", sample, "_", colnames(counts_merged))
    
    counts_merged[is.na(counts_merged)] <- 0
  
    # Create and save the SCE object
    sce <- SingleCellExperiment(assays = list(counts = counts_merged))
    colData(sce)$'Sample' <- paste0("Sample", sample)
    colData(sce)$'GSE' <- "GSE199243"
    saveRDS(sce, file.path(output_dir, paste0(sample_prefix, "_SCE.rds")))
    
  } else {
    message(paste("Archivos no encontrados para", sample_prefix))
  }
}

###############################################################################
# 3. Generate the final combined SingleCellExperiment
###############################################################################

# Obtain a list of all the SCE filenames (21 samples)
sce_filenames <- list.files(path = ".", pattern="_SCE\\.rds$", full.names=TRUE)

# Load all the SCE objects in a list
sce_objects <- lapply(sce_filenames, readRDS)

# Obtain the complete gene list for each sample
gene_lists <- lapply(sce_objects, rownames)

# Find the common genes for further comparisons
common_genes <- Reduce(intersect, gene_lists)

# Filter each SCE to keep just the common genes
filtered_sce_objects <- lapply(sce_objects, function(sce) sce[common_genes, ])

# Combine the count matrices
combined_counts <- do.call(cbind, lapply(filtered_sce_objects, assay))

# Combine the cell metadata (empty for now)
combined_colData <- do.call(rbind, lapply(filtered_sce_objects, colData))

# Create the combined SCE object
combined_sce <- SingleCellExperiment(
  assays = list(counts = combined_counts),
  colData = combined_colData
)

saveRDS(combined_sce, "combined_sce.rds")

###############################################################################
# 4. Add cell metadata to the combined SCE
###############################################################################
# We will use as metadata the original article's Supplementary Table 1,
# which contains patients information.

metadata(combined_sce)$gene_filtering <- "Common_genes_only"
metadata(combined_sce)$original_genes_per_sample <- sapply(sce_objects, function(x) length(rownames(x)))
sce_filenames <- gsub("^\\./", "", sce_filenames)
sce_filenames <- gsub("\\.rds$", "", sce_filenames)
metadata(combined_sce)$sample <- sce_filenames

# Supplementary Table 1 file
patients_metadata = readxl::read_xlsx("mmc2.xlsx", col_names = TRUE)
colnames(patients_metadata) <- patients_metadata[1,]
patients_metadata <- patients_metadata[-1,]
patients_metadata$'Specimen ID' <- paste0("Sample", patients_metadata$'Specimen ID')
# Substitute the 'NA' for 'Control' in the condition column for more clarity
patients_metadata$'Note'[is.na(patients_metadata$'Note')] <- "Control" 

# Separate 'Note' column into two: condition and Braak Stage level
split_Note <- strsplit(as.character(patients_metadata$'Note'), ",")
patients_metadata$'Note' <- sapply(split_Note, '[', 1)
patients_metadata$'Braak Stage' <- sapply(split_Note, '[', 2)

# Make new covariable group_sex
patients_metadata$Group_sex <- paste(patients_metadata$'Note',
                                     patients_metadata$'Gender', sep = "_")

combined_metadata <- merge(colData(combined_sce), patients_metadata,
                           by.x="Sample", by.y="Specimen ID", sort = FALSE)

# Coincident order because of the previous merge parameter 'sort = FALSE'
rownames(combined_metadata) <- rownames(colData(combined_sce))
# Verification
all.equal(colnames(combined_sce), rownames(colData(combined_sce)))

# Update the colData (cell metadata) of the combined SCE
colData(combined_sce) <- DataFrame(combined_metadata)

# Save the final SCE, ready for the scRNA-seq analysis pipeline
saveRDS(combined_sce,"sce_final_raw.rds")




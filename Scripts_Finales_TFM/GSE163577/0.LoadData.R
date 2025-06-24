#################################################
# Dataset: GSE163577
# Article: A human brain vascular atlas reveals
# diverse mediators of Alzheimer's risk
# PMID: 35165441
#################################################

# Author: Patricia del Carmen Vaca Rubio
# Date: 30/04/25

setwd("/clinicfs/projects/i63/tfm_hipocampo/GSE163577/")
system("mkdir -p scripts")
system("mkdir -p figures")
system("mkdir -p processed_data")

library(GEOquery)
library(DropletUtils)

########################
# 1. Download data     #
########################
# Prevent the download process from crashing
options(timeout = max(5000, getOption("timeout")))

# Donwload the desired files of the study (supp files in this case)
getGEOSuppFiles("GSE163577")

# Rename the directory created by getGEOSuppFiles
system("mv GSE163577/ raw_data/")

# Uncompress the downloaded file and keep the output in the same folder, not current dir
system("tar -xvf raw_data/GSE163577_RAW.tar -C raw_data/")

# Move the original raw tar to current dir
system("mv raw_data/GSE163577_RAW.tar .")

###########################################
# 2. Creation of the SingleCellExperiment #
###########################################
# The study includes samples from hippocampus (Hp) but also from 
# superior frontal cortex (SFC). As we are focusing on hippocampus in this
# thesis, SCE from SFC samples won't be constructed.

# Control hippocampus Patients ID: NCI 1-8 (No Cognitive Impairment)
# AD hippocampus Patients ID: AD 1-9 (Alzheimer's Disease)

# Move away all the SFC samples (containing Ctx in their filenames)
system("mkdir -p raw_ctx_data/")
system("mv raw_data/*_Ctx_* raw_ctx_data/")

# Unzip each of the files (one per sample)
files <- list.files("raw_data", full.names = TRUE) # To save the whole path from current dir
for (f in files) {
  system(paste("tar -xvzf ", f, " -C raw_data/"))
}

# Remove the compressed files left
system("rm raw_data/*.tar.gz")

samples_dir <- list.files("raw_data", full.names = TRUE)
sce_objects <- list()

for (sample in samples_dir) {
   # In this study we have the matrix, features and barcode files
   sce <- DropletUtils::read10xCounts(sample, col.names = TRUE)
   
   sample_name <- paste("Sample", gsub("_filtered.*", "", basename(sample)), sep="_")
   colnames(sce) <- paste(sample_name, colnames(sce), sep="_")
     
   sce_objects[[sample_name]] <- sce
}

# Obtain the gene list of each sample
gene_lists <- lapply(sce_objects, rownames)

# Identify the genes common to all samples to ensure comparability in downstream analysis
common_genes <- Reduce(intersect, gene_lists)

# Filter each SCE to keep just the common genes (rows filtering)
filtered_sce_objects <- lapply(sce_objects, function(sce) sce[common_genes, ])

# Combine all the count matrices
combined_counts <- do.call(cbind, lapply(filtered_sce_objects, assay))

# Combine the cell metadata
combined_colData <- do.call(rbind, lapply(filtered_sce_objects, colData))

# Create the combined SCE
combined_sce <- SingleCellExperiment(
  assays = list(counts = combined_counts),
  colData = combined_colData
)

saveRDS(combined_sce, "sce_gse163577_raw.rds")
saveRDS(sce_objects, "sce_objects.rds")

# combined_sce <- readRDS("sce_gse163577_raw.rds")
# sce_objects <- readRDS("sce_objects.rds")

# Add original metadata as a reference
metadata(combined_sce)$gene_filtering <- "Common_genes_only"
metadata(combined_sce)$original_genes_per_sample <- sapply(sce_objects, function(x) length(rownames(x)))
metadata(combined_sce)$sample <- names(sce_objects)

## Patients metadata ##
# Patients metadata file we will use: Supplementary Table 1 file of the article
# system("wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-04369-3/MediaObjects/41586_2021_4369_MOESM2_ESM.xlsx")
patients_metadata = readxl::read_xlsx("41586_2021_4369_MOESM2_ESM.xlsx", col_names = TRUE)
patients_metadata <- patients_metadata[which(patients_metadata$`Brain region profiled`=="Hippocampus"),]

patients_metadata$'Clinical AD or not'[which((patients_metadata$'Clinical AD or not') == "No")] <- "Control"

# Make new covariable group_sex
patients_metadata$Group_sex <- paste(patients_metadata$'Clinical AD or not', patients_metadata$'Sex', sep = "_")

patients_metadata$`Patient ID` <- gsub(" ", "_", patients_metadata$`Patient ID`)

## colData ##
# Create a covariable to store the batch (treat the H2004 exception)
colData(combined_sce)$batch <- gsub(".*_(\\d{2}_\\d{2}|H2004)_.*", "\\1", rownames(colData(combined_sce)))

# Change row names (filenames) according to their corresponding GSM in GEO
# for example, GSM4982083_01_10_C_... -> Control 1 Hpc
rownames(colData(combined_sce)) <- gsub("^Sample_01_10_C", "NCI_1", rownames(colData(combined_sce)))
rownames(colData(combined_sce)) <- gsub("^Sample_01_05_AD", "AD_1", rownames(colData(combined_sce)))
rownames(colData(combined_sce)) <- gsub("^Sample_01_06_C", "NCI_2", rownames(colData(combined_sce)))
rownames(colData(combined_sce)) <- gsub("^Sample_01_09_C", "NCI_3", rownames(colData(combined_sce)))
rownames(colData(combined_sce)) <- gsub("^Sample_02_11_C", "NCI_4", rownames(colData(combined_sce)))
rownames(colData(combined_sce)) <- gsub("^Sample_02_01_C", "NCI_5", rownames(colData(combined_sce)))
rownames(colData(combined_sce)) <- gsub("^Sample_02_04", "AD_2", rownames(colData(combined_sce)))
rownames(colData(combined_sce)) <- gsub("^Sample_03_02_AD", "AD_3", rownames(colData(combined_sce)))
rownames(colData(combined_sce)) <- gsub("^Sample_03_03_C", "NCI_6", rownames(colData(combined_sce)))
rownames(colData(combined_sce)) <- gsub("^Sample_03_09_C", "NCI_7", rownames(colData(combined_sce)))
rownames(colData(combined_sce)) <- gsub("^Sample_04_11_AD", "AD_4", rownames(colData(combined_sce)))
rownames(colData(combined_sce)) <- gsub("^Sample_04_13_AD", "AD_5", rownames(colData(combined_sce)))
rownames(colData(combined_sce)) <- gsub("^Sample_04_15_AD", "AD_6", rownames(colData(combined_sce)))
rownames(colData(combined_sce)) <- gsub("^Sample_04_09_C", "NCI_8", rownames(colData(combined_sce)))
rownames(colData(combined_sce)) <- gsub("^Sample_05_10_AD", "AD_7", rownames(colData(combined_sce)))
rownames(colData(combined_sce)) <- gsub("^Sample_05_04_AD", "AD_8", rownames(colData(combined_sce)))
rownames(colData(combined_sce)) <- gsub("^Sample_H2004_AD", "AD_9", rownames(colData(combined_sce)))

# Make a new patient ID column in colData for later merge
colData(combined_sce)$Patient_ID <- gsub("^(NCI_[0-9]{1}|AD_[0-9]{1}|H2004)_.*", "\\1", rownames(colData(combined_sce)))

combined_metadata <- merge(colData(combined_sce), patients_metadata,
                           by.x="Patient_ID", by.y="Patient ID", sort = FALSE)
# by.x to indicate the first df's column and by.y to indicate the second df's column to merge by

# Coincident order because of the previous merge parameter 'sort = FALSE'
rownames(combined_metadata) <- rownames(colData(combined_sce))
# Verification
all.equal(colnames(combined_sce), rownames(colData(combined_sce)))

# Update the colData (cell metadata) of the combined SCE
colData(combined_sce) <- DataFrame(combined_metadata)

# Save the final SCE, ready for the scRNA-seq analysis pipeline
saveRDS(combined_sce,"sce_final_raw.rds")
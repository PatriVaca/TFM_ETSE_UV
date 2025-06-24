#################################################
# Dataset: syn52293442
# Article: Single-cell multiregion dissection of
# Alzheimer’s disease
# DOI: https://doi.org/10.1038/s41586-024-07606-7
#################################################

# Author: Patricia del Carmen Vaca Rubio
# Date: 05/05/25

# Load packages
library(SeuratObject)
library(Seurat)
library(SingleCellExperiment)

setwd("/clinicfs/projects/i63/tfm_hipocampo/syn52293442_hc/")
system("mkdir -p scripts")
system("mkdir -p figures")
system("mkdir -p processed_data")

## Read data
# Convert Seurat original object to SingleCellExperiment
Hippocampus <- readRDS("data/Hippocampus.rds")
sce <- as.SingleCellExperiment(Hippocampus)
coldata <- as.data.frame(colData(sce))
dim(sce)

metadata <- read.table("metadata/Supplementary_Table_1_sample_metadata.txt", sep = "\t", header = T)

# In this thesis, we will focus just on hippocampus
metadata <- metadata[metadata$region == "HC",]
# 48 HC samples, 8 columns

# Explore all the articles' metadata
metadata2 <- read.table("metadata/MIT_ROSMAP_Multiomics_assay_multiome_metadata.csv", sep = ",", header = T)
metadata3 <- read.table("metadata/MIT_ROSMAP_Multiomics_assay_RNAseq_metadata.csv", sep = ",", header = T)
metadata4 <- read.table("metadata/MIT_ROSMAP_Multiomics_assay_snATACseq_metadata.csv", sep = ",", header = T)
metadata5 <- read.table("metadata/MIT_ROSMAP_Multiomics_assay_snRNAseq_metadata.csv", sep = ",", header = T)
metadata6 <- read.table("metadata/MIT_ROSMAP_Multiomics_biospecimen_metadata.csv", sep = ",", header = T)
metadata7 <- read.table("metadata/MIT_ROSMAP_Multiomics_individual_metadata.csv", sep = ",", header = T)
rosmap_metadata <- read.table("metadata/ROSMAP_clinical.csv", sep = ",", header = T)
# rosmap_metadata <- rosmap_metadata[rosmap_metadata$individualID %in% metadata6$individualID,]

## Metadata creation
metadata <-  merge(metadata, unique(metadata7[,c(1,19)]), 
                    by.x = "subject", 
                    by.y = "subject",
                    all.x = TRUE)

# sum(duplicated(metadata$subject))
# sum(duplicated(metadata$individualID))

metadata <-  merge(metadata, rosmap_metadata, 
                    by.x = "individualID", 
                    by.y = "individualID")

# Transform msex covariable: 1 = Male, 0 = Female
# (https://www.synapse.org/Synapse:syn3191087)
# We will leave msex.y with the original values
metadata$msex.x[metadata$msex.x=="1"] <- "Male"
metadata$msex.x[metadata$msex.x=="0"] <- "Female"

# Transform pathAD variable
metadata$pathAD[metadata$pathAD=="non-AD"] <- "Control"

# Make new covariable Group_sex
metadata$Group_sex <- paste(metadata$pathAD,
                            metadata$msex.x, sep = "_")

write.table(metadata, file = "processed_data/final_metadata_HC.txt")

# They have projid column in common
combined_metadata <- merge(colData(sce), metadata, sort = FALSE)

# Coincident order because of the previous merge parameter 'sort = FALSE'
rownames(combined_metadata) <- rownames(colData(sce))

# Verification
all.equal(colnames(sce), rownames(colData(sce)))

# Update the colData (cell metadata) of the combined SCE
colData(sce) <- DataFrame(combined_metadata)

saveRDS(sce, "processed_data/sce.rds")




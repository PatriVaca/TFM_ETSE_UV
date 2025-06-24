# Single cell RNA-seq data analysis
# 8.Differential Gene Expression Analysis
# Author: Patricia del Carmen Vaca Rubio
# Date: 05/25

#----------- LOAD PACKAGES AND DATA -----------#
pacman::p_load(
  # Core
  SingleCellExperiment,
  SummarizedExperiment,
  S4Vectors,
  # DEG analysis
  MAST,
  # Data manipulation
  dplyr,
  tidyverse,
  # Disperse matrices
  Matrix,
  DelayedArray,
  # Annotation
  org.Hs.eg.db
)

setwd("/clinicfs/projects/i63/tfm_hipocampo/GSE163577/")
system("mkdir -p processed_data/8.DEG")

# Load the data
sce.annot <- readRDS("processed_data/scegv_dw_clusters_filtered_annot1.rds")

# Check if there is any Gene Symbol or ENSEMBLID = NA
table(is.na(rownames(sce.annot))) 
table(is.na(rowData(sce.annot)$gene_symbol))

# Filter out genes with NA as ENSEMBLID
sce.annot <- sce.annot[!(is.na(rownames(sce.annot))), ]

# Filter out genes with NA as Gene Symbol
sce.annot <- sce.annot[!is.na(rowData(sce.annot)$gene_symbol), ]

sum(duplicated(rowData(sce.annot)$gene_symbol))
sum(duplicated(rownames(sce.annot)))

# Filter out genes with duplicated ENSEMBLIDS (same id for a single gene symbol) or Gene symbols 
sce.annot <- sce.annot[!duplicated(rowData(sce.annot)$gene_symbol), ]
sce.annot <- sce.annot[!duplicated(rownames(sce.annot)), ]

rowData(sce.annot)$ENSEMBL <- rownames(sce.annot)

saveRDS(sce.annot, "processed_data/scegv_dw_clusters_filtered_annot1.rds")

#----------- FUNCTIONS -----------#

# MAST requires TPM normalization
## Calculate log2(TPM + 1) counts
tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

#--------- Select cell type ---------#

cell_types <- unique(sce.annot$SCINA_label)
for (i in 1:length(cell_types)){
  
  cell_type <- cell_types[i]
  print(paste0("Vamos con ",cell_type))
  
  sce.dge <- sce.annot[,sce.annot$SCINA_label == cell_type]
  
  dir.create(paste0("./processed_data/8.DEG/", cell_type,"/"), showWarnings=F)
  
  
  #----------- Obtain each transcript length for TPM normalization -----------#
  print("Genes iniciales")
  print(dim(sce.dge))
  
  rownames(sce.dge) <- unlist(sapply(strsplit(rownames(sce.dge), ".", fixed = TRUE),
                                     `[`, 1, simplify=FALSE))
  
  ## Load the transcript lengths file
  transcript <- read.delim("../MAST/transcript_length_borja.txt", header = TRUE, sep = "\t")
  # Remove duplicates
  transcript <- transcript[!duplicated(transcript$Gene.stable.ID),]
  
  # Keep only the common genes between our sce and our lengths file
  # Keep the common genes IN the sce
  x <- rownames(sce.dge) %in% transcript$Gene.stable.ID
  sce.dge  <- sce.dge[x,]
  
  # Keep the common genes IN the lengths file
  x <- transcript$Gene.stable.ID %in% rownames(sce.dge)
  transcript <- transcript[x,]
  
  # Verification
  transcript <- transcript[match(rownames(sce.dge), transcript$Gene.stable.ID),]
  all.equal(rownames(sce.dge), transcript$Gene.stable.ID)
  
  # Add the transcripts length information to our sce
  rowData(sce.dge)[,"gene.start"] <- transcript$Gene.start..bp.
  rowData(sce.dge)[,"gene.end"] <- transcript$Gene.end..bp.
  rowData(sce.dge)[,"gene.length"] <- transcript$Gene.end..bp. - transcript$Gene.start..bp. + 1
  
  x <- as.data.frame(rowData(sce.dge))
  combination <- merge(x, transcript, by.x = "ENSEMBL", by.y = "Gene.stable.ID")
  rowData(sce.dge) <- DataFrame(combination)
  
  print("Genes que se incluirán en el analisis DEG")
  print(dim(sce.dge))
  
  #----------- TPM normalization -----------#
  
  ### Compute log2(TPM + 1) values
  # TPM is very similar to RPKM and FPKM. The only difference is the order of operations. Here’s how you calculate TPM:
  ## Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
  ## Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
  ## Divide the RPK values by the “per million” scaling factor. This gives you TPM.
  
  # table(is.na(rowData(sce.dge)$"gene.length"))
  
  # Operate over the raw counts matrix (counts, not logcounts)
  assay(sce.dge, "tpm") <- tpm(counts = assay(sce.dge, "counts"), len = rowData(sce.dge)[,"gene.length"])
  sce.dge@assays@data[["log2tpm"]] <- log1p(assay(sce.dge, "tpm")) / log(2)
  
  # Scale genes
  cdr2_cells <- colSums(assay(sce.dge, "log2tpm")>0)
  colData(sce.dge)$cngeneson <- scale(cdr2_cells)  
  
  #----------- MAST modelling -----------#
  
  print(paste0("Vamos a modelar para ", cell_type))
  
  # 1.- Build the MAST object
  sce.assay_cells <- SceToSingleCellAssay(sce.dge , check_sanity = FALSE)
  rowData(sce.assay_cells)$primerid <- rownames(sce.assay_cells)
  
  # 2.- Prepare the Hurdle model and adjust it to the data
  zlmCond <- zlm(~0+Group_sex+cngeneson+Patient_ID, sce.assay_cells, exprs_values = "log2tpm", method = "bayesglm", ebayes = TRUE, parallel = TRUE)
  # Save the model
  saveRDS(zlmCond, file = paste0("./processed_data/8.DEG/", cell_type,"/zlmCond.rds"))
  # Save the dataset
  saveRDS(sce.dge, file = paste0("./processed_data/8.DEG/", cell_type,"/sce_dge.rds"))
}

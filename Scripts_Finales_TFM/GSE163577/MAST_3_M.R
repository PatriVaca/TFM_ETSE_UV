# Single cell RNA-seq data analysis
# 8.Differential Gene Expression Analysis. Male contrast.
# Author: Patricia del Carmen Vaca Rubio
# Date: 05/25

#----------- LOAD PACKAGES AND DATA -----------#
# Load packages
pacman::p_load(SingleCellExperiment, scater, scran, dplyr, BiocParallel, 
               scDblFinder, AnnotationHub, tidyverse, patchwork, ggvenn, broom,
               kableExtra,HDF5Array, viridis, DropletUtils,
               devtools, celda, ggvenn, bluster, cluster, ggpubr, MAST, limma, edgeR)

setwd("/clinicfs/projects/i63/tfm_hipocampo/GSE163577/")

# Load the data
sce.annot <- readRDS("processed_data/scegv_dw_clusters_filtered_annot1.rds")

#----------- FUNCTIONS -----------#

# MAST requires TPM normalization
## Calculate log2(TPM + 1) counts
tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

## Extract results from lrtest
extract_results <- function(zz) {
  
  results <- data.frame(rownames(zz))
  colnames(results) <- "gene_id"
  
  # Selection of
  results[,"lambda"] <- data.frame(zz[,3,1]) #Statistic
  results[,"p.value"] <- data.frame(zz[,3,3]) #p.value hurdle
  results[,"p.adjusted"] <- p.adjust(results[,"p.value"], "BH") #p.adjusted BH
  
  return(results)
}

## logFC standard calculation
LFC.standard <- function(data, results, case, control) {
  
  # Calculate cpms and select data of each group. Sum pseudocount factor (1)
  cpm_case <- calculateCPM(counts(data))[results[,"gene_id"],data$Group_sex==case]  + 1
  cpm_control <- calculateCPM(counts(data))[results[,"gene_id"],data$Group_sex==control]  + 1
  
  # Calculate logFC
  LFC <- as.data.frame(log2(rowMeans(cpm_case)/rowMeans(cpm_control)))
  colnames(LFC) <- "logFC.standard"
  LFC$gene_id <- rownames(cpm_case)
  
  return(LFC)
}

## logFC MAST calculation
LFC.MAST <- function(zlm, control, case, results) {
  coefnames <- colnames(zlm@LMlike@modelMatrix)
  
  # Constrast 0
  contrast0 <- setNames(rep(0, length(coefnames)), coefnames)
  contrast0[control] <- 1
  contrast1 <- data.frame("case"=setNames(rep(0, length(coefnames)), coefnames))
  contrast1[case,] <- 1
  
  control <- Hypothesis(control, colnames(zlm@LMlike@modelMatrix))
  case <-Hypothesis(case, colnames(zlm@LMlike@modelMatrix))
  
  fc <- getLogFC(zlm, contrast0, contrast1)
  sum <-summary(zlm, logFC=fc)$datatable
  sum <- sum[contrast == "case" & component=='logFC',-c(2,3)]
  sum$SE <- (sum[, "ci.hi"] - sum[, "ci.lo"])/ 3.92
  colnames(sum)[c(1,4)] <- c("gene_id","logFC")
  sum[is.na(sum)] <- 0
  
  results <- merge(results, sum, by ="gene_id", sort=F)
  
  return(results)
}

#--------- Select cell type ---------#

cell_types = unique(sce.annot$SCINA_label)

for (i in 1:length(cell_types)){
  
  cell_type <- cell_types[i]
  print(paste0("Vamos con ",cell_type))
  
  dir.create(paste0("./processed_data/8.DEG/", cell_type,"/"), showWarnings=F)
  path <- paste0("./processed_data/8.DEG/", cell_type,"/")
  
  # Load the model
  zlmCond <- readRDS(paste0("./processed_data/8.DEG/", cell_type,"/zlmCond.rds"))
  # Load the dataset
  sce.dge <- readRDS(paste0("./processed_data/8.DEG/", cell_type,"/sce_dge.rds"))
  
  #----------- Hypothesis contrast -----------#
  
  ### CaseMale - ControlMale ###
  zlmCond_cells_M <- lrTest(zlmCond, Hypothesis("Group_sexAD_M - Group_sexControl_M"))
  saveRDS(zlmCond_cells_M, file = paste0("./processed_data/8.DEG/", cell_type,"/zlmCond_cellsM.rds"))
  results <- extract_results(zz = zlmCond_cells_M)
  
  # Calculate LFCs
  LFC.sM <- LFC.standard(data = sce.dge, results = results, case = "AD_M", control = "Control_M")
  
  LFC.mM <- LFC.MAST(zlm = zlmCond, control = "Group_sexControl_M", case = "Group_sexAD_M", results = results)
  
  results <- LFC.mM
  results$logFC.CPM <- LFC.sM$logFC.standard
  
  results$SYMBOL <- rowData(sce.dge)$gene_symbol
  results <- results[order(results$p.adjusted), ]
  
  write.table(results, file=paste0(path,cell_type,"res_M.tsv"), col.names=TRUE, row.names=FALSE, sep="\t")
  write.table(results[results$p.adjusted < 0.05,], file=paste0(path,cell_type,"res_M_sig.tsv"), col.names=TRUE, row.names=FALSE, sep="\t")
  print("Male completado")
}

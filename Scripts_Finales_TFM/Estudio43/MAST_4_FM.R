# Single cell RNA-seq data analysis
# 8.Differential Gene Expression Analysis. Sex Contrast
# Author: Patricia del Carmen Vaca Rubio
# Date: 05/25

#----------- LOAD PACKAGES AND DATA -----------#
# Load the packages
pacman::p_load(SingleCellExperiment, scater, scran, dplyr, BiocParallel, 
               scDblFinder, AnnotationHub, tidyverse, patchwork, ggvenn, broom,
               kableExtra,HDF5Array, viridis, DropletUtils,
               devtools, celda, ggvenn, bluster, cluster, ggpubr, MAST, limma, edgeR)

setwd("/clinicfs/projects/i63/tfm_hipocampo/Study43/")

# Load the data
sce.annot <- readRDS("processed_data/scegv_dw_clusters_filtered_annot2.rds")

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

cell_types <- unique(sce.annot$tipo.celular)

for (i in 1:length(cell_types)){
  
  cell_type <- cell_types[i]
  
  print(paste0("Vamos con ",cell_type))
  
  dir.create(paste0("./processed_data/8.DEG/", cell_type,"/"), showWarnings=F)
  path <- paste0("./processed_data/8.DEG/", cell_type,"/")
  
  # Load the model
  zlmCond <- readRDS(paste0("./processed_data/8.DEG/", cell_type,"/zlmCond.rds"))
  # Load the dataset
  sce.dge <- readRDS(paste0("./processed_data/8.DEG/", cell_type,"/sce_dge.rds"))
  
  print("Comenzamos")
  
  #----------- Hypothesis contrast -----------#
  
  ### (CaseFemale - ControlFemale) - (CaseMale - ControlMale)  ###
  
  # Test hypothesis
  zlmCond_cells_FM <- lrTest(zlmCond, Hypothesis("(Group_sexAlzheimers_Disease_Female - Group_sexControl_Female) -
                                                 (Group_sexAlzheimers_Disease_Male - Group_sexControl_Male)"))
  saveRDS(zlmCond_cells_FM, file = paste0("./processed_data/8.DEG/", cell_type,"/zlmCond_cellsFM.rds"))

  # Load each contrast model
  # zlmCond_cells_FM = readRDS(paste0("./processed_data/8.DEG/", cell_type,"/zlmCond_cellsFM.rds"))
  results_2 <- extract_results(zz = zlmCond_cells_FM)

  zlmCond_cellsF <- readRDS(paste0("./processed_data/8.DEG/", cell_type,"/zlmCond_cellsF.rds"))
  resultsF <- extract_results(zz = zlmCond_cellsF)
  LFC.sF <- LFC.standard(data = sce.dge, results = resultsF, case = "Alzheimers_Disease_Female", control = "Control_Female")

  zlmCond_cellsM <- readRDS(paste0("./processed_data/8.DEG/", cell_type,"/zlmCond_cellsM.rds"))
  resultsM <- extract_results(zz = zlmCond_cellsM)
  
  # Calculate LFCs
  LFC.sM <- LFC.standard(data = sce.dge, results = resultsM, case = "Alzheimers_Disease_Male", control = "Control_Male")

  LFC.mFM <- LFC.MAST(zlm = zlmCond,
                     control = c("Group_sexControl_Female","Group_sexAlzheimers_Disease_Male"),
                     case =  c("Group_sexControl_Male","Group_sexAlzheimers_Disease_Female"),
                     results = results_2)

  logFC_FM <- merge(LFC.sF, LFC.sM, by ="gene_id")
  logFC_FM$logFC.cpm <- logFC_FM[,"logFC.standard.x"] - logFC_FM[,"logFC.standard.y"]
  logFC_FM <- logFC_FM[,c("gene_id","logFC.cpm")]

  results <- LFC.mFM
  results$logFC.CPM <- logFC_FM$logFC.cpm

  results$SYMBOL <- rowData(sce.dge)$SYMBOL
  results <- results[order(results$p.adjusted), ]

  write.table(results, file=paste0(path,cell_type,"res_FM.tsv"), col.names=TRUE, row.names=FALSE, sep="\t")
  write.table(results[results$p.adjusted < 0.05,], file=paste0(path,cell_type,"res_FM_sig.tsv"), col.names=TRUE, row.names=FALSE, sep="\t")

  print("Sex completado")
}


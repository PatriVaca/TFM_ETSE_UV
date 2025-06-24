# Single cell RNA-seq data analysis
# 1.Quality Control (QC)
# Author: Patricia del Carmen Vaca Rubio
# Date: 04/25

setwd("/clinicfs/projects/i63/tfm_hipocampo/GSE163577/")
system("mkdir -p figures/1.QC")

library(Matrix)
library(scuttle)
library(scater)
library(scDblFinder)
library(dplyr)
library(tidyr)
# options(warn = -1) # Avoid warnings to interrupt the execution

# Load the combined SingleCellExperiment
sce.raw <- readRDS("sce_final_raw.rds")

#--------------------- QUALITY FILTERING ------------------------------------#
## Quality analysis and filtering (thresholds according to the original paper) ##

# Remove undetected genes
sce.raw <- sce.raw[Matrix::rowSums(counts(sce.raw) > 0) > 0, ]

# Doublets identification per patient
print("Número de doublets identificados")
assay(sce.raw, "counts") <- as(Matrix(as.matrix(counts(sce.raw))), "dgCMatrix")
set.seed(100)
sce.raw <- scDblFinder(sce.raw, samples = "Sample")

# Mitochondrial genes
is.mito <-   grepl("^MT-", rownames(sce.raw))
mito_genes <- rownames(sce.raw[is.mito,]) 

# Ribosomal genes
is.ribo <-   grepl("^RP[SL]", rownames(sce.raw))
ribo_genes <- rownames(sce.raw[is.ribo,]) 

# Spike-ins (ERCC)
is.ercc <-   grepl("^ERCC", rownames(sce.raw))
ercc_genes <- rownames(sce.raw[is.ercc,]) 

# QC metrics
sce.raw <- scuttle::addPerCellQC(sce.raw, 
                                 subsets=list(ERCC=ercc_genes, 
                                              Mito=mito_genes,
                                              Ribo=ribo_genes))

saveRDS(sce.raw, "sce_postdblfinder.rds")

##### MANUAL QC FOR CELLS #####

sce.qc <- sce.raw

# No library size correction because the original authors don't do it

# Feature detection correction
discard_ngenes <- (sce.qc$detected < 200)
sce.qc$discard_ngenes <- discard_ngenes
print("Células eliminadas por número de genes: ")
table(discard_ngenes)

# %Mitochondrial correction
discard_mito <- (sce.qc$subsets_Mito_percent > 5)
sce.qc$discard_mito <- discard_mito
print("Células eliminadas por porcentaje mitocondrial: ")
table(discard_mito)

# Doublets correction
discard_doublet <- (sce.qc$scDblFinder.class == "doublet")
sce.qc$discard_doublet <- discard_doublet
print("Células eliminadas por dobletes: ")
table(discard_doublet)

# Discard the cells
discard <-  (discard_ngenes | discard_mito | discard_doublet)
sce.qc$discard <- discard
coldata <- as.data.frame(colData(sce.qc))
print("Células a descartar: ")
table(discard) # Just based on number of genes and doublets is discarded


## Non-coding genes (https://www.genenames.org/about/guidelines/) ##
# Pseudogenes are identified with the last consonant 'P'
pseudogenes <- rownames(sce.qc)[grep("P$",rownames(sce.qc))]
# Processed pseudogenes
pseudogenes1 <- rownames(sce.qc)[grep("P[0-9]+$",rownames(sce.qc))]
# tRNA
trna <- rownames(sce.qc)[grep("^TR.-",rownames(sce.qc))]
# Small nuclear RNAs
snrna <- rownames(sce.qc)[grep("^RNU",rownames(sce.qc))]
# Small nucleolar RNAs
snorna <- rownames(sce.qc)[grep("^SNORD",rownames(sce.qc))] # SNORD# for “small nucleolar RNA, C/D box” genes
snorna1 <- rownames(sce.qc)[grep("^SNORA",rownames(sce.qc))] # SNORA# for “small nucleolar RNA, H/ACA box” genes
snorna2 <- rownames(sce.qc)[grep("^SCARNA",rownames(sce.qc))] # SCARNA# for “small Cajal body‐specific RNA” genes
# Ribosomal RNAs
rrna <- rownames(sce.qc)[grep("^RNA[0-9]",rownames(sce.qc))]
# (pseudo?)genes not described https://www.biostars.org/p/51456/
rprna <- rownames(sce.qc)[grep("^RP[0-9]",rownames(sce.qc))]

# Discard genes not detected in X cells (not made in the original paper)
# print("Genes a descartar por aparecer en menos de 10 celulas: ")
# selected_features <- scater::nexprs(sce.qc, byrow = TRUE, detection_limit = 1) >= 10
# In this case, the discarded genes are the FALSE ones
# table(selected_features)

# Combine all the genes to discard
genes.descarte <- c(pseudogenes, pseudogenes1, trna, snrna, snorna, snorna1,
                    snorna2, rrna, rprna, mito_genes, ribo_genes, ercc_genes)

# In this study, we see there are no genes to discard by these filters
print("Genes a descartar por pertenecer a mito, ribo, etc: ")
length(genes.descarte) # = 0

# So we keep all genes/features

saveRDS(sce.qc, "processed_data/sce_qc.rds")
##### QC FILTERED DATASET ##### 
sce.filt <- sce.qc[, !discard]

dim(sce.qc)
dim(sce.filt)

saveRDS(sce.filt, "processed_data/sce_filt.rds")

#############################################################
#                   DIAGNOSTIC QC PLOTS                     #
#############################################################

#-------- Comparative plots (before and after in a single plot) ---------------

path <- "./figures/1.QC/"

plot <- scater::plotColData(sce.qc, x="Patient_ID", y="detected",
                            colour_by="discard_ngenes",
                            other_fields = "Clinical.AD.or.not") +
  facet_wrap(~Clinical.AD.or.not, scales = "free_x") +
  labs(title = "Genes por paciente", y = "Genes por célula", x = "Paciente") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_colour_manual(name="", values = c("FALSE" = "gray", "TRUE" = "red"))
ggsave(filename=paste0(path, "ComparativeGenes.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path, "ComparativeGenes.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")


# This plot doesn't make sense this time because no cells were discarded
# by their %mitochondrial
plot <- scater::plotColData(sce.qc, x="Patient_ID", y="subsets_Mito_percent",
                            colour_by="discard_mito",
                            other_fields = "Clinical.AD.or.not") +
  facet_wrap(~Clinical.AD.or.not, scales = "free_x") +
  labs(title = "% Mitocondrial por paciente", y = "%Mitocondrial por célula", x = "Paciente") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_colour_manual(name="", values = c("FALSE" = "gray", "TRUE" = "red"))
ggsave(filename=paste0(path, "ComparativeMitochondrial.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path, "ComparativeMitochondrial.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")


#------------------------- Before QC plots --------------------------------

plot <- scater::plotColData(sce.qc, x="Patient_ID", y="sum", 
                            other_fields = "Clinical.AD.or.not",
                            colour_by = "Group_sex") + 
  facet_wrap(~Clinical.AD.or.not, scales = "free_x") +
  labs(title = "Conteos por paciente", y = "Conteos por célula", x = "Paciente") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_viridis_d(name="Condición y sexo", labels=c("AD_F"="AD_Mujer",
                                                          "AD_M"="AD_Hombre",
                                                          "Control_F"="Control_Mujer",
                                                          "Control_M"="Control_Hombre"))
ggsave(filename=paste0(path, "BeforeQC_Counts.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path, "BeforeQC_Counts.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")


plot <- scater::plotColData(sce.qc, x="Patient_ID", y="detected",
                            colour_by = "Group_sex",
                            other_fields = "Clinical.AD.or.not") +
  facet_wrap(~Clinical.AD.or.not, scales = "free_x") +
  labs(title = "Genes por paciente antes del filtrado", y = "Genes por célula", x = "Paciente") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_colour_viridis_d(name="Condición y sexo", labels=c("AD_F"="AD_Mujer",
                                                           "AD_M"="AD_Hombre",
                                                           "Control_F"="Control_Mujer",
                                                           "Control_M"="Control_Hombre"))
ggsave(filename=paste0(path, "BeforeQC_Genes.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path, "BeforeQC_Genes.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")


plot <- scater::plotColData(sce.qc, x="Patient_ID", y="subsets_Mito_percent",
                            colour_by = "Group_sex",
                            other_fields = "Clinical.AD.or.not") +
  facet_wrap(~Clinical.AD.or.not, scales = "free_x") +
  labs(title = "% Mito por paciente antes del filtrado", y = "% Mitocondrial", x = "Paciente") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_colour_viridis_d(name="Condición y sexo", labels=c("AD_F"="AD_Mujer",
                                                           "AD_M"="AD_Hombre",
                                                           "Control_F"="Control_Mujer",
                                                           "Control_M"="Control_Hombre"))
ggsave(filename=paste0(path, "BeforeQC_Mitochondrial.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path, "BeforeQC_Mitochondrial.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")


plot <- ggplot(as.data.frame(colData(sce.qc)), aes(x = sum)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Distribución de la densidad de conteos (antes del QC)",
       x = "Conteos totales (UMIs)", y = "Densidad") +
  theme_minimal()
ggsave(filename=paste0(path, "BeforeQC_Counts_DensityDistrib.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path, "BeforeQC_Counts_DensityDistrib.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

#------------------------- After QC plots --------------------------------
# To plot the cell discard criteria
sce.filtdetected <- sce.qc[, !discard_ngenes]

plot <- scater::plotColData(sce.filtdetected, x="Patient_ID", y="detected",
                            colour_by = "Group_sex",
                            other_fields = "Clinical.AD.or.not") +
  facet_wrap(~Clinical.AD.or.not, scales = "free_x") +
  labs(title = "Genes por paciente tras el filtrado", y = "Genes por célula", x = "Paciente") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_viridis_d(name="Condición y sexo", labels=c("AD_F"="AD_Mujer",
                                                            "AD_M"="AD_Hombre",
                                                            "Control_F"="Control_Mujer",
                                                            "Control_M"="Control_Hombre"))
ggsave(filename=paste0(path, "AfterQC_Genes.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path, "AfterQC_Genes.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

#counts_postQC <- Matrix::colSums(counts(sce.filt)) # equivalent to 'sum'
plot <- ggplot(as.data.frame(colData(sce.filt)), aes(x = sum)) +
    geom_density(fill = "orange", alpha = 0.5) +
    labs(title = "Distribución de la densidad de conteos (después del QC)",
         x = "Conteos totales (UMIs)", y = "Densidad") +
    theme_minimal()
ggsave(filename=paste0(path, "AfterQC_Counts_DensityDistrib.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path, "AfterQC_Counts_DensityDistrib.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

#-------------------------- Extra plots -----------------------------------

## Count density distribution comparison before and after QC
# Extract counts before and after QC
counts_preQC <- sce.qc$sum
counts_postQC <- sce.filt$sum

# Create a combined df
df <- data.frame(
  Counts = c(counts_preQC, counts_postQC),
  QC_Status = rep(c("Pre-QC", "Post-QC"), 
                  times = c(length(counts_preQC), length(counts_postQC)))
  )

plot <- ggplot(df, aes(x = Counts, fill = QC_Status)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Pre-QC" = "skyblue", "Post-QC" = "orange")) +
  labs(
    title = "Distribución de conteos antes y después del QC",
    x = "Conteos totales (UMIs)",
    y = "Densidad",
    fill = "Fase"
  ) +
  theme_minimal()
ggsave(filename=paste0(path, "Density_Comparison_Pre_vs_Post_QC.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path, "Density_Comparison_Pre_vs_Post_QC.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")


# Number of total discarded cells per sample
coldata.qc <- as.data.frame(colData(sce.qc))
grouped_df <- coldata.qc %>%
  group_by(Patient_ID, Group_sex) %>%
  summarise(Num_cels_discarded = sum(discard))

plot <- ggplot(grouped_df, aes(x = as.factor(Patient_ID), y = Num_cels_discarded)) +
    geom_col(aes(fill=Group_sex)) +
    labs(title = "Células descartadas por paciente", y = "Número de células", x = "Paciente") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    scale_fill_viridis_d(name="Condición y sexo", labels=c("AD_F"="AD_Mujer",
                                                           "AD_M"="AD_Hombre",
                                                           "Control_F"="Control_Mujer",
                                                           "Control_M"="Control_Hombre"))
ggsave(filename=paste0(path, "NumCellsDiscardedPerPatient.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path, "NumCellsDiscardedPerPatient.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

# Proportion of discarded cells per patient
grouped_df2 <- as.data.frame(colData(sce.qc)) %>%
  group_by(Patient_ID) %>%
  summarise(Totales = n(),
            Descartadas = sum(discard)
  ) %>%
  pivot_longer(cols = c(Totales, Descartadas),
               names_to = "Tipo", 
               values_to = "Cantidad")

plot <- ggplot(grouped_df2, aes(x = Patient_ID, y = Cantidad)) +
    geom_col(aes(fill=Tipo)) +
    labs(title = "Células por paciente", y = "Número de células", x = "Paciente") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.text = element_text(size = 7)) +
    scale_fill_manual(values = c("Totales" = "pink", "Descartadas" = "red"))
ggsave(filename=paste0(path, "ProportionDiscardedCellsPerPatient.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path, "ProportionDiscardedCellsPerPatient.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")


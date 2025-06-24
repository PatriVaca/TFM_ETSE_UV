# Single cell RNA-seq data analysis
# 4.Dimension Reduction
# Author: Patricia del Carmen Vaca Rubio
# Date: 05/25

# After analyzing 3.FeatureSelection script output plots, we will continue
# our analysis just with the HVGs selected via modelGeneVar with the
# density.weights=FALSE option

# Load packages
library(scran)
library(SingleCellExperiment)
library(scater)

setwd("/clinicfs/projects/i63/tfm_hipocampo/syn52293442_hc/")
system("mkdir -p figures/4.DimensionReduction/4.2.ModelGeneVar_dw")

sce.featuresel <- readRDS("processed_data/sce_featuresel.rds")

#############################################################################
# HVGs obtained via modelGeneVar non-overfitted (density.weigths = FALSE)   #
#############################################################################

path2 <- "./figures/4.DimensionReduction/4.2.ModelGeneVar_dw/"

# ------------------------------- PCA ----------------------------------------
set.seed(100)
scegv_dw <- fixedPCA(sce.featuresel,
                     subset.row = rowSubset(sce.featuresel, "HVG.genevar.weigths")) 
percent.vargv_dw <- attr(reducedDim(scegv_dw), "percentVar")
cum.sumgv_dw <- cumsum(percent.vargv_dw) 
elbowgv_dw <- PCAtools::findElbowPoint(percent.vargv_dw)
print(paste("Número de PCs elegidos por método del codo:", elbowgv_dw))

# ----------- PCA plots -----------
jpeg(file=paste0(path2, "PCA1.jpeg", sep=""),
     width = 7, height = 5, units = "in", res = 300)
plot(percent.vargv_dw, xlab="Componente Principal", ylab="Varianza explicada (%)")
abline(v=elbowgv_dw, col="red")
text(elbowgv_dw, (elbowgv_dw-0.8), labels = paste("Codo: ", elbowgv_dw, sep = ""),
     cex = 1, pos = 4, col = "red")
dev.off()

jpeg(file=paste0(path2, "PCA2.jpeg", sep=""),
     width = 7, height = 5, units = "in", res = 300)
plot(percent.vargv_dw, log="y", xlab="Componente Principal", ylab="Log - Varianza explicada (%)")
abline(v=elbowgv_dw, col="red")
text(elbowgv_dw, (elbowgv_dw-5), labels = paste("Codo: ", elbowgv_dw, sep = ""),
     cex = 1, pos = 4, col = "red")
dev.off()

jpeg(file=paste0(path2, "PCA3.jpeg", sep=""),
     width = 7, height = 5, units = "in", res = 300)
plot(cum.sumgv_dw, xlab="PC", ylab="Varianza acumulada explicada (%)")
abline(v=elbowgv_dw, col="red")
text(elbowgv_dw, (cum.sumgv_dw[elbowgv_dw]-0.8),
     labels = paste("Varianza total expl.: ",
                    format(cum.sumgv_dw[elbowgv_dw], 2),"%", sep=" "),
     cex = 0.8, pos = 4, col = "red")
dev.off()

# Keep only the PCs before the elbow
reducedDim(scegv_dw, "PCA") <- reducedDim(scegv_dw, "PCA")[, 1:elbowgv_dw]


plot <- plotReducedDim(scegv_dw, dimred="PCA", colour_by="individualID") +
  scale_color_viridis_d("Paciente")
ggsave(filename=paste0(path2, "PCAplotpersample.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "PCAplotpersample.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")


plot <- plotReducedDim(scegv_dw, dimred="PCA", colour_by="Group_sex") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre",
                                   "Control_Mujer", "Control_Hombre"))
ggsave(filename=paste0(path2, "PCAplotperCondSex.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "PCAplotperCondSex.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")


plot <- plotReducedDim(scegv_dw, dimred="PCA", ncomponents=4,
                       colour_by="Group_sex") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre",
                                   "Control_Mujer", "Control_Hombre"))
ggsave(filename=paste0(path2, "MultiplePCAgraph.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "MultiplePCAgraph.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

# ------------------------------ t-SNE plots ---------------------------------
# Seed 100 - Perplexity 5
set.seed(100)
scegv_dw <- runTSNE(scegv_dw, dimred = "PCA", perplexity = 5, name = "TSNE5")
out5scegv_dw <- plotReducedDim(scegv_dw, dimred="TSNE5", colour_by="Group_sex") +
  ggtitle("t-SNE perplexity = 5") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre",
                                   "Control_Mujer", "Control_Hombre")) +
  labs(x = "TSNE1", y = "TSNE2")

# Seed 100 - Perplexity 80
set.seed(100)
scegv_dw <- runTSNE(scegv_dw, dimred = "PCA", perplexity = 80, name = "TSNE80")
out80scegv_dw <- plotReducedDim(scegv_dw, dimred="TSNE80", colour_by="Group_sex") +
  ggtitle("t-SNE perplexity = 80") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre",
                                   "Control_Mujer", "Control_Hombre")) +
  labs(x = "TSNE1", y = "TSNE2")

# Seed 100 - Perplexity 120
set.seed(100)
scegv_dw <- runTSNE(scegv_dw, dimred = "PCA", perplexity = 120, name = "TSNE120")
out120scegv_dw <- plotReducedDim(scegv_dw, dimred="TSNE120", colour_by="Group_sex") +
  ggtitle("t-SNE perplexity = 120") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre",
                                   "Control_Mujer", "Control_Hombre")) +
  labs(x = "TSNE1", y = "TSNE2")

# Seed 100 - Perplexity 200
set.seed(100)
scegv_dw <- runTSNE(scegv_dw, dimred = "PCA", perplexity = 200, name = "TSNE200")
out200scegv_dw <- plotReducedDim(scegv_dw, dimred="TSNE200", colour_by="Group_sex") +
  ggtitle("t-SNE perplexity = 200") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre",
                                   "Control_Mujer", "Control_Hombre")) +
  labs(x = "TSNE1", y = "TSNE2")


# -------------------------------- UMAP plots ---------------------------------
# Seed 100 - neighbor 5
set.seed(100)
scegv_dw <- runUMAP(scegv_dw, dimred = "PCA", n_neighbors = 5, name = "UMAP5")
umap5scegv_dw <- plotReducedDim(scegv_dw, dimred="UMAP5", colour_by="Group_sex") +
  ggtitle("UMAP n_neighbors = 5") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre",
                                   "Control_Mujer", "Control_Hombre")) +
  labs(x = "UMAP1", y = "UMAP2")

# Seed 100 - neighbor 80
set.seed(100)
scegv_dw <- runUMAP(scegv_dw, dimred = "PCA", n_neighbors = 80, name = "UMAP80")
umap80scegv_dw <- plotReducedDim(scegv_dw, dimred="UMAP80", colour_by="Group_sex") +
  ggtitle("UMAP n_neighbors = 80") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre",
                                   "Control_Mujer", "Control_Hombre")) +
  labs(x = "UMAP1", y = "UMAP2")

# Seed 100 - neighbor 300
set.seed(100)
scegv_dw <- runUMAP(scegv_dw, dimred = "PCA", n_neighbors = 300, name = "UMAP300")
umap300scegv_dw <- plotReducedDim(scegv_dw, dimred="UMAP300", colour_by="Group_sex") +
  ggtitle("UMAP n_neighbors = 300") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre",
                                   "Control_Mujer", "Control_Hombre")) +
  labs(x = "UMAP1", y = "UMAP2")

# Seed 100 - neighbor 600
set.seed(100)
scegv_dw <- runUMAP(scegv_dw, dimred = "PCA", n_neighbors = 600, name = "UMAP600")
umap600scegv_dw <- plotReducedDim(scegv_dw, dimred="UMAP600", colour_by="Group_sex") +
  ggtitle("UMAP n_neighbors = 600") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre",
                                   "Control_Mujer", "Control_Hombre")) +
  labs(x = "UMAP1", y = "UMAP2")

saveRDS(scegv_dw, "processed_data/scegv_dw.rds")

######################
### SAVE THE PLOTS ###
######################

########## JPEG ##########
#--------GV DW=FALSE---------
# Save t-SNE gv_dw
jpeg(file=paste0(path2, "tSNE_grid.jpeg", sep=""),
     width = 7, height = 5, units = "in", res = 300)
gridExtra::grid.arrange(out5scegv_dw, out80scegv_dw, out120scegv_dw,
                        out200scegv_dw, nrow=2, ncol=2)
dev.off()

# Save UMAP gv_dw
jpeg(file=paste0(path2, "UMAP_grid.jpeg", sep=""),
     width = 7, height = 5, units = "in", res = 300)
gridExtra::grid.arrange(umap5scegv_dw, umap80scegv_dw, umap300scegv_dw,
                        umap600scegv_dw, nrow=2, ncol=2)
dev.off()

########## SVG ##########
#--------GV DW=FALSE---------
# Save t-SNE gv_dw
svg(file=paste0(path2, "tSNE_grid.svg", sep=""),
    width = 7, height = 5)
gridExtra::grid.arrange(out5scegv_dw, out80scegv_dw, out120scegv_dw,
                        out200scegv_dw, nrow=2, ncol=2)
dev.off()

# Save UMAP gv_dw
svg(file=paste0(path2, "UMAP_grid.svg", sep=""),
    width = 7, height = 5)
gridExtra::grid.arrange(umap5scegv_dw, umap80scegv_dw, umap300scegv_dw,
                        umap600scegv_dw, nrow=2, ncol=2)
dev.off()



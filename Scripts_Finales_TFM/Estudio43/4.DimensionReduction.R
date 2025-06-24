# Single cell RNA-seq data analysis
# 4.Dimension Reduction
# Author: Patricia del Carmen Vaca Rubio
# Date: 04/25

setwd("/clinicfs/projects/i63/tfm_hipocampo/Study43/")
system("mkdir -p figures/4.DimensionReduction")
system("mkdir -p figures/4.DimensionReduction/4.1.ModelGeneVar")
system("mkdir -p figures/4.DimensionReduction/4.2.ModelGeneVar_dw")
system("mkdir -p figures/4.DimensionReduction/4.3.ModelGeneVarByPoisson")

library(scran)
library(SingleCellExperiment)
library(scater)

sce.featuresel <- readRDS("processed_data/sce_featuresel.rds")


##############################################################
# HVGs obtained via modelGeneVar
##############################################################

path1 <- "./figures/4.DimensionReduction/4.1.ModelGeneVar/"

# ------------------------------- PCA ----------------------------------------
set.seed(100)
scegv <- fixedPCA(sce.featuresel,
                  subset.row = rowSubset(sce.featuresel, "HVGs.genevar"))
percent.var <- attr(reducedDim(scegv), "percentVar")
cum.sum <- cumsum(percent.var) 
elbow <- PCAtools::findElbowPoint(percent.var)
print(paste("Número de PCs elegidos por método del codo:", elbow))

# ----------- PCA plots -----------
jpeg(file=paste0(path1, "PCA1.jpeg", sep=""),
     width = 7, height = 5, units = "in", res = 300)
plot(percent.var, xlab="Componente Principal", ylab="Varianza explicada (%)")
abline(v=elbow, col="blue")
text(elbow, (elbow-0.8), labels = paste("Codo: ", elbow, sep = ""),
     cex = 1, pos = 4, col = "blue")
dev.off()

jpeg(file=paste0(path1, "PCA2.jpeg", sep=""),
     width = 7, height = 5, units = "in", res = 300)
plot(percent.var, log="y", xlab="Componente Principal", ylab="Log - Varianza explicada (%)")
abline(v=elbow, col="blue")# Log representation for visualization purposes
text(elbow, (elbow-5), labels = paste("Codo: ", elbow, sep = ""),
     cex = 1, pos = 4, col = "blue")
dev.off()

jpeg(file=paste0(path1, "PCA3.jpeg", sep=""),
     width = 7, height = 5, units = "in", res = 300)
plot(cum.sum, xlab="PC", ylab="Varianza acumulada explicada (%)")
abline(v=elbow, col="blue")
text(elbow, (cum.sum[elbow]-0.8),
     labels = paste("Varianza total expl.: ",
                    format(cum.sum[elbow], 2),"%", sep=" "),
                    cex = 0.8, pos = 4, col = "blue")
dev.off()

# Keep only the PCs before the elbow
reducedDim(scegv, "PCA") <- reducedDim(scegv, "PCA")[, 1:elbow]

plot <- plotReducedDim(scegv, dimred="PCA", colour_by="Sample") +
  scale_color_viridis_d("Muestra")
ggsave(filename=paste0(path1, "PCAplotpersample.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path1, "PCAplotpersample.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")


plot <- plotReducedDim(scegv, dimred="PCA", colour_by="Group_sex") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre",
                                   "Control_Mujer", "Control_Hombre"))
ggsave(filename=paste0(path1, "PCAplotperCondSex.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path1, "PCAplotperCondSex.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")


plot <- plotReducedDim(scegv, dimred="PCA", ncomponents=4,
               colour_by="Group_sex") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre",
                                   "Control_Mujer", "Control_Hombre"))
ggsave(filename=paste0(path1, "MultiplePCAgraph.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path1, "MultiplePCAgraph.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")
  

# ------------------------------ t-SNE plots ---------------------------------
## Try different perplexities for different granularity in the visualization ##
# Seed 100 - Perplexity 5
set.seed(100)
scegv <- runTSNE(scegv, dimred = "PCA", perplexity = 5, name = "TSNE5")
out5 <- plotReducedDim(scegv, dimred="TSNE5", colour_by="Group_sex") +
  ggtitle("t-SNE perplexity = 5") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre",
                                   "Control_Mujer", "Control_Hombre")) +
  labs(x = "TSNE1", y = "TSNE2")


# Seed 100 - Perplexity 80
set.seed(100)
scegv <- runTSNE(scegv, dimred = "PCA", perplexity = 80, name = "TSNE80")
out80 <- plotReducedDim(scegv, dimred="TSNE80", colour_by="Group_sex") +
  ggtitle("t-SNE perplexity = 80") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre", "Control_Mujer", "Control_Hombre")) +
  labs(x = "TSNE1", y = "TSNE2")


# Seed 100 - Perplexity 120
set.seed(100)
scegv <- runTSNE(scegv, dimred = "PCA", perplexity = 120, name = "TSNE120")
out120 <- plotReducedDim(scegv, dimred="TSNE120", colour_by="Group_sex") +
  ggtitle("t-SNE perplexity = 120") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre", "Control_Mujer", "Control_Hombre")) +
  labs(x = "TSNE1", y = "TSNE2")


# Seed 100 - Perplexity 200
set.seed(100)
scegv <- runTSNE(scegv, dimred = "PCA", perplexity = 200, name = "TSNE200")
out200 <- plotReducedDim(scegv, dimred="TSNE200", colour_by="Group_sex") +
  ggtitle("t-SNE perplexity = 200") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre", "Control_Mujer", "Control_Hombre")) +
  labs(x = "TSNE1", y = "TSNE2")


# -------------------------------- UMAP plots ---------------------------------
# Seed 100 - neighbor 5
set.seed(100)
scegv <- runUMAP(scegv, dimred = "PCA", n_neighbors = 5, name = "UMAP5")
umap5 <- plotReducedDim(scegv, dimred="UMAP5", colour_by="Group_sex") +
  ggtitle("UMAP n_neighbors = 5") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre", "Control_Mujer", "Control_Hombre")) +
  labs(x = "UMAP1", y = "UMAP2")


# Seed 100 - neighbor 80
set.seed(100)
scegv <- runUMAP(scegv, dimred = "PCA", n_neighbors = 80, name = "UMAP80")
umap80 <- plotReducedDim(scegv, dimred="UMAP80", colour_by="Group_sex") +
  ggtitle("UMAP n_neighbors = 80") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre", "Control_Mujer", "Control_Hombre")) +
  labs(x = "UMAP1", y = "UMAP2")


# Seed 100 - neighbor 300
set.seed(100)
scegv <- runUMAP(scegv, dimred = "PCA", n_neighbors = 300, name = "UMAP300")
umap300 <- plotReducedDim(scegv, dimred="UMAP300", colour_by="Group_sex") +
  ggtitle("UMAP n_neighbors = 300") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre", "Control_Mujer", "Control_Hombre")) +
  labs(x = "UMAP1", y = "UMAP2")


# Seed 100 - neighbor 600
set.seed(100)
scegv <- runUMAP(scegv, dimred = "PCA", n_neighbors = 600, name = "UMAP600")
umap600 <- plotReducedDim(scegv, dimred="UMAP600", colour_by="Group_sex") +
  ggtitle("UMAP n_neighbors = 600") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre", "Control_Mujer", "Control_Hombre")) +
  labs(x = "UMAP1", y = "UMAP2")

saveRDS(scegv, "processed_data/scegv.rds")


##############################################################
# HVGs obtained via modelGeneVar non-overfitted (density.weigths = FALSE)
##############################################################

path2 <- "./figures/4.DimensionReduction/4.2.ModelGeneVar_dw/"

set.seed(100)
# Hacemos la PCA con los HVGs de genevar
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

plot <- plotReducedDim(scegv_dw, dimred="PCA", colour_by="Sample") +
  scale_color_viridis_d("Muestra")
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

##############################################################
# HVGs calculados a partir de modelGeneVarByPoisson
##############################################################

path3 <- "./figures/4.DimensionReduction/4.3.ModelGeneVarByPoisson/"
  
# ------------------------------- PCA ----------------------------------------
set.seed(100)
# Hacemos la PCA con los HVGs de genevar
sceP <- fixedPCA(sce.featuresel,
                 subset.row = rowSubset(sce.featuresel, "HVGs.poisson"))
percent.varP <- attr(reducedDim(sceP), "percentVar")
cum.sumP <- cumsum(percent.varP) 
elbowP <- PCAtools::findElbowPoint(percent.varP)
print(paste("Número de PCs elegidos por método del codo:", elbowP))

# ----------- PCA plots -----------

jpeg(file=paste0(path3, "PCA1.jpeg", sep=""),
    width = 7, height = 5, units = "in", res = 300)
plot(percent.varP, xlab="Componente Principal", ylab="Varianza explicada (%)")
abline(v=elbowP, col="green")
text(elbowP, (elbowP-0.8), labels = paste("Codo: ", elbowP, sep = ""),
     cex = 1, pos = 4, col = "green")
dev.off()

jpeg(file=paste0(path3, "PCA2.jpeg", sep=""),
     width = 7, height = 5, units = "in", res = 300)
plot(percent.varP, log="y", xlab="Componente Principal", ylab="Log - Varianza explicada (%)")
abline(v=elbowP, col="green")
text(elbow, (elbowP-5), labels = paste("Codo: ", elbowP, sep = ""),
     cex = 1, pos = 4, col = "green")
dev.off()

jpeg(file=paste0(path3, "PCA3.jpeg", sep=""),
     width = 7, height = 5, units = "in", res = 300)
plot(cum.sumP, xlab="PC", ylab="Varianza acumulada explicada (%)")
abline(v=elbowP, col="green")
text(elbowP, (cum.sumP[elbowP]-0.8),
     labels = paste("Varianza total expl.: ",
                    format(cum.sumP[elbowP], 2),"%", sep=" "),
     cex = 0.8, pos = 4, col = "green")
dev.off()

# Keep only the PCs before the elbow
reducedDim(sceP, "PCA") <- reducedDim(sceP, "PCA")[, 1:elbowP]

plot <- plotReducedDim(sceP, dimred="PCA", colour_by="Sample") +
  scale_color_viridis_d("Muestra")
ggsave(filename=paste0(path3, "PCAplotpersample.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path3, "PCAplotpersample.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

plot <- plotReducedDim(sceP, dimred="PCA", colour_by="Group_sex") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre",
                                   "Control_Mujer", "Control_Hombre"))
ggsave(filename=paste0(path3, "PCAplotperCondSex.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path3, "PCAplotperCondSex.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

plot <- plotReducedDim(sceP, dimred="PCA", ncomponents=4,
               colour_by="Group_sex") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre", 
                                   "Control_Mujer", "Control_Hombre"))
ggsave(filename=paste0(path3, "MultiplePCAgraph.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path3, "MultiplePCAgraph.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

# ------------------------------ t-SNE plots ---------------------------------
# Seed 100 - Perplexity 5
set.seed(100)
sceP <- runTSNE(sceP, dimred = "PCA", perplexity = 5, name = "TSNE5")
out5P <- plotReducedDim(sceP, dimred="TSNE5", colour_by="Group_sex") +
  ggtitle("t-SNE perplexity = 5") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre",
                                   "Control_Mujer", "Control_Hombre")) +
  labs(x = "TSNE1", y = "TSNE2")

# Seed 100 - Perplexity 80
set.seed(100)
sceP <- runTSNE(sceP, dimred = "PCA", perplexity = 80, name = "TSNE80")
out80P <- plotReducedDim(sceP, dimred="TSNE80", colour_by="Group_sex") +
  ggtitle("t-SNE perplexity = 80") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre",
                                   "Control_Mujer", "Control_Hombre")) +
  labs(x = "TSNE1", y = "TSNE2")

# Seed 100 - Perplexity 120
set.seed(100)
sceP <- runTSNE(sceP, dimred = "PCA", perplexity = 120, name = "TSNE120")
out120P <- plotReducedDim(sceP, dimred="TSNE120", colour_by="Group_sex") +
  ggtitle("t-SNE perplexity = 120") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre",
                                   "Control_Mujer", "Control_Hombre")) +
  labs(x = "TSNE1", y = "TSNE2")

# Seed 100 - Perplexity 200
set.seed(100)
sceP <- runTSNE(sceP, dimred = "PCA", perplexity = 200, name = "TSNE200")
out200P <- plotReducedDim(sceP, dimred="TSNE200", colour_by="Group_sex") +
  ggtitle("t-SNE perplexity = 200") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre",
                                   "Control_Mujer", "Control_Hombre")) +
  labs(x = "TSNE1", y = "TSNE2")


# -------------------------------- UMAP plots ---------------------------------
# Seed 100 - neighbor 5
set.seed(100)
sceP <- runUMAP(sceP, dimred = "PCA", n_neighbors = 5, name = "UMAP5")
umap5P <- plotReducedDim(sceP, dimred="UMAP5", colour_by="Group_sex") +
  ggtitle("UMAP n_neighbors = 5") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre",
                                   "Control_Mujer", "Control_Hombre")) +
  labs(x = "UMAP1", y = "UMAP2")

# Seed 100 - neighbor 80
set.seed(100)
sceP <- runUMAP(sceP, dimred = "PCA", n_neighbors = 80, name = "UMAP80")
umap80P <- plotReducedDim(sceP, dimred="UMAP80", colour_by="Group_sex") +
  ggtitle("UMAP n_neighbors = 80") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre",
                                   "Control_Mujer", "Control_Hombre")) +
  labs(x = "UMAP1", y = "UMAP2")

# Seed 100 - neighbor 300
set.seed(100)
sceP <- runUMAP(sceP, dimred = "PCA", n_neighbors = 300, name = "UMAP300")
umap300P <- plotReducedDim(sceP, dimred="UMAP300", colour_by="Group_sex") +
  ggtitle("UMAP n_neighbors = 300") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre",
                                   "Control_Mujer", "Control_Hombre")) +
  labs(x = "UMAP1", y = "UMAP2")

# Seed 100 - neighbor 600
set.seed(100)
sceP <- runUMAP(sceP, dimred = "PCA", n_neighbors = 600, name = "UMAP600")
umap600P <- plotReducedDim(sceP, dimred="UMAP600", colour_by="Group_sex") +
  ggtitle("UMAP n_neighbors = 600") +
  scale_color_viridis_d(name = "Condición y sexo",
                        labels = c("AD_Mujer", "AD_Hombre",
                                   "Control_Mujer", "Control_Hombre")) +
  labs(x = "UMAP1", y = "UMAP2")

saveRDS(sceP, "processed_data/sceP.rds")

######################
### SAVE THE PLOTS ###
######################

########## JPEG ##########
#--------GV---------
# Save t-SNE gv
jpeg(file=paste0(path1, "tSNE_grid.jpeg", sep=""),
     width = 7, height = 5, units = "in", res = 300)
gridExtra::grid.arrange(out5, out80, out120, out200, nrow=2, ncol=2)
dev.off()

# Save UMAP gv
jpeg(file=paste0(path1, "UMAP_grid.jpeg", sep=""),
     width = 7, height = 5, units = "in", res = 300)
gridExtra::grid.arrange(umap5, umap80, umap300, umap600, nrow=2, ncol=2)
dev.off()


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

#--------Poisson---------
# Save t-SNE poisson
jpeg(file=paste0(path3, "tSNE_grid.jpeg", sep=""),
     width = 7, height = 5, units = "in", res = 300)
gridExtra::grid.arrange(out5P, out80P, out120P, out200P, nrow=2, ncol=2)
dev.off()

# Save UMAP poisson
jpeg(file=paste0(path3, "UMAP_grid.jpeg", sep=""),
     width = 7, height = 5, units = "in", res = 300)
gridExtra::grid.arrange(umap5P, umap80P, umap300P, umap600P, nrow=2, ncol=2)
dev.off()

########## SVG ##########
#--------GV---------
# Save t-SNE gv
svg(file=paste0(path1, "tSNE_grid.svg", sep=""),
    width = 7, height = 5)
gridExtra::grid.arrange(out5scegv_dw, out80scegv_dw, out120scegv_dw,
                        out200scegv_dw, nrow=2, ncol=2)
dev.off()

# Save UMAP gv_dw
svg(file=paste0(path1, "UMAP_grid.svg", sep=""),
    width = 7, height = 5)
gridExtra::grid.arrange(umap5scegv_dw, umap80scegv_dw, umap300scegv_dw,
                        umap600scegv_dw, nrow=2, ncol=2)
dev.off()

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

#--------Poisson---------
# Save t-SNE gv_dw
svg(file=paste0(path3, "tSNE_grid.svg", sep=""),
    width = 7, height = 5)
gridExtra::grid.arrange(out5scegv_dw, out80scegv_dw, out120scegv_dw,
                        out200scegv_dw, nrow=2, ncol=2)
dev.off()

# Save UMAP gv_dw
svg(file=paste0(path3, "UMAP_grid.svg", sep=""),
    width = 7, height = 5)
gridExtra::grid.arrange(umap5scegv_dw, umap80scegv_dw, umap300scegv_dw,
                        umap600scegv_dw, nrow=2, ncol=2)
dev.off()

# Single cell RNA-seq data analysis
# 2.Normalization
# Author: Patricia del Carmen Vaca Rubio
# Date: 03/25

setwd("/clinicfs/projects/i63/tfm_hipocampo/Study43/")
system("mkdir -p figures/2.Normalization")

library(SingleCellExperiment)
library(scran)
library(ggplot2)

# Load dataset
sce.filt <- readRDS("processed_data/sce_filt.rds")

path <- "./figures/2.Normalization/"

# ----------------------- Normalization methods -------------------------------

### LIBRARY SIZE FACTORS ###
print("Resumen de los size factors de la normalización por tamaño de la  librería: ")
# Calculate the library size factor for each cell
lib.sce.filt <- scuttle::librarySizeFactors(sce.filt) 
sce.filt@colData$library_size_factors <- lib.sce.filt
summary(lib.sce.filt)

jpeg(file=paste0(path, "LibSizeFactorsDistribution.jpeg", sep=""),
     width = 7, height = 5, units = "in", res = 300)
hist(log10(lib.sce.filt), xlab="Log10[Size factor]", col='grey80')
dev.off()

svg(file=paste0(path, "LibSizeFactorsDistribution.svg", sep=""),
    width = 7, height = 5)
hist(log10(lib.sce.filt), xlab="Log10[Size factor]", col='grey80')
dev.off()


##### DECONVOLUTION METHOD ##### 
set.seed(100)
quick_clusters <- scran::quickCluster(sce.filt,
                                      method = "igraph")

# Calculate the SizeFactors for each cell
sce.filt <- scran::computeSumFactors(sce.filt, 
                                     clusters = quick_clusters)

sce.norm <- scuttle::logNormCounts(sce.filt, 
                                   assay.type = "counts", 
                                   log=TRUE, 
                                   pseudo.count=1)

print("Resumen de la normalización por deconvolución")
summary(sizeFactors(sce.norm))

# Save the normalized data
saveRDS(sce.norm, file = "processed_data/sce_norm.rds")

# ------------------------ Normalization plots ------------------------------

colData(sce.filt)$total_filt_counts <- Matrix::colSums(counts(sce.filt))
colData(sce.norm)$total_norm_counts <- Matrix::colSums(logcounts(sce.norm))
sce.filt.df <- as.data.frame(colData(sce.filt))
sce.norm.df <- as.data.frame(colData(sce.norm))

# Before normalization #
plot <- ggplot(sce.filt.df, aes(x = Sample, y = total_filt_counts, fill = Group_sex)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge()) +
  labs(title = "Conteos por paciente antes de la normalización",
       y = "Total de conteos por célula",
       x = "Paciente") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_d(name="Condición y sexo", labels=c("Alzheimer's Disease_Female"="AD_Mujer",
                                                          "Alzheimer's Disease_Male"="AD_Hombre",
                                                          "Control_Female"="Control_Mujer",
                                                          "Control_Male"="Control_Hombre"))
ggsave(filename=paste0(path, "BeforeNorm_Counts.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path, "BeforeNorm_Counts.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

# After normalization #
plot <- ggplot(sce.norm.df, aes(x = Sample, y = total_norm_counts, fill = Group_sex)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge()) +
  labs(title = "Conteos por paciente después de la normalización",
       y = "Total de conteos por célula",
       x = "Paciente") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_d(name="Condición y sexo", labels=c("Alzheimer's Disease_Female"="AD_Mujer",
                                                          "Alzheimer's Disease_Male"="AD_Hombre",
                                                          "Control_Female"="Control_Mujer",
                                                          "Control_Male"="Control_Hombre"))
ggsave(filename=paste0(path, "AfterNorm_Counts.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path, "AfterNorm_Counts.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")


## Library size vs number of genes per cell ##
# Before Normalization #
plot <- ggplot(sce.filt.df, aes(x = detected, y = total_filt_counts)) +
  geom_point(alpha = 0.4, color = "steelblue") +
  labs(x = "Número de genes detectados por célula",
       y = "Tamaño de la librería",
       title = "Antes de la normalización") +
  theme_minimal()
ggsave(filename=paste0(path, "BeforeNorm_LibSizeVSNumGenes.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path, "BeforeNorm_LibSizeVSNumGenes.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

# After Normalization #
plot <- ggplot(sce.norm.df, aes(x = detected, y = total_norm_counts)) +
  geom_point(alpha = 0.4, color = "steelblue") +
  labs(x = "Número de genes detectados por célula",
       y = "Tamaño de la librería",
       title = "Después de la normalización") +
  theme_minimal()
ggsave(filename=paste0(path, "AfterNorm_LibSizeVSNumGenes.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path, "AfterNorm_LibSizeVSNumGenes.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")


## Normalization methods comparison ##
plot <- ggplot(sce.filt.df, aes(x=library_size_factors, y=sizeFactor)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, color="red") +
  theme_classic() +
  xlab("Factor de normalización por tamaño de librería")+
  ylab("Factor de normalización por deconvolución")
ggsave(filename=paste0(path, "NormMethodsComparison.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path, "NormMethodsComparison.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

## LibrarySizeFactor vs Counts ##
plot <- ggplot(sce.filt.df, aes(x=total_filt_counts, y=library_size_factors)) +
  geom_point() +
  theme_classic() +
  xlab("Número de conteos") +
  ylab("Factor de normalización por tamaño de librería")
ggsave(filename=paste0(path, "LibSizeFactorsVSCounts.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path, "LibSizeFactorsVSCounts.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

## SizeFactors (deconvolution method) vs Counts ##
plot <- ggplot(sce.filt.df, aes(x=total_filt_counts, y=sizeFactor)) +
  geom_point() +
  theme_classic() +
  xlab("Número de conteos") +
  ylab("Factor de normalización por deconvolución")
ggsave(filename=paste0(path, "SizeFactorsVSCounts.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path, "SizeFactorsVSCounts.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")


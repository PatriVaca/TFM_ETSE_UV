# Single cell RNA-seq data analysis
# 3.Feature selection (selection of highly variable genes, HVGs)
# Author: Patricia del Carmen Vaca Rubio
# Date: 04/25

setwd("/clinicfs/projects/i63/tfm_hipocampo/Study43/")
system("mkdir -p figures/3.FeatureSelection")

library(scater)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)

sce.norm <- readRDS("processed_data/sce_norm.rds")

path <- "./figures/3.FeatureSelection/"

# ------------- ModelGeneVar() without overfitting correction -----------------
set.seed(100)
dec.genevar <- modelGeneVar(sce.norm)

# Visualizing the fit
plot_data <- data.frame(
  mean = metadata(dec.genevar)$mean,
  var = metadata(dec.genevar)$var
)

plot <- ggplot(plot_data, aes(x = mean, y = var)) +
  geom_point() +
  stat_function(fun = metadata(dec.genevar)$trend, 
                color = "dodgerblue", size = 1.5) +
  labs(x = "Media de la expresión logarítmica",
       y = "Varianza de la expresión logarítmica") +
  theme_minimal()
ggsave(filename=paste0(path, "ModelGeneVar.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path, "ModelGeneVar.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

# Ordering by most interesting genes for inspection
dec.genevar[order(dec.genevar$bio, decreasing=TRUE),]

# ------ ModelGeneVar with overfitting correction (density.weights = FALSE)----
set.seed(100)
dec.genevarweights <- modelGeneVar(sce.norm, density.weights = FALSE)

# Visualizing the fit
plot_data <- data.frame(
  mean = metadata(dec.genevarweights)$mean,
  var = metadata(dec.genevarweights)$var
)

# Create plot
plot <- ggplot(plot_data, aes(x = mean, y = var)) +
  geom_point() +
  stat_function(fun = metadata(dec.genevarweights)$trend, 
                color = "dodgerblue", size = 1.5) +
  labs(x = "Media de la expresión logarítmica",
       y = "Varianza de la expresión logarítmica") +
  theme_minimal()
ggsave(filename=paste0(path, "ModelGeneVar_dw.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path, "ModelGeneVar_dw.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

# Ordering by most interesting genes for inspection.
dec.genevarweights[order(dec.genevarweights$bio, decreasing=TRUE),]

# ----------------------- ModelGeneVarByPoisson() -----------------------------
set.seed(100)
dec.pois <- modelGeneVarByPoisson(sce.norm)

# Visualizing the fit
plot_data <- data.frame(
  mean = metadata(dec.pois)$mean,
  var = metadata(dec.pois)$var
)

# Create plot
plot <- ggplot(plot_data, aes(x = mean, y = var)) +
  geom_point() +
  stat_function(fun = metadata(dec.pois)$trend, 
                color = "dodgerblue", size = 1.5) +
  labs(x = "Media de la expresión logarítmica",
       y = "Varianza de la expresión logarítmica") +
  theme_minimal()
ggsave(filename=paste0(path, "ModelGeneVarByPoisson.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path, "ModelGeneVarByPoisson.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

# Ordering by most interesting genes for inspection
dec.pois <- dec.pois[order(dec.pois$bio, decreasing=TRUE),]
head(dec.pois)

## No batch effect in this studies (just one sample per patient and no batch
# specified in the original metadata)

# -------------------------- Top 2000 HVGs ------------------------------------
## HVGs ModelGeneVar ##
chosen_genevar <- getTopHVGs(dec.genevar, n = 2000)

## HVGs ModelGeneVar density.weights=FALSE ##
chosen_genevar_weights <- getTopHVGs(dec.genevarweights, n = 2000)

## HVGs ModelGeneVarByPoisson ##
chosen_pois <- getTopHVGs(dec.pois, n = 2000)

# -------------------------- Venn Diagram ------------------------------------
# Even though variance vs average plots look different in the three cases,
# let's compare how many HVGs are shared and unique for the three methods
# (HVGs) via a Venn Diagram
myCol <- brewer.pal(3, "Pastel2")

# Chart
venn.diagram(
  x = list(chosen_genevar, chosen_genevar_weights, chosen_pois),
  category.names = c("Overfit", "No overfit", "Poisson"),
  filename = 'figures/3.FeatureSelection/VennDiagram.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1500,
  width = 1500, 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = 1,               
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 1.2,         
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-20, 20, 180),
  cat.dist = c(0.05, 0.05, 0.06),
  cat.fontfamily = "sans",
  rotation = 1
)

rowSubset(sce.norm, "HVGs.genevar") <- chosen_genevar 
rowSubset(sce.norm, "HVGs.poisson") <- chosen_pois 
rowSubset(sce.norm, "HVG.genevar.weigths") <- chosen_genevar_weights 

colnames(rowData(sce.norm))

saveRDS(sce.norm, "processed_data/sce_featuresel.rds")

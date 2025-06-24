# Single cell RNA-seq data analysis
# 6.Cell Type Annotation
# Author: Patricia del Carmen Vaca Rubio
# Date: 05/25

setwd("/clinicfs/projects/i63/tfm_hipocampo/GSE163577/")
# system("mkdir -p figures/6.CellTypeAnnot/6.1.ClusteringGV/6.1.1.SCINA")
system("mkdir -p figures/6.CellTypeAnnot/6.2.ClusteringGV_dw/6.2.1.SCINA")
# system("mkdir -p figures/6.CellTypeAnnot/6.3.ClusteringPoisson/6.3.1.SCINA")

library("dplyr")
library("HGNChelper")
library("SingleCellExperiment")
library("scater")
library("SCINA")
library("SingleR")
library("CHETAH")
library("AnnotationDbi")
library("org.Hs.eg.db")

# Due to row names of the SCE being Ensembl IDs, we need to seek for their
# corresponding Gene Symbols to match the genetic markers list for
# gene annotation. We will just annotate and focus from now on on the HVGs
# obtained via modelGeneVar with the overfitting correction. 

# In this study we will just annotate cell types with SCINA

###############################################################################
#                      ENSEMBLID TO GENE SYMBOL                               #
###############################################################################

# ---------------------------- modelGeneVar -----------------------------------
scegv_clusters <- readRDS("processed_data/scegv_clusters.rds")

# Obtain the Ensembl IDs
ensembl_ids <- rownames(scegv_clusters)

# Map Ensembl IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")  # In case of multiple gene symbols, keep the first one

# Add the gene symbols column
rowData(scegv_clusters)$gene_symbol <- gene_symbols

saveRDS(scegv_clusters, "processed_data/scegv_clusters.rds")

# -------------------- modelGeneVar density.weights = FALSE -------------------

scegv_dw_clusters <- readRDS("processed_data/scegv_dw_clusters.rds")

# Obtain the Ensembl IDs
ensembl_ids <- rownames(scegv_dw_clusters)

# Map Ensembl IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Add the gene symbols column
rowData(scegv_dw_clusters)$gene_symbol <- gene_symbols

saveRDS(scegv_dw_clusters, "processed_data/scegv_dw_clusters.rds")

# ------------------------ modelGeneVarByPoisson ----------------------------

sceP_clusters <- readRDS("processed_data/sceP_clusters.rds")

# Obtain the Ensembl IDs
ensembl_ids <- rownames(sceP_clusters)

# Map Ensembl IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Add the gene symbols column
rowData(sceP_clusters)$gene_symbol <- gene_symbols

saveRDS(sceP_clusters, "processed_data/sceP_clusters.rds")

###############################################################################
#  GENETIC MARKERS (GM) LIST PREPARATION FOR CELL TYPE ANNOTATION WITH SCINA  #
###############################################################################

### 1. INTERSECTION BETWEEN CELLMARKER + SCTYPE GMs WITH THE SCE GENES ###
## 1.1 CELLMARKER (only brain GMs) ##
cellmarkerdb <- read.csv("../notacion_celular/genes_marcadores/cellmarker_experiment.csv")

# Keep only the healthy condition genetic markers (not cancer)
cellmarkerdb <- cellmarkerdb[cellmarkerdb$Cancer=="Normal cell",]

# Keep the tissues of interest (add in this study the vascular-related cells)
cellmarkerdb <- cellmarkerdb[cellmarkerdb$Tissue.Type %in% c("Brain",
                                                             "Hippocampus",
                                                             "Blood vessel",
                                                             "Vessel",
                                                             "Microvascular endothelium"),]

# Get each brain cell types GMs
cellmarker_gmdf <- cellmarkerdb %>%  
  group_by(`Cell.name`) %>%  
  summarise(geneset = list(`Cell.marker`))

# Conversion of df into a list
cellmarker_gmlist = setNames(cellmarker_gmdf$geneset, cellmarker_gmdf$`Cell.name`)

gm.cellmarker <- list(
  "Astrocitos" = c(
    cellmarker_gmlist[["Astrocyte"]],
    cellmarker_gmlist[["Fibrous astrocyte"]],
    cellmarker_gmlist[["Mature astrocyte"]]
  ),
  "Oligodendrocitos" = c(
    cellmarker_gmlist[["Oligodendrocyte"]],
    cellmarker_gmlist[["Mature oligodendrocyte"]],
    cellmarker_gmlist[["Oligodendrocyte-like cell"]]
  ),
  "Microglia" = c(
    cellmarker_gmlist[["Microglial cell"]],
    cellmarker_gmlist[["Homeostatic microglial cell"]],
    cellmarker_gmlist[["Activated microglial cell"]],
    cellmarker_gmlist[["Disease-associated microglial cell"]],
    cellmarker_gmlist[["Demyelinating microglial cell"]],
    cellmarker_gmlist[["Microglia-like cell"]],
    cellmarker_gmlist[["Resting microglial cell"]]
  ),
  "OPC" = c(
    cellmarker_gmlist[["Oligodendrocyte precursor cell"]],
    cellmarker_gmlist[["Oligodendrocyte progenitor cell"]],
    cellmarker_gmlist[["Oligodendrocyte-like cell"]]
  ),
  "Excitadoras" = c(
    cellmarker_gmlist[["Excitatory neuron"]],
    cellmarker_gmlist[["Glutamatergic neuron"]],
    cellmarker_gmlist[["Deep layer neuron"]],
    cellmarker_gmlist[["Upper layer cortical neuron"]]
  ),
  "Inhibidoras" = c(
    cellmarker_gmlist[["Inhibitory neuron"]],
    cellmarker_gmlist[["GABAergic neuron"]],
    cellmarker_gmlist[["Interneuron"]],
    cellmarker_gmlist[["Cholinergic neuron"]],
    cellmarker_gmlist[["Dopaminergic neuron"]],
    cellmarker_gmlist[["Serotonergic neuron"]]
  ),
  "Celulas_epiteliales" = c(
    cellmarker_gmlist[["Endothelial cell"]],
    cellmarker_gmlist[["Pericyte"]],
    cellmarker_gmlist[["Mural cell"]],
    cellmarker_gmlist[["Smooth muscle cell"]],
    cellmarker_gmlist[["Vascular smooth muscle cell(VSMC)"]],
    cellmarker_gmlist[["Ependymal cell"]],
    cellmarker_gmlist[["Leptomeningeal cell"]]
  )
)

saveRDS(cellmarker_gmlist, "processed_data/cellmarker_gmlist.rds")

###############################################################################
#                           OUTPUT DIRECTORIES                                #
###############################################################################

# path1_1 <- "./figures/6.CellTypeAnnot/6.1.ClusteringGV/6.1.1.SCINA/"

path2_1 <- "./figures/6.CellTypeAnnot/6.2.ClusteringGV_dw/6.2.1.SCINA/"

# path3_1 <- "./figures/6.CellTypeAnnot/6.3.ClusteringPoisson/6.3.1.SCINA/"

###############################################################################
#                SCINA ANALYSIS: CELLULAR LEVEL, UNKNOWNS: YES                #
###############################################################################

# # ---------------------------- modelGeneVar -----------------------------------
# ### Data preprocessing
# # Filter each cell type to keep only the common GMs
# gm.cellmarker.present <- lapply(gm.cellmarker, function(gm) {
#   unique(gm[gm %in% rowData(scegv_clusters)$gene_symbol])
# })
# 
# # Limit the number of GMs per cell type to avoid later error annotating with SCINA
# subset_gm <- sapply(gm.cellmarker.present, function(x) head(x,130))
# 
# sce.filtgv <- scegv_clusters[rowData(scegv_clusters)$gene_symbol %in% unlist(subset_gm),]
# sce.matrix.filtgv <- as.matrix(assay(sce.filtgv, "logcounts"))
# rownames(sce.matrix.filtgv) <- rowData(sce.filtgv)$gene_symbol
# 
# ### Annotation
# resultsSCINA = SCINA(sce.matrix.filtgv, subset_gm, max_iter = 100,
#                      convergence_n = 10, convergence_rate = 0.999,
#                      sensitivity_cutoff = 0.9, rm_overlap=TRUE,
#                      allow_unknown=TRUE, log_file='SCINA.log')
# 
# colData(scegv_clusters)$SCINA_label <- resultsSCINA$cell_labels
# 
# ### Plots
# # Obtain the cell types
# cell_types <- unique(scegv_clusters$SCINA_label)
# 
# # Assign the default viridis colors
# personal_colors <- viridis::viridis(length(cell_types))
# names(personal_colors) <- cell_types
# 
# # Modify the color of epitelial cells
# personal_colors[["Celulas_epiteliales"]] <- "#F77FBE" 
# 
# ## TSNE
# # Cluster40
# set.seed(100)
# plot <- plotReducedDim(scegv_clusters, dimred="TSNE120", colour_by="SCINA_label", text_by="cluster40") +
#   ggtitle("t-SNE perplexity = 120 y k = 40") +
#   scale_color_manual(name="Tipo celular", values=personal_colors) +
#   labs(x = "TSNE1", y = "TSNE2")
# ggsave(filename=paste0(path1_1, "TSNE120_k40.jpeg", sep=""),
#        plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
# ggsave(filename=paste0(path1_1, "TSNE120_k40.svg", sep=""),
#        plot=plot, width=10, height=6, units="in", device="svg")
# 
# ## UMAP
# # Cluster40
# set.seed(100)
# plot <- plotReducedDim(scegv_clusters, dimred="UMAP300", colour_by="SCINA_label", text_by="cluster40") +
#   ggtitle("UMAP n_neighbors = 300 y k = 40") +
#   scale_color_manual(name="Tipo celular", values=personal_colors) +
#   labs(x = "UMAP1", y = "UMAP2")
# ggsave(filename=paste0(path1_1, "UMAP300_k40.jpeg", sep=""),
#        plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
# ggsave(filename=paste0(path1_1, "UMAP300_k40.svg", sep=""),
#        plot=plot, width=10, height=6, units="in", device="svg")
# 
# ## Heatmap
# # Cluster40
# chosen_cluster <- "cluster40"
# tab <- table(Assigned = scegv_clusters$SCINA_label,
#              Cluster = scegv_clusters[[chosen_cluster]])
# plot <- pheatmap::pheatmap(log2(tab + 10),
#                            cluster_rows = FALSE,
#                            cluster_cols = FALSE,
#                            color = colorRampPalette(c("white", "blue"))(101))
# ggsave(filename=paste0(path1_1, "Heatmap_k40.jpeg", sep=""),
#        plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
# ggsave(filename=paste0(path1_1, "Heatmap_k40.svg", sep=""),
#        plot=plot, width=10, height=6, units="in", device="svg")

# -------------------- modelGeneVar density.weights = FALSE -------------------

### Data preprocessing
# Filter each cell type to keep only the common GMs
gm.cellmarker.present <- lapply(gm.cellmarker, function(gm) {
  unique(gm[gm %in% rowData(scegv_dw_clusters)$gene_symbol])
})

# Limit the number of GMs per cell type to avoid later error annotating with SCINA
subset_gm <- sapply(gm.cellmarker.present, function(x) head(x,130))

sce.filtgv_dw <- scegv_dw_clusters[rowData(scegv_dw_clusters)$gene_symbol %in% unlist(subset_gm),]
sce.matrix.filtgv_dw <- as.matrix(assay(sce.filtgv_dw, "logcounts"))
rownames(sce.matrix.filtgv_dw) <- rowData(sce.filtgv_dw)$gene_symbol

### Annotation
resultsSCINA = SCINA(sce.matrix.filtgv_dw, subset_gm, max_iter = 100,
                     convergence_n = 10, convergence_rate = 0.999,
                     sensitivity_cutoff = 0.9, rm_overlap=TRUE,
                     allow_unknown=TRUE, log_file='SCINA.log')

colData(scegv_dw_clusters)$SCINA_label <- resultsSCINA$cell_labels

# Obtain the cell types
cell_types <- unique(scegv_dw_clusters$SCINA_label)

# Assign the default viridis colors
personal_colors <- viridis::viridis(length(cell_types))
names(personal_colors) <- cell_types

# Modify the color of epitelial cells
personal_colors["Celulas_epiteliales"] <- "#F77FBE" 


### Plots 
## TSNE
# Cluster40
set.seed(100)
plot <- plotReducedDim(scegv_dw_clusters, dimred="TSNE120", colour_by="SCINA_label", text_by="cluster40") +
  ggtitle("t-SNE perplexity = 120 y k = 40") +
  scale_color_manual(name="Tipo celular", values=personal_colors) +
  labs(x = "TSNE1", y = "TSNE2")
ggsave(filename=paste0(path2_1, "TSNE120_k40.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2_1, "TSNE120_k40.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

## UMAP
# Cluster40
set.seed(100)
plot <- plotReducedDim(scegv_dw_clusters, dimred="UMAP300", colour_by="SCINA_label", text_by="cluster40") +
  ggtitle("UMAP n_neighbors = 300 y k = 40") +
  scale_color_manual(name="Tipo celular", values=personal_colors) +
  labs(x = "UMAP1", y = "UMAP2")
ggsave(filename=paste0(path2_1, "UMAP300_k40.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2_1, "UMAP300_k40.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

## Heatmap
# Cluster40
chosen_cluster <- "cluster40"
tab <- table(Assigned = scegv_dw_clusters$SCINA_label,
             Cluster = scegv_dw_clusters[[chosen_cluster]])
plot <- pheatmap::pheatmap(log2(tab + 10),
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           color = colorRampPalette(c("white", "blue"))(101))
ggsave(filename=paste0(path2_1, "Heatmap_k40.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2_1, "Heatmap_k40.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

# ------------------------ modelGeneVarByPoisson ----------------------------
# 
# ### Data preprocessing
# # Filter each cell type to keep only the common GMs
# gm.cellmarker.present <- lapply(gm.cellmarker, function(gm) {
#   unique(gm[gm %in% rowData(sceP_clusters)$gene_symbol])
# })
# 
# # Limit the number of GMs per cell type to avoid later error annotating with SCINA
# subset_gm <- sapply(gm.cellmarker.present, function(x) head(x,130))
# 
# sce.filtP <- sceP_clusters[rowData(sceP_clusters)$gene_symbol %in% unlist(subset_gm),]
# sce.matrix.filtP <- as.matrix(assay(sce.filtP, "logcounts"))
# rownames(sce.matrix.filtP) <- rowData(sce.filtP)$gene_symbol
# 
# ### Annotation
# resultsSCINA = SCINA(sce.matrix.filtP, subset_gm, max_iter = 100,
#                      convergence_n = 10, convergence_rate = 0.999,
#                      sensitivity_cutoff = 0.9, rm_overlap=TRUE,
#                      allow_unknown=TRUE, log_file='SCINA.log')
# 
# colData(sceP_clusters)$SCINA_label <- resultsSCINA$cell_labels
# 
# ### Plots
# ## TSNE
# # Cluster40
# set.seed(100)
# plot <- plotReducedDim(sceP_clusters, dimred="TSNE120", colour_by="SCINA_label", text_by="cluster40") +
#   ggtitle("t-SNE perplexity = 120 y k = 40") +
#   scale_color_manual(name="Tipo celular", values=personal_colors) +
#   labs(x = "TSNE1", y = "TSNE2")
# ggsave(filename=paste0(path3_1, "TSNE120_k40.jpeg", sep=""),
#        plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
# ggsave(filename=paste0(path3_1, "TSNE120_k40.svg", sep=""),
#        plot=plot, width=10, height=6, units="in", device="svg")
# 
# ## UMAP
# # Cluster40
# set.seed(100)
# plot <- plotReducedDim(sceP_clusters, dimred="UMAP300", colour_by="SCINA_label", text_by="cluster40") +
#   ggtitle("UMAP n_neighbors = 300 y k = 40") +
#   scale_color_manual(name="Tipo celular", values=personal_colors) +
#   labs(x = "UMAP1", y = "UMAP2")
# ggsave(filename=paste0(path3_1, "UMAP300_k40.jpeg", sep=""),
#        plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
# ggsave(filename=paste0(path3_1, "UMAP300_k40.svg", sep=""),
#        plot=plot, width=10, height=6, units="in", device="svg")
# 
# ## Heatmap
# # Cluster40
# chosen_cluster <- "cluster40"
# tab <- table(Assigned = sceP_clusters$SCINA_label,
#              Cluster = sceP_clusters[[chosen_cluster]])
# plot <- pheatmap::pheatmap(log2(tab + 10),
#                            cluster_rows = FALSE,
#                            cluster_cols = FALSE,
#                            color = colorRampPalette(c("white", "blue"))(101))
# ggsave(filename=paste0(path3_1, "Heatmap_k40.jpeg", sep=""),
#        plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
# ggsave(filename=paste0(path3_1, "Heatmap_k40.svg", sep=""),
#        plot=plot, width=10, height=6, units="in", device="svg")

# -------- Save the data ------------
# saveRDS(scegv_clusters, "processed_data/scegv_clusters_annot1.rds")
saveRDS(scegv_dw_clusters, "processed_data/scegv_dw_clusters_annot1.rds")
# saveRDS(sceP_clusters, "processed_data/sceP_clusters_annot1.rds")
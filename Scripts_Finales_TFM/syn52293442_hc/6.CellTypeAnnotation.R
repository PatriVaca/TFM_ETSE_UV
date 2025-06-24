# Single cell RNA-seq data analysis
# 6.Cell Type Annotation
# Author: Patricia del Carmen Vaca Rubio
# Date: 05/25

setwd("/clinicfs/projects/i63/tfm_hipocampo/syn52293442_hc/")
system("mkdir -p figures/6.CellTypeAnnot/6.2.ClusteringGV_dw/6.2.1.SCINA")
system("mkdir -p figures/6.CellTypeAnnot/6.2.ClusteringGV_dw/6.2.2.SingleR")
system("mkdir -p figures/6.CellTypeAnnot/6.2.ClusteringGV_dw/6.2.3.CHETAH")

library("dplyr")
library("HGNChelper")
library("SingleCellExperiment")
library("scater")
library("SCINA")
library("SingleR")
library("CHETAH")

# There are 3 genetic markers (GMs) datasets to annotate with: The Human
# Protein Atlas dataset (genes_marcadores.rds), the experimentally annotated
# csv (cellmarker_experiment.csv) and the sctype database.

# We will perform the cell type annotation with 3 packages: SCINA, singleR and
# CHETAH

###############################################################################
#  GENETIC MARKERS (GM) LIST PREPARATION FOR CELL TYPE ANNOTATION WITH SCINA  #
###############################################################################

### 1. INTERSECTION BETWEEN CELLMARKER + SCTYPE GMs WITH THE SCE GENES ###
## 1.1 CELLMARKER (only brain GMs) ##
cellmarkerdb <- read.csv("../notacion_celular/genes_marcadores/cellmarker_experiment.csv")

# Keep only the healthy condition genetic markers (not cancer)
cellmarkerdb <- cellmarkerdb[cellmarkerdb$Cancer=="Normal cell",]

# Keep the tissues of interest
cellmarkerdb <- cellmarkerdb[cellmarkerdb$Tissue.Type %in% c("Brain", "Hippocampus"),]

# Get each brain cell types GMs
cellmarker_gmdf <- cellmarkerdb %>%  
  group_by(`Cell.name`) %>%  
  summarise(geneset = list(`Cell.marker`))

# Conversion of df into a list
cellmarker_gmlist <- setNames(cellmarker_gmdf$geneset, cellmarker_gmdf$`Cell.name`)

# Select the cell types of interest's GMs to annotate the clusters
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
    cellmarker_gmlist[["Microglia-like cell"]]
  ),
  "OPC" = c(
    cellmarker_gmlist[["Oligodendrocyte precursor cell"]],
    cellmarker_gmlist[["Oligodendrocyte progenitor cell"]],
    cellmarker_gmlist[["Oligodendrocyte-like cell"]]
  ),
  "Neuronas excitadoras" = c(
    cellmarker_gmlist[["Excitatory neuron"]],
    cellmarker_gmlist[["Glutamatergic neuron"]],
    cellmarker_gmlist[["Deep layer neuron"]],
    cellmarker_gmlist[["Upper layer cortical neuron"]]
  ),
  "Neuronas inhibidoras" = c(
    cellmarker_gmlist[["Inhibitory neuron"]],
    cellmarker_gmlist[["GABAergic neuron"]],
    cellmarker_gmlist[["Interneuron"]],
    cellmarker_gmlist[["Cholinergic neuron"]],
    cellmarker_gmlist[["Dopaminergic neuron"]],
    cellmarker_gmlist[["Serotonergic neuron"]]
  )
)

## 1.2 SCTYPE (https://github.com/IanevskiAleksandr/sc-type) ##
# Load libraries and functions
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# Get cell-type-specific gene sets from their in-built database (DB)
sctypedb <- gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Brain")[[1]]

# Select the cell types of interest's GMs to annotate the clusters
gm.sctype <- list(
  "Astrocitos" = sctypedb[["Astrocytes"]],
  "Oligodendrocitos" = c(
    sctypedb[["Oligodendrocytes"]],
    sctypedb[["Myelinating Schwann cells"]]),
  "Microglia" = sctypedb[["Microglial cells"]],
  "OPC" = sctypedb[["Oligodendrocyte precursor cells"]],
  "Neuronas excitadoras" = sctypedb[["Glutamatergic neurons"]],
  "Neuronas inhibidoras" = c(
    sctypedb[["GABAergic neurons"]],
    sctypedb[["Cholinergic neurons"]],
    sctypedb[["Dopaminergic neurons"]],
    sctypedb[["Serotonergic neurons"]]
  )
)

## 1.3 INTERSECTION ##
sce.norm <- readRDS("processed_data/sce_norm.rds")

intersection_list <- lapply(names(gm.sctype), function(cell_type) {
  intersect(gm.sctype[[cell_type]], gm.cellmarker[[cell_type]])
})

names(intersection_list) <- names(gm.sctype)
# names(intersection_list) = c("astrocitos", "oligodendrocitos", "microglia", "opc", "excitadoras", "inhibidoras")
# gm.combinacion = c(gm.cellmarker, gm.sctype)

# names(intersection_list)

excitadoras <- unlist(intersection_list[["Neuronas excitadoras"]])
excitadoras <- intersect(excitadoras, rownames(sce.norm)) # Only keep the GMs present in our dataset (sce)

inhibidoras <- unlist(intersection_list[["Neuronas inhibidoras"]])
inhibidoras <- intersect(inhibidoras, rownames(sce.norm))

astrocitos <- unlist(intersection_list[["Astrocitos"]])
astrocitos <- intersect(astrocitos, rownames(sce.norm))

oligodendrocitos <- unlist(intersection_list[["Oligodendrocitos"]])
oligodendrocitos <- intersect(oligodendrocitos, rownames(sce.norm))

opc <- unlist(intersection_list[["OPC"]])
opc <- intersect(opc, rownames(sce.norm))

microglia <- unlist(intersection_list["Microglia"])
microglia <- intersect(microglia, rownames(sce.norm))

gm.manual <- list(excitadoras, inhibidoras, astrocitos, oligodendrocitos, opc, microglia)
names(gm.manual) <- c("Neuronas excitadoras", "Neuronas inhibidoras", "Astrocitos",
                     "Oligodendrocitos","OPC","Microglia")
saveRDS(gm.manual, "processed_data/gm_manual.rds")

### 2. INTERSECTION BETWEEN HUMAN PROTEIN ATLAS (HPA) GMs WITH THE SCE GENES ###
gm.manual2 <- readRDS("../notacion_celular/genes_marcadores/genes_marcadores.rds")

excitadoras <- unlist(gm.manual2[["Ex"]])
excitadoras <- intersect(excitadoras, rownames(sce.norm))

inhibidoras <- unlist(gm.manual2[["Inh"]])
inhibidoras <- intersect(inhibidoras, rownames(sce.norm))

astrocitos <- unlist(gm.manual2[["Astrocitos"]])
astrocitos <- intersect(astrocitos, rownames(sce.norm))

oligodendrocitos <- unlist(gm.manual2["Oligodendrocito"])
oligodendrocitos <- intersect(oligodendrocitos, rownames(sce.norm))

opc <- unlist(gm.manual2["OPC"])
opc <- intersect(opc, rownames(sce.norm))

microglia <- unlist(gm.manual2["Microglia"])
microglia <- intersect(microglia, rownames(sce.norm))

endotelio <- unlist(gm.manual2["Endothelial"])
endotelio <- intersect(endotelio, rownames(sce.norm))

fibroblasto <- unlist(gm.manual2["Fibroblastos"])
fibroblasto <- intersect(fibroblasto, rownames(sce.norm))

# HPA combination
gm.manual.combinado <- list(excitadoras, inhibidoras,
                           astrocitos, oligodendrocitos,opc,microglia,
                           endotelio,fibroblasto)
names(gm.manual.combinado) <- c("Neuronas excitadoras", "Neuronas inhibidoras",
                               "Astrocitos", "Oligodendrocitos", "OPC", "Microglia",
                               "Endotelio", "Fibroblasto")
saveRDS(gm.manual.combinado, "processed_data/gm_manual_combinado.rds")
# gm.manual.combinado --> contains common GMs between cellmarker_experiment.csv,
# sctype, rownames(sce.norm) and HPA from genes_marcadores.rds


# Combine the mannual annotation from the 2 dbs (cellmarker + sctype + genes present
# in sce.norm) with the complete combined annotation (gm.manual.combi_completo: HPA + 
# genes present in sce.norm)

# gm.manual <- readRDS("processed_data/gm_manual.rds")
# gm.manual.combinado <- readRDS("processed_data/gm_manual_combinado.rds")

gm.manual.3 <- lapply(names(gm.manual), function(cell_type) {
  c(gm.manual[[cell_type]], gm.manual.combinado[[cell_type]])
})
names(gm.manual.3) = c("Excitadoras", "Inhibidoras", "Astrocitos", "Oligodendrocitos","OPC","Microglia")

saveRDS(gm.manual.3, "processed_data/gm_manual_3.rds")

###############################################################################
#                            DATA PREPROCESSING                               #
###############################################################################

scegv_dw_clusters <- readRDS("processed_data/scegv_dw_clusters.rds")
sce.filtgv_dw <- scegv_dw_clusters[rownames(scegv_dw_clusters) %in% unlist(gm.manual.3),]
sce.matrix.filtgv_dw <- as.matrix(assay(sce.filtgv_dw, "logcounts"))

###############################################################################
#                           OUTPUT DIRECTORIES                                #
###############################################################################

path2_1 <- "./figures/6.CellTypeAnnot/6.2.ClusteringGV_dw/6.2.1.SCINA/"
path2_2 <- "./figures/6.CellTypeAnnot/6.2.ClusteringGV_dw/6.2.2.SingleR/"
path2_3 <- "./figures/6.CellTypeAnnot/6.2.ClusteringGV_dw/6.2.3.CHETAH/"

###############################################################################
#                      SCINA ANALYSIS -> UNKNOWNS: YES                        #
###############################################################################
# Cell type annotation with SCINA using gm.manual.3 (cellmarkercsv + sctype
# + hpa + rownames(scenorm))

# gm.manual.3 <- readRDS("processed_data/gm_manual_3.rds")

# ---------------------- modelGeneVar density.weights = FALSE -----------------

resultsSCINA <- SCINA(sce.matrix.filtgv_dw, gm.manual.3, max_iter = 100,
                     convergence_n = 10, convergence_rate = 0.999,
                     sensitivity_cutoff = 0.9, rm_overlap=TRUE,
                     allow_unknown=TRUE, log_file='SCINA.log')

# View(resultsSCINA$cell_labels)
# View(resultsSCINA$probabilities)

colData(scegv_dw_clusters)$SCINA_label <- resultsSCINA$cell_labels

## TSNE
# Cluster40
set.seed(100)
plot <- plotReducedDim(scegv_dw_clusters, dimred="TSNE120", colour_by="SCINA_label", text_by="cluster40") +
  ggtitle("t-SNE perplexity = 120 y k = 40") +
  scale_color_viridis_d(name="Tipo celular") +
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
  scale_color_viridis_d(name="Tipo celular") +
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

# -------- Save the data ------------
saveRDS(scegv_dw_clusters, "processed_data/scegv_dw_clusters_annot1.rds")

###############################################################################
#                 ANALISIS CON SINGLER: NIVEL CELULAR, UNKNOWN: NO            #
###############################################################################

sce.downsample <- readRDS("../notacion_celular/human_atlas_dataset/sce_allan_red_downsample.rds")

sce.downsample$tipos_celulares[sce.downsample$tipos_celulares == "Exc"] <- "Excitadoras"
sce.downsample$tipos_celulares[sce.downsample$tipos_celulares == "Inh"] <- "Inhibidoras"
sce.downsample$tipos_celulares[sce.downsample$tipos_celulares == "Astro"] <- "Astrocitos"
sce.downsample$tipos_celulares[sce.downsample$tipos_celulares == "Oligo"] <- "Oligodendrocitos"
sce.downsample$tipos_celulares[sce.downsample$tipos_celulares == "Micro"] <- "Microglia"

# -------------------- modelGeneVar density.weights = FALSE -------------------

resultsSingleR <- SingleR(test=sce.matrix.filtgv_dw,
                          ref=sce.downsample,
                          labels=sce.downsample$tipos_celulares,
                          assay.type.test="logcounts",
                          assay.type.ref="logcounts")

colData(scegv_dw_clusters)$singleR_label <- resultsSingleR$labels

## TSNE
# Cluster40
set.seed(100)
plot <- plotReducedDim(scegv_dw_clusters, dimred="TSNE120", colour_by="singleR_label", text_by="cluster40") +
  ggtitle("t-SNE perplexity = 120 y k = 40") +
  scale_color_viridis_d(name="Tipo celular") +
  labs(x = "TSNE1", y = "TSNE2")
ggsave(filename=paste0(path2_2, "TSNE120_k40.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2_2, "TSNE120_k40.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")


## UMAP
# Cluster40
set.seed(100)
plot <- plotReducedDim(scegv_dw_clusters, dimred="UMAP300", colour_by="singleR_label", text_by="cluster40") +
  ggtitle("UMAP n_neighbors = 300 y k = 40") +
  scale_color_viridis_d(name="Tipo celular") +
  labs(x = "UMAP1", y = "UMAP2")
ggsave(filename=paste0(path2_2, "UMAP300_k40.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2_2, "UMAP300_k40.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

## ScoreHeatmap
plot <- plotScoreHeatmap(resultsSingleR)
ggsave(filename=paste0(path2_2, "ScoreHeatmap.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2_2, "ScoreHeatmap.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

## Heatmap
# Cluster40
chosen_cluster <- "cluster40"
tab <- table(Assigned = scegv_dw_clusters$singleR_label,
             Cluster = scegv_dw_clusters[[chosen_cluster]])
plot <- pheatmap::pheatmap(log2(tab + 10),
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           color = colorRampPalette(c("white", "blue"))(101))
ggsave(filename=paste0(path2_2, "Heatmap_k40.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2_2, "Heatmap_k40.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

# -------- Save the data ------------
saveRDS(scegv_dw_clusters, "processed_data/scegv_dw_clusters_annot2.rds")

###############################################################################
#          CHETAH ANALYSIS -> UNKNOWNS: NO, INTERMEDIUM NODES: YES            #
###############################################################################
# -------------------- modelGeneVar density.weights = FALSE -------------------

resultsCHETAHgv_dw <- CHETAHclassifier(input=scegv_dw_clusters,
                                       ref_cells=sce.downsample,
                                       ref_ct="tipos_celulares",
                                       input_c="logcounts",
                                       ref_c="logcounts")

colData(scegv_dw_clusters)$CHETAH_label <- resultsCHETAHgv_dw$celltype_CHETAH

## TSNE
# Cluster40
set.seed(100)
plot <- plotReducedDim(scegv_dw_clusters, dimred="TSNE120", colour_by="CHETAH_label", text_by="cluster40") +
  ggtitle("t-SNE perplexity = 120 y k = 40") +
  scale_color_viridis_d(name="Tipo celular") +
  labs(x = "TSNE1", y = "TSNE2")
ggsave(filename=paste0(path2_3, "TSNE120_k40.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2_3, "TSNE120_k40.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

## UMAP
# Cluster40
set.seed(100)
plot <- plotReducedDim(scegv_dw_clusters, dimred="UMAP300", colour_by="CHETAH_label", text_by="cluster40") +
  ggtitle("UMAP n_neighbors = 300 y k = 40") +
  scale_color_viridis_d(name="Tipo celular") +
  labs(x = "UMAP1", y = "UMAP2")
ggsave(filename=paste0(path2_3, "UMAP300_k40.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2_3, "UMAP300_k40.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

## Heatmap
# Cluster40
chosen_cluster <- "cluster40"
tab <- table(Assigned = scegv_dw_clusters$CHETAH_label,
             Cluster = scegv_dw_clusters[[chosen_cluster]])
plot <- pheatmap::pheatmap(log2(tab + 10),
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           color = colorRampPalette(c("white", "blue"))(101))
ggsave(filename=paste0(path2_3, "Heatmap_k40.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2_3, "Heatmap_k40.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

# -------- Save the data ------------
saveRDS(scegv_dw_clusters, "processed_data/scegv_dw_clusters_annot3.rds")

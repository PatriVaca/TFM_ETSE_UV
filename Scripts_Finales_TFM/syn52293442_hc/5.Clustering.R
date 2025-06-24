# Single cell RNA-seq data analysis
# 5.Clustering
# Author: Patricia del Carmen Vaca Rubio
# Date: 05/25

setwd("/clinicfs/projects/i63/tfm_hipocampo/syn52293442_hc/")
system("mkdir -p figures/5.Clustering/5.2.ClusteringGV_dw")

# Load packages
library(scran)
library(bluster)
library(scater)
library(dplyr)
library(ggplot2)

# ------------------------------ FUNCTIONS ------------------------------------

## Clustering function ##
do_clustering <- function(sce, knum, louvain_list) {
  louvain <- scran::clusterCells(sce,
                                 use.dimred = "PCA",
                                 BLUSPARAM=NNGraphParam(k=knum, cluster.fun="louvain"))
  
  # Save the clusters assignations in a list
  louvain_list[[paste0("k", knum)]] <- louvain
  
  # Save the clusters assignations also in the colData
  colData(sce)[[paste0("cluster", knum)]] <- factor(louvain)
  
  return(list(sce, louvain_list))
}

## Silhouette width plot function ##
SilhouetteWidthEvaluation <- function(sce, clusters_k, knum, tablesilh_list) {
  # Silhouette width calculation
  sil.approx <- approxSilhouette(reducedDim(sce, "PCA"), clusters_k)
  sil.data <- as.data.frame(sil.approx)
  sil.data$best.choice <- factor(ifelse(sil.data$width > 0, clusters_k, sil.data$other))
  sil.data$clusters <- factor(clusters_k)
  
  # Save the stats
  tablek <- table(Assigned=sil.data$clusters, Closest=sil.data$best.choice)
  tablesilh_list[[paste0("k", knum)]] <- tablek
  
  # Plot Silhouette plot
  silhplot <- ggplot(sil.data, aes(x=clusters, y=width, colour=best.choice)) +
    ggbeeswarm::geom_quasirandom(method="smiley") +
    ggtitle(paste("Silhouette plot - k =", knum)) +
    ylab("Width") +
    xlab("Cluster")
  
  return(list(tablesilh_list, silhplot))
}

## Cluster Purity plot function ##
ClusterPurityEvaluation <- function(sce, clusters_k, knum, tablespure_list) {
  ## Cluster purity
  pure <- neighborPurity(reducedDim(sce, "PCA"), clusters_k)
  pure.data <- as.data.frame(pure)
  pure.data$maximum <- factor(pure.data$maximum)
  pure.data$cluster <- factor(clusters_k)
  
  # Save the statss
  tablek <- table(Cluster=pure.data$cluster, pure.data$maximum)
  tablespure_list[[paste0("k", knum)]] <- tablek
  
  # Plot Cluster Purity
  purityplot <- ggplot(pure.data, aes(x=cluster, y=purity, colour=maximum)) +
    ggbeeswarm::geom_quasirandom(method="smiley") +
    ggtitle(paste("Purity plot - k =", knum)) +
    ylab("Purity") +
    xlab("Cluster")
  
  return (list(tablespure_list, purityplot))
}

## RMSD intra-cluster ##
RMSDcluster <- function(sce, clusters_k, knum, rmsd_list) {
  rmsd <- clusterRMSD(reducedDim(sce, "PCA"), clusters_k)
  rmsd.df <- data.frame(rmsd = rmsd, cluster = names(rmsd))
  rmsd.df$cluster <- factor(rmsd.df$cluster, levels = sort(unique(as.numeric(rmsd.df$cluster))))
  
  # Save the values
  rmsd_list[[paste0("k", knum)]] <- rmsd.df
  
  # Plot RMSD intra-cluster
  prmsd <- ggplot(rmsd.df, aes(x=cluster, y=rmsd, fill = cluster))+
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 5, linetype = "dashed", color = "grey") +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
    theme_minimal() +
    ggtitle(paste("RMSD intra-cluster - k =", knum)) +
    ylab("RMSD") +
    xlab("Cluster")
  
  return(list(rmsd_list, prmsd, rmsd))
}

## RMSD vs Size Factors Relation ##
RMSD_vs_SF <- function(sce, clusters_k, knum, sf_list, rmsd) {
  by.clust <- split(sizeFactors(sce), clusters_k)
  sf.by.clust <- vapply(by.clust, mean, 0)
  point.df <- data.frame(rmsd = rmsd, sf = sf.by.clust, label = names(rmsd))
  
  # Save the values
  sf_list[[paste0("k", knum)]] <- point.df
  
  # Plot the RMSD vs Size Factors Relation
  pointplot <- ggplot(point.df, aes(x = rmsd, y = sf)) +
    geom_point() +
    ggrepel::geom_text_repel(aes(label = label)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
    geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
    geom_vline(xintercept = 5, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
    theme_minimal() +
    ggtitle(paste("RMSD vs Size Factors - k =", knum)) +
    xlab("RMSD") +
    ylab("Mean Size Factor")
  
  return(list(sf_list, pointplot))
}

## PCA ##
PCA_evaluation <- function(sce, knum, folder) {
  jpeg(file=paste0("./figures/5.Clustering/", folder, "/PCA_by_clusters", knum, ".jpeg", sep=""),
       width=10, height=6, units="in", res=300)
  print(
    plotReducedDim(sce,
                   dimred="PCA",
                   colour_by=paste0("cluster", knum, sep=""),
                   text_by=paste0("cluster", knum, sep="")) +
      ggtitle(paste0("PCA 8 componentes y k = ", knum, sep="")) +
      xlab("PCA1")+
      ylab("PCA2")+
      scale_color_viridis_d(name = "Grupos")
  )
  dev.off()
}

## tSNE
tSNE_evaluation <- function(sce, knum, folder) {
  jpeg(file=paste0("./figures/5.Clustering/", folder, "/TSNE120_by_clusters", knum, ".jpeg", sep=""),
       width=10, height=6, units="in", res=300)
  print(
    plotReducedDim(sce,
                   dimred="TSNE120",
                   colour_by=paste0("cluster", knum, sep=""),
                   text_by=paste0("cluster", knum, sep="")) +
      ggtitle(paste0("t-SNE perplexity = 120 y k = ", knum, sep="")) +
      xlab("tSNE1")+
      ylab("tSNE2")+
      scale_color_viridis_d(name = "Grupos")
  )
  dev.off()
}

## UMAP
UMAP_evaluation <- function(sce, knum, folder) {
  jpeg(file=paste0("./figures/5.Clustering/", folder, "/UMAP300_by_clusters", knum, ".jpeg", sep=""),
       width=10, height=6, units="in", res=300)
  
  print(
    plotReducedDim(sce,
                   dimred="UMAP300",
                   colour_by=paste0("cluster", knum, sep=""),
                   text_by=paste0("cluster", knum, sep="")) +
      ggtitle(paste0("UMAP n_neighbors = 300 y k = ", knum, sep="")) +
      xlab("UMAP1")+
      ylab("UMAP2")+
      scale_color_viridis_d(name = "Grupos")
  )
  dev.off()
}


# ---------------------------------- MAIN -------------------------------------

k_list <- list(5, 15, 20, 30, 40)
# k=25 triggers an error in this study

########################################
## ModelGeneVar density.weights=FALSE ##
########################################

scegv_dw <- readRDS("processed_data/scegv_dw.rds")

louvain_list_gv_dw <- list()
tablesilh_list_gv_dw <- list()
tablespure_list_gv_dw <- list()
rmsd_list_gv_dw <- list()
sf_list_gv_dw <- list()

for (knum in k_list) {
  ## Clustering ##
  set.seed(100)
  clustering_results <- do_clustering(scegv_dw, knum, louvain_list_gv_dw)
  scegv_dw <- clustering_results[[1]]
  louvain_list_gv_dw <- clustering_results[[2]]
  
  ## Quantifying clustering behavior ##
  # Current cluster (of scegv_dw in this case)
  clusters_k <- colData(scegv_dw)[[paste0("cluster", knum)]]
  
  # Silhouette width
  silh <- SilhouetteWidthEvaluation(scegv_dw, clusters_k, knum, tablesilh_list_gv_dw)
  tablesilh_list_gv_dw <- silh[[1]]
  silhplot <- silh[[2]]
  # Save the plot
  ggsave(paste0("./figures/5.Clustering/5.2.ClusteringGV_dw/Silhouette_k", knum, ".jpeg", sep=""),
         plot=silhplot, width=10, height=6, units="in", device="jpeg", dpi=300)
  ggsave(paste0("./figures/5.Clustering/5.2.ClusteringGV_dw/Silhouette_k", knum, ".svg", sep=""),
         plot=silhplot, width=10, height=6, units="in", device="svg")
  
  # Purity
  pure <- ClusterPurityEvaluation(scegv_dw, clusters_k, knum, tablespure_list_gv_dw)
  tablespure_list_gv_dw <- pure[[1]]
  purityplot <- pure[[2]]
  # Save the plot
  ggsave(paste0("./figures/5.Clustering/5.2.ClusteringGV_dw/Purity_k", knum, ".jpeg", sep=""),
         plot=purityplot, width=10, height=6, units="in", device="jpeg", dpi=300)
  ggsave(paste0("./figures/5.Clustering/5.2.ClusteringGV_dw/Purity_k", knum, ".svg", sep=""),
         plot=purityplot, width=10, height=6, units="in", device="svg")
  
  # RMSD intra-cluster
  rmsd_results <- RMSDcluster(scegv_dw, clusters_k, knum, rmsd_list_gv_dw)
  rmsd_list_gv_dw <- rmsd_results[[1]]
  rmsd_plot <- rmsd_results[[2]]
  rmsd <- rmsd_results[[3]]
  # Save the plot
  ggsave(paste0("./figures/5.Clustering/5.2.ClusteringGV_dw/RMSD_k", knum, ".jpeg", sep=""),
         plot=rmsd_plot, width=10, height=6, units="in", device="jpeg", dpi=300)
  ggsave(paste0("./figures/5.Clustering/5.2.ClusteringGV_dw/RMSD_k", knum, ".svg", sep=""),
         plot=rmsd_plot, width=10, height=6, units="in", device="svg")
  
  # RMSD vs Size Factors Relation
  sfresults <- RMSD_vs_SF(scegv_dw, clusters_k, knum, sf_list_gv_dw, rmsd)
  sf_list_gv_dw <- sfresults[[1]]
  pointplot <- sfresults[[2]]
  # Save the plot
  ggsave(paste0("./figures/5.Clustering/5.2.ClusteringGV_dw/RMSD_vs_SF_k", knum, ".jpeg", sep=""),
         plot=pointplot, width=10, height=6, units="in", device="jpeg", dpi=300)
  ggsave(paste0("./figures/5.Clustering/5.2.ClusteringGV_dw/RMSD_vs_SF_k", knum, ".svg", sep=""),
         plot=pointplot, width=10, height=6, units="in", device="svg")
  
  folder <- "5.2.ClusteringGV_dw"
  # PCA
  PCA_evaluation(scegv_dw, knum, folder)
  
  # tSNE
  tSNE_evaluation(scegv_dw, knum, folder)
  
  # UMAP
  UMAP_evaluation(scegv_dw, knum, folder)
}

# ---- Save the data ------
saveRDS(scegv_dw, "processed_data/scegv_dw_clusters.rds")


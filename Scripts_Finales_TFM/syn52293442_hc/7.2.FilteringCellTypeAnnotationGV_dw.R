# Single cell RNA-seq data analysis
# 7.FilteringCellTypeAnnotationGV_DW
# Author: Patricia del Carmen Vaca Rubio
# Date: 05/25

# pacman::p_load(SingleCellExperiment, scater, scran, dplyr, BiocParallel, 
#                scDblFinder, AnnotationHub, tidyverse, patchwork, ggvenn, broom,
#                kableExtra,HDF5Array, viridis, DropletUtils,
#                devtools, celda, ggvenn, bluster, cluster, ggpubr, SCINA, SingleR,
#                CHETAH, pheatmap, ggplot2)

pacman::p_load(
  # single nucleus/single cells experiments 
  SingleCellExperiment,
  scater,
  # Data manipulation
  dplyr,
  tidyverse,
  # Visualization
  ggplot2,
  viridis,
  pheatmap,
  # Cell type annotation
  SingleR,
  CHETAH,
  SCINA,
  # Export
  openxlsx
)

setwd("/clinicfs/projects/i63/tfm_hipocampo/syn52293442_hc/")
system("mkdir -p figures/7.FilteredCellTypeAnnotation/7.2.ClusteringGV_dw/prefiltering")
system("mkdir -p figures/7.FilteredCellTypeAnnotation/7.2.ClusteringGV_dw/filtered")

scegv_dw_clusters <- readRDS("processed_data/scegv_dw_clusters_annot3.rds")

# ---------------------------- Prefiltering ----------------------------------

# CONSENSUAL RESULTS GENERATION #
path1 <- "figures/7.FilteredCellTypeAnnotation/7.2.ClusteringGV_dw/prefiltering/"

# Convert intermedium nodes into unknowns to obtain the 3 annotation methods consensum
scegv_dw_clusters$CHETAH_label <- str_replace(scegv_dw_clusters$CHETAH_label, "Node1", "unknown")
scegv_dw_clusters$CHETAH_label <- str_replace(scegv_dw_clusters$CHETAH_label, "Node2", "unknown")
scegv_dw_clusters$CHETAH_label <- str_replace(scegv_dw_clusters$CHETAH_label, "Node3", "unknown")
scegv_dw_clusters$CHETAH_label <- str_replace(scegv_dw_clusters$CHETAH_label, "Node4", "unknown")
scegv_dw_clusters$CHETAH_label <- str_replace(scegv_dw_clusters$CHETAH_label, "Unassigned", "unknown")
scegv_dw_clusters$SCINA_label <- str_replace(scegv_dw_clusters$SCINA_label, "endotelio", "unknown")
scegv_dw_clusters$SCINA_label <- str_replace(scegv_dw_clusters$SCINA_label, "fibroblasto", "unknown")

# Combine the 3 annotation methods results
all_res <-data.frame(
  scegv_dw_clusters$SCINA_label,
  scegv_dw_clusters$CHETAH_label,
  scegv_dw_clusters$singleR_label
)

# Obtain the consensus and generate an output tsv file
all_res$cell_type <- apply(all_res, 1, function(x) names(which.max(table(as.character(x)))))
all_res$n_suport <- apply(all_res[- dim(all_res)[2]], 1, function(x) max(table(as.character(x))))
all_res$cell_type[which(all_res$n_suport==1)] <-"unknown"
apply(all_res, 2, table)
write.table(all_res,file=paste0("processed_data/matrix_annot_gv_dw.tsv"), sep="\t", row.names = FALSE)

# Annotated SCE creation
sce.annot  <- scegv_dw_clusters
scegv_dw_clusters$tipo.celular  <- all_res$cell_type
sce.annot$tipo.celular  <- all_res$cell_type
saveRDS(scegv_dw_clusters, file  <- "processed_data/scegv_dw_clusters_filtered_annot.rds")

# FINAL RESULTS BEFORE FILTERING #
# Select a cluster
cluster <- "cluster40"

tab <- table(label = sce.annot$tipo.celular, cluster = colData(sce.annot)$cluster40) 
y  <- as.data.frame(tab)

# CELL TYPES PER CLUSTER
f  <- ggplot(y, aes(x=cluster, y=Freq, fill=label))
plot <- f + 
  geom_bar(stat = "identity") +
  scale_fill_viridis_d(name="Tipo celular") +
  ggtitle("Num cel por tipo celular") +
  xlab("Cluster") + 
  ylab("Núm células")
ggsave(filename=paste0(path1, "barplot_preanotacion.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path1, "barplot_preanotacion.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

# NUMBER OF CELLS PER CELL TYPE
f  <- ggplot(y, aes(x=label, y=Freq, fill=label)) 
plot <- f + geom_bar(stat="identity") +
  scale_fill_viridis_d(name="Tipo celular") +
  ggtitle("Num cel por tipo celular") +
  xlab("Tipo celular") + 
  ylab("Núm celulas")
ggsave(filename=paste0(path1, "barplot_preanotacion_abundancia.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path1, "barplot_preanotacion_abundancia.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

# HEATMAP
plot <- pheatmap::pheatmap(log2(tab+10), cluster_rows = F, cluster_cols = F,
                           color=colorRampPalette(c("white", "blue"))(101))
ggsave(filename=paste0(path1, "preanotacion_heatmap.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path1, "preanotacion_heatmap.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

# DELETED CELLS (UNKNOWNS) BARPLOT
delete_unknown <- grep("unknown", sce.annot$tipo.celular)  # Select the cells tagged as unknown
cluster_unknown <-  sce.annot@colData[delete_unknown,cluster] # Define what clusters contain those cells

deleted  <- table(cluster_unknown)  # Cells of the cluster that will be deleted
total  <- table(sce.annot$cluster40) # Cells assigned to each cluster
kept  <- total - deleted # Kept cells after discarding the unknowns

df.deleted  <- data.frame(deleted = deleted)
df.deleted$tipo  <- "eliminado"
colnames(df.deleted)  <- c("Cluster","Frecuencia","Tipo")

df.kept  <- data.frame(kept = kept)
df.kept$tipo  <- "mantenido"
colnames(df.kept)  <- c("Cluster","Frecuencia","Tipo")

combined.df  <- as.data.frame(rbind(df.deleted, df.kept))

f  <- ggplot(combined.df, aes(x=Cluster, y=Frecuencia, fill=Tipo)) 
plot <- f + geom_bar(stat="identity") +
  ggtitle("Celulas consenso") +
  xlab("Cluster") + 
  ylab("Núm celulas")
ggsave(filename=paste0(path1, "barplot_preanotacion_eliminados.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path1, "barplot_preanotacion_eliminados.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")


# FINAL PCA
a  <- plotReducedDim((sce.annot), dimred="PCA" , colour_by="tipo.celular", text_by= cluster)+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Anotación Final") 
a  <- a +  scale_color_viridis_d(name="Tipo celular")
ggsave(filename=paste0(path1, "pca_preanotacion.jpeg", sep=""),
       plot=a, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path1, "pca_preanotacion.svg", sep=""),
       plot=a, width=10, height=6, units="in", device="svg")


# FINAL TSNE
a  <- plotReducedDim((sce.annot), dimred="TSNE120", colour_by="tipo.celular", text_by= cluster)+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Anotación Final") +
  labs(x = "TSNE1", y = "TSNE2")
a <- a +  scale_color_viridis_d(name="Tipo celular")
ggsave(filename=paste0(path1, "tsne_preanotacion.jpeg", sep=""),
       plot=a, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path1, "tsne_preanotacion.svg", sep=""),
       plot=a, width=10, height=6, units="in", device="svg")

# FINAL TSNE 
a  <- plotReducedDim((sce.annot), dimred="TSNE120", colour_by="Group_sex", text_by= cluster)+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Anotación Final") +
  labs(x = "TSNE1", y = "TSNE2")
a <- a + scale_color_viridis_d(name="Condición y sexo")
ggsave(filename=paste0(path1, "tsne_preanotacion2.jpeg", sep=""),
       plot=a, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path1, "tsne_preanotacion2.svg", sep=""),
       plot=a, width=10, height=6, units="in", device="svg")

# FINAL UMAP
a <- plotReducedDim((sce.annot), dimred="UMAP300", colour_by="tipo.celular", text_by= cluster)+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Anotación Final")  +
  labs(x = "UMAP1", y = "UMAP2")
a <- a +  scale_color_viridis_d(name="Tipo celular")
ggsave(filename=paste0(path1, "umap_preanotacion.jpeg", sep=""),
       plot=a, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path1, "umap_preanotacion.svg", sep=""),
       plot=a, width=10, height=6, units="in", device="svg")

# Check if the "unknown" cell types are disease exclusive # 
# Extract UMAP coordinates and variables
umap_data <- as.data.frame(reducedDim(sce.annot, "UMAP300"))
umap_data$Group_sex <- colData(sce.annot)$Group_sex
umap_data$cell_type <- colData(sce.annot)$tipo.celular

a <- ggplot(umap_data, aes(x = V1, y = V2, color = cell_type)) +
  geom_point() +
  facet_wrap(~ Group_sex,
             labeller = labeller(Group_sex = c(
               "AD_Female" = "AD_Mujer",
               "AD_Male"   = "AD_Hombre",
               "Control_Female" = "Control_Mujer",
               "Control_Male"   = "Control_Hombre"
             ))) +
  theme(legend.key.size = unit(0.25, 'cm'))
a <- a +  scale_color_viridis_d(name="Tipo celular")
ggsave(filename=paste0(path1, "disease_cell_unknown.jpeg", sep=""),
       plot=a, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path1, "disease_cell_unknown.svg", sep=""),
       plot=a, width=10, height=6, units="in", device="svg")

# ---------------------------- Filtering ----------------------------------
table(sce.annot$cluster40, sce.annot$tipo.celular)
sce.annot =  sce.annot[,!sce.annot@colData[,"tipo.celular"] == "unknown"]

table(sce.annot$cluster40, sce.annot$tipo.celular)
# Manual filtering to remove the outliers of each cluster
sce.annot <- sce.annot[,!(sce.annot$tipo.celular != "Astrocitos" & sce.annot$cluster40 == 1)]
sce.annot <- sce.annot[,!(sce.annot$tipo.celular != "Astrocitos" & sce.annot$cluster40 == 2)]
sce.annot <- sce.annot[,!sce.annot$cluster40 == 3] # Delete whole cluster because the cell type is not clear
sce.annot <- sce.annot[,!(sce.annot$tipo.celular != "OPC" & sce.annot$cluster40 == 4)]
sce.annot <- sce.annot[,!(sce.annot$tipo.celular != "Oligodendrocitos" & sce.annot$cluster40 == 5)]
sce.annot <- sce.annot[,!(sce.annot$tipo.celular != "Excitadoras" & sce.annot$cluster40 == 6)]
sce.annot <- sce.annot[,!(sce.annot$tipo.celular != "Excitadoras" & sce.annot$cluster40 == 7)]
sce.annot <- sce.annot[,!(sce.annot$tipo.celular != "Excitadoras" & sce.annot$cluster40 == 8)]
sce.annot <- sce.annot[,!(sce.annot$tipo.celular != "Inhibidoras" & sce.annot$cluster40 == 9)]
sce.annot <- sce.annot[,!(sce.annot$tipo.celular != "Excitadoras" & sce.annot$cluster40 == 10)]
sce.annot <- sce.annot[,!(sce.annot$tipo.celular != "Excitadoras" & sce.annot$cluster40 == 11)]
sce.annot <- sce.annot[,!(sce.annot$tipo.celular != "Oligodendrocitos" & sce.annot$cluster40 == 12)]
sce.annot <- sce.annot[,!(sce.annot$tipo.celular != "Microglia" & sce.annot$cluster40 == 13)]
sce.annot <- sce.annot[,!(sce.annot$tipo.celular != "Microglia" & sce.annot$cluster40 == 14)]
sce.annot <- sce.annot[,!(sce.annot$tipo.celular != "OPC" & sce.annot$cluster40 == 15)]
sce.annot <- sce.annot[,!(sce.annot$tipo.celular != "Oligodendrocitos" & sce.annot$cluster40 == 16)]
sce.annot <- sce.annot[,!(sce.annot$tipo.celular != "Oligodendrocitos" & sce.annot$cluster40 == 17)]
sce.annot <- sce.annot[,!(sce.annot$tipo.celular != "Oligodendrocitos" & sce.annot$cluster40 == 18)]


table(sce.annot$cluster40, sce.annot$tipo.celular)

as.data.frame(table(sce.annot$cluster40, sce.annot$tipo.celular)) %>% 
  group_by(Var2)  %>% 
  summarise(cond_disp = sum(Freq))

path2 <- "figures/7.FilteredCellTypeAnnotation/7.2.ClusteringGV_dw/filtered/"

tab <- table(label=sce.annot$tipo.celular, cluster=sce.annot$cluster40) 
y = as.data.frame(tab)

# HEATMAP
plot <- pheatmap::pheatmap(log2(tab+10), cluster_rows = F, cluster_cols = F, color=colorRampPalette(c("white", "blue"))(101))
ggsave(filename=paste0(path2, "anotacion_heatmap.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "anotacion_heatmap.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")


# CELL TYPE PER CLUSTER BARPLOT
f <- ggplot(y, aes(x=cluster, y=Freq, fill=label)) 
plot <- f + geom_bar(stat="identity") +
  scale_fill_viridis_d(name="Tipo celular") +
  ggtitle("Num cel por tipo celular") +
  xlab("Cluster") +
  ylab("Núm celulas")
ggsave(filename=paste0(path2, "barplot_anotacion.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "barplot_anotacion.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

# DELETED CELLS BARPLOT
delete_unknown <- grep("unknown", sce.annot$tipo.celular)
cluster_unknown <- sce.annot@colData[delete_unknown,cluster]
deleted <- table(cluster_unknown)
total <- table(sce.annot@colData[,cluster])
kept <- total - deleted

df.deleted <- data.frame(deleted = deleted)
df.deleted$tipo <- "eliminado"
colnames(df.deleted) <- c("Cluster","Frecuencia","Tipo")

df.kept <- data.frame(kept = kept)
df.kept$tipo <- "mantenido"
colnames(df.kept) <- c("Cluster","Frecuencia","Tipo")
combined.df <- as.data.frame(rbind(df.deleted, df.kept))

f <- ggplot(combined.df, aes(x=Cluster, y=Frecuencia, fill=Tipo)) 
plot <- f + geom_bar(stat="identity") +
  ggtitle("Celulas consenso") +
  xlab("Cluster") + 
  ylab("Núm celulas")
ggsave(filename=paste0(path2, "barplot_anotacion_eliminados.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "barplot_anotacion_eliminados.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

# FINAL PCA
a <- plotReducedDim((sce.annot), dimred="PCA" , colour_by="tipo.celular", text_by= cluster)+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Anotación Final")
a <- a +  scale_color_viridis_d(name="Tipo celular")
ggsave(filename=paste0(path2, "pca-anotacion.jpeg", sep=""),
       plot=a, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "pca-anotacion.svg", sep=""),
       plot=a, width=10, height=6, units="in", device="svg")


# FINAL TSNE
a <- plotReducedDim((sce.annot), dimred="TSNE120", colour_by="tipo.celular", text_by= cluster)+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Anotación Final") +
  labs(x = "TSNE1", y = "TSNE2")
a <- a +  scale_color_viridis_d(name="Tipo celular")
ggsave(filename=paste0(path2, "tsne-anotacion.jpeg", sep=""),
       plot=a, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "tsne-anotacion.svg", sep=""),
       plot=a, width=10, height=6, units="in", device="svg")

a <- plotReducedDim((sce.annot), dimred="TSNE120", colour_by="Group_sex", text_by= cluster)+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Anotación Final") +
  labs(color = "Tipo celular") +
  labs(x = "TSNE1", y = "TSNE2")
a <- a + scale_color_viridis_d(name="Condición y sexo", labels=c("AD_Mujer", "AD_Hombre",
                                                                 "Control_Mujer", "Control_Hombre"))
ggsave(filename=paste0(path2, "tsne-anotacion_2.jpeg", sep=""),
       plot=a, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "tsne-anotacion_2.svg", sep=""),
       plot=a, width=10, height=6, units="in", device="svg")


# FINAL UMAP
a <- plotReducedDim((sce.annot), dimred="UMAP300", colour_by="tipo.celular", text_by= cluster)+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Anotación Final") +
  labs(color = "Tipo celular") +
  labs(x = "UMAP1", y = "UMAP2")
a <- a +  scale_color_viridis_d(name="Tipo celular")
ggsave(filename=paste0(path2, "umap_anotacion.jpeg", sep=""),
       plot=a, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "umap_anotacion.svg", sep=""),
       plot=a, width=10, height=6, units="in", device="svg")

# ---------------------------- EXTRA PLOTS ----------------------------------

x <- as.data.frame(table(condicion = sce.annot$pathAD, cluster = sce.annot@colData[,cluster]))
f <- ggplot(x, aes(x=cluster, y=Freq, fill=condicion))  
plot <- f + geom_bar(stat="identity") +
  ggtitle("Enfermos vs controles por cluster") +
  xlab("Cluster") + 
  ylab("Núm celulas")
ggsave(filename=paste0(path2, "barplot_enfermos_controles.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "barplot_enfermos_controles.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

x <- as.data.frame(table(condicion = sce.annot$Group_sex, cluster = sce.annot@colData[,cluster]))
f <- ggplot(x, aes(x=cluster, y=Freq, fill=condicion))  
plot <- f + geom_bar(stat="identity") +
  scale_fill_viridis_d(name="Condición y sexo", labels=c("AD_Mujer", "AD_Hombre",
                                                         "Control_Mujer", "Control_Hombre")) +
  ggtitle("Grupos por cluster") +
  xlab("Cluster") + 
  ylab("Núm celulas")
ggsave(filename=paste0(path2, "barplot_grupos_cluster.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "barplot_grupos_cluster.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

x <- as.data.frame(table(condicion = sce.annot$Group_sex, cluster = sce.annot$tipo.celular))
f <- ggplot(x, aes(x=cluster, y=Freq, fill=condicion))  
plot <- f + geom_bar(stat="identity") +
  scale_fill_viridis_d(name="Condición y sexo", labels=c("AD_Mujer", "AD_Hombre",
                                                         "Control_Mujer", "Control_Hombre")) +
  ggtitle("Grupos por tipo celular") +
  xlab("Tipo celular") + 
  ylab("Núm celulas")
ggsave(filename=paste0(path2, "barplot_tipocelular_grupo.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "barplot_tipocelular_grupo.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")


x <- as.data.frame(table(Tipo_celular = sce.annot$tipo.celular, Condicion_Sexo = sce.annot$Group_sex))
x <- x %>%
  group_by(Condicion_Sexo) %>%
  mutate(Freq_relativa = Freq / sum(Freq) * 100)
f <- ggplot(x, aes(x=Condicion_Sexo, y=Freq_relativa, fill=Tipo_celular))  
plot <- f + geom_bar(stat="identity") +
  scale_fill_viridis_d(name="Tipo celular") +
  ggtitle("Tipos celulares por grupo") +
  ylab("Núm celulas") +
  scale_x_discrete(name = "Grupo", labels = c("AD_Mujer", "AD_Hombre",
                                              "Control_Mujer", "Control_Hombre"))
ggsave(filename=paste0(path2, "barplot_grupo_tipocelular.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "barplot_grupo_tipocelular.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")


conteo <- as.data.frame(colData(sce.annot)) %>% 
  group_by(tipo.celular, Group_sex)  %>% 
  dplyr::summarize(count = n())

conteo$CellTypeCondSex = paste0(conteo$tipo.celular,"_",conteo$Group_sex)

axis_spanish_labels <- c(
  "OPC_Control_Male" = "OPC_Control_Hombre",
  "OPC_Control_Female" = "OPC_Control_Mujer",
  "OPC_AD_Male" = "OPC_AD_Hombre",
  "OPC_AD_Female" = "OPC_AD_Mujer",
  "Oligodendrocitos_Control_Male" = "Oligodendrocitos_Control_Hombre",
  "Oligodendrocitos_Control_Female" = "Oligodendrocitos_Control_Mujer",
  "Oligodendrocitos_AD_Male" = "Oligodendrocitos_AD_Hombre",
  "Oligodendrocitos_AD_Female" = "Oligodendrocitos_AD_Mujer",
  "Microglia_Control_Male" = "Microglía_Control_Hombre",
  "Microglia_Control_Female" = "Microglía_Control_Mujer",
  "Microglia_AD_Male" = "Microglía_AD_Hombre",
  "Microglia_AD_Female" = "Microglía_AD_Mujer",
  "Inhibidoras_Control_Male" = "NeuronasInhibidoras_Control_Hombre",
  "Inhibidoras_Control_Female" = "NeuronasInhibidoras_Control_Mujer",
  "Inhibidoras_AD_Male" = "NeuronasInhibidoras_AD_Hombre",
  "Inhibidoras_AD_Female" = "NeuronasInhibidoras_AD_Mujer",
  "Excitadoras_Control_Male" = "NeuronasExcitadoras_Control_Hombre",
  "Excitadoras_Control_Female" = "NeuronasExcitadoras_Control_Mujer",
  "Excitadoras_AD_Male" = "NeuronasExcitadoras_AD_Hombre",
  "Excitadoras_AD_Female" = "NeuronasExcitadoras_AD_Mujer",
  "Astrocitos_Control_Male" = "Astrocitos_Control_Hombre",
  "Astrocitos_Control_Female" = "Astrocitos_Control_Mujer",
  "Astrocitos_AD_Male" = "Astrocitos_AD_Hombre",
  "Astrocitos_AD_Female" = "Astrocitos_AD_Mujer"
)


plot <- ggplot(conteo, aes(x=CellTypeCondSex, y=count, fill = tipo.celular, label = count)) + 
  geom_bar(stat = "identity")+
  theme(axis.text.x = 
          element_text(angle = 90, vjust = 1, hjust = 1)) +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  coord_flip() +
  scale_x_discrete(name="Condición y sexo por tipo celular", labels=axis_spanish_labels) +
  ylab("Número de células") +
  scale_fill_viridis_d(name="Tipo celular")
ggsave(filename=paste0(path2, "barplot_chuli.jpeg", sep=""), 
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "barplot_chuli.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

plot <- ggplot(conteo, aes(x=CellTypeCondSex, y=count, fill = tipo.celular, label = count)) + 
  geom_bar(stat = "identity")+
  theme(axis.text.x = 
          element_text(angle = 45, vjust = 1, hjust = 1)) +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_discrete(name="Tipo celular") +
  scale_x_discrete(name = "Condición y sexo por tipo celular", labels=axis_spanish_labels) +
  ylab("Número de células") +
  scale_fill_viridis_d(name="Tipo celular")
ggsave(filename=paste0(path2, "barplot_chuli_2.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "barplot_chuli_2.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

plot <- ggplot(conteo, aes(x=CellTypeCondSex, y=count, fill = Group_sex, label = count)) + 
  geom_bar(stat = "identity")+
  theme(axis.text.x = 
          element_text(angle = 90, vjust = 1, hjust = 1)) +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  coord_flip() +
  scale_fill_discrete(name="Condición y sexo", labels=c("AD_Mujer", "AD_Hombre",
                                                        "Control_Mujer", "Control_Hombre")) +
  scale_x_discrete(name="Condición y sexo por tipo celular") +
  ylab("Número de células") +
  scale_fill_viridis_d(name="Condición y sexo", labels=axis_spanish_labels)
ggsave(filename=paste0(path2, "barplot_chuli_3.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "barplot_chuli_3.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

plot <- ggplot(conteo, aes(x=CellTypeCondSex, y=count, fill = Group_sex, label = count)) + 
  geom_bar(stat = "identity")+
  theme(axis.text.x = 
          element_text(angle = 45, vjust = 1, hjust = 1)) +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_discrete(name="Condición y sexo", labels=c("AD_Mujer", "AD_Hombre",
                                                        "Control_Mujer", "Control_Hombre")) +
  scale_x_discrete(name = "Condición y sexo por tipo celular", labels=axis_spanish_labels) +
  ylab("Número de células") +
  scale_fill_viridis_d(name="Condición y sexo")
ggsave(filename=paste0(path2, "barplot_chuli_4.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "barplot_chuli_4.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

saveRDS(sce.annot, file = "processed_data/scegv_dw_clusters_filtered_annot2.rds")

x <- as.data.frame(table(sce.annot$tipo.celular))
y <- as.vector(x$Freq)
names(y) <- x$Var1
z <- (y*100)/sum(y)

proporciones <- data.frame(tipo.celular= names(z),
                           porcentaje = z,
                           total = y)

require(openxlsx)
write.xlsx(proporciones, file = paste0("processed_data/abundancia_celular_gv_dw.xlsx"),rowNames=TRUE)


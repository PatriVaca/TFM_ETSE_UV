# Single cell RNA-seq data analysis
# 7.FilteringCellTypeAnnotationGV_DW
# Author: Patricia del Carmen Vaca Rubio
# Date: 05/25

pacman::p_load(SingleCellExperiment, scater, scran, dplyr, BiocParallel, 
               scDblFinder, AnnotationHub, tidyverse, patchwork, ggvenn, broom,
               kableExtra,HDF5Array, viridis, DropletUtils,
               devtools, celda, ggvenn, bluster, cluster, ggpubr, SCINA, SingleR,
               CHETAH, pheatmap, ggplot2)

setwd("/clinicfs/projects/i63/tfm_hipocampo/GSE163577/")
system("mkdir -p figures/7.FilteredCellTypeAnnotation/7.2.ClusteringGV_dw/prefiltering")
system("mkdir -p figures/7.FilteredCellTypeAnnotation/7.2.ClusteringGV_dw/filtered")

scegv_dw_clusters <- readRDS("processed_data/scegv_dw_clusters_annot1.rds")

# Obtain the cell types
cell_types <- unique(scegv_dw_clusters$SCINA_label)

# Assign the default viridis colors
personal_colors <- viridis::viridis(length(cell_types))
names(personal_colors) <- cell_types

# Modify the color of epithelial cells
personal_colors["Celulas_epiteliales"] <- "#F77FBE" 

# ---------------------------- Prefiltering ----------------------------------

# CONSENSUAL RESULTS GENERATION #
path1 <- "figures/7.FilteredCellTypeAnnotation/7.2.ClusteringGV_dw/prefiltering/"

sce.annot  <- scegv_dw_clusters

# ---------------------------- Prefiltering ----------------------------------

# FINAL RESULTS BEFORE FILTERING #
# Select a cluster
cluster <- "cluster40"

tab <- table(label = sce.annot$SCINA_label, cluster = colData(sce.annot)$cluster40) 
y  <- as.data.frame(tab)

# CELL TYPES PER CLUSTER
f  <- ggplot(y, aes(x=cluster, y=Freq, fill=label))
plot <- f + 
  geom_bar(stat = "identity") +
  scale_fill_manual(name="Tipo celular", values=personal_colors) +
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
  scale_fill_manual(name="Tipo celular", values=personal_colors) +
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
a  <- plotReducedDim((sce.annot), dimred="PCA" , colour_by="SCINA_label", text_by= cluster)+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Anotación Final") 
a  <- a +  scale_color_manual(name="Tipo celular", values=personal_colors)
ggsave(filename=paste0(path1, "pca_preanotacion.jpeg", sep=""),
       plot=a, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path1, "pca_preanotacion.svg", sep=""),
       plot=a, width=10, height=6, units="in", device="svg")


# FINAL TSNE
a  <- plotReducedDim((sce.annot), dimred="TSNE120", colour_by="SCINA_label", text_by= cluster)+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Anotación Final") +
  labs(x = "TSNE1", y = "TSNE2")
a <- a +  scale_color_manual(name="Tipo celular", values=personal_colors)
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
a <- plotReducedDim((sce.annot), dimred="UMAP300", colour_by="SCINA_label", text_by= cluster)+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Anotación Final")  +
  labs(x = "UMAP1", y = "UMAP2")
a <- a +  scale_color_manual(name="Tipo celular", values=personal_colors)
ggsave(filename=paste0(path1, "umap_preanotacion.jpeg", sep=""),
       plot=a, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path1, "umap_preanotacion.svg", sep=""),
       plot=a, width=10, height=6, units="in", device="svg")

# Check if the "unknown" cell types are disease exclusive # 
# Extract UMAP coordinates and variables
umap_data <- as.data.frame(reducedDim(sce.annot, "UMAP300"))
umap_data$Group_sex <- colData(sce.annot)$Group_sex
umap_data$cell_type <- colData(sce.annot)$SCINA_label

a <- ggplot(umap_data, aes(x = V1, y = V2, color = cell_type)) +
  geom_point() +
  facet_wrap(~ Group_sex,
             labeller = labeller(Group_sex = c(
               "AD_F" = "AD_Mujer",
               "AD_M"   = "AD_Hombre",
               "Control_F" = "Control_Mujer",
               "Control_M"   = "Control_Hombre"
             ))) +
  theme(legend.key.size = unit(0.25, 'cm'))
a <- a +  scale_color_manual(name="Tipo celular", values=personal_colors)
ggsave(filename=paste0(path1, "disease_cell_unknown.jpeg", sep=""),
       plot=a, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path1, "disease_cell_unknown.svg", sep=""),
       plot=a, width=10, height=6, units="in", device="svg")

# ---------------------------- Filtering ----------------------------------
table(sce.annot$cluster40, sce.annot$SCINA_label)
sce.annot =  sce.annot[,!sce.annot@colData[,"SCINA_label"] == "unknown"]

table(sce.annot$cluster40, sce.annot$SCINA_label)
# Manual filtering to remove the outliers of each cluster
sce.annot <- sce.annot[,!(sce.annot$SCINA_label != "Celulas_epiteliales" & sce.annot$cluster40 == 1)]
sce.annot <- sce.annot[,!(sce.annot$SCINA_label != "Celulas_epiteliales" & sce.annot$cluster40 == 2)]
sce.annot <- sce.annot[,!(sce.annot$SCINA_label != "Excitadoras" & sce.annot$cluster40 == 3)]
sce.annot <- sce.annot[,!(sce.annot$SCINA_label != "Inhibidoras" & sce.annot$cluster40 == 4)]
sce.annot <- sce.annot[,!(sce.annot$SCINA_label != "Astrocitos" & sce.annot$cluster40 == 5)]
sce.annot <- sce.annot[,!(sce.annot$SCINA_label != "Excitadoras" & sce.annot$cluster40 == 6)]
sce.annot <- sce.annot[,!(sce.annot$SCINA_label != "Oligodendrocitos" & sce.annot$cluster40 == 7)]
sce.annot <- sce.annot[,!(sce.annot$SCINA_label != "Astrocitos" & sce.annot$cluster40 == 8)]
sce.annot <- sce.annot[,!(sce.annot$SCINA_label != "Oligodendrocitos" & sce.annot$cluster40 == 9)]
sce.annot <- sce.annot[,!(sce.annot$SCINA_label != "Microglia" & sce.annot$cluster40 == 10)]
sce.annot <- sce.annot[,!(sce.annot$SCINA_label != "Celulas_epiteliales" & sce.annot$cluster40 == 11)]
sce.annot <- sce.annot[,!(sce.annot$SCINA_label != "Oligodendrocitos" & sce.annot$cluster40 == 12)]
sce.annot <- sce.annot[,!(sce.annot$SCINA_label != "Inhibidoras" & sce.annot$cluster40 == 13)]
sce.annot <- sce.annot[,!(sce.annot$SCINA_label != "OPC" & sce.annot$cluster40 == 14)]
sce.annot <- sce.annot[,!(sce.annot$SCINA_label != "Celulas_epiteliales" & sce.annot$cluster40 == 15)]
sce.annot <- sce.annot[,!(sce.annot$SCINA_label != "Celulas_epiteliales" & sce.annot$cluster40 == 16)]

table(sce.annot$cluster40, sce.annot$SCINA_label)

as.data.frame(table(sce.annot$cluster40, sce.annot$SCINA_label)) %>% 
  group_by(Var2)  %>% 
  summarise(cond_disp = sum(Freq))

path2 <- "figures/7.FilteredCellTypeAnnotation/7.2.ClusteringGV_dw/filtered/"

tab <- table(label=sce.annot$SCINA_label, cluster=sce.annot$cluster40) 
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
  scale_fill_manual(name="Tipo celular", values=personal_colors) +
  ggtitle("Num cel por tipo celular") +
  xlab("Cluster") + 
  ylab("Núm celulas")
ggsave(filename=paste0(path2, "barplot_anotacion.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "barplot_anotacion.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

# DELETED CELLS BARPLOT
delete_unknown <- grep("unknown", sce.annot$SCINA_label)
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
a <- plotReducedDim((sce.annot), dimred="PCA" , colour_by="SCINA_label", text_by= cluster)+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Anotación Final")
a <- a +  scale_color_manual(name="Tipo celular", values=personal_colors)
ggsave(filename=paste0(path2, "pca-anotacion.jpeg", sep=""),
       plot=a, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "pca-anotacion.svg", sep=""),
       plot=a, width=10, height=6, units="in", device="svg")


# FINAL TSNE
a <- plotReducedDim((sce.annot), dimred="TSNE120", colour_by="SCINA_label", text_by= cluster)+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Anotación Final") +
  labs(x = "TSNE1", y = "TSNE2")
a <- a +  scale_color_manual(name="Tipo celular", values=personal_colors)
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
a <- plotReducedDim((sce.annot), dimred="UMAP300", colour_by="SCINA_label", text_by= cluster)+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Anotación Final") +
  labs(color = "Tipo celular") +
  labs(x = "UMAP1", y = "UMAP2")
a <- a +  scale_color_manual(name="Tipo celular", values=personal_colors)
ggsave(filename=paste0(path2, "umap_anotacion.jpeg", sep=""),
       plot=a, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "umap_anotacion.svg", sep=""),
       plot=a, width=10, height=6, units="in", device="svg")

# ---------------------------- EXTRA PLOTS ----------------------------------

x <- as.data.frame(table(condicion = sce.annot$Clinical.AD.or.not, cluster = sce.annot@colData[,cluster]))
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

x <- as.data.frame(table(condicion = sce.annot$Group_sex, cluster = sce.annot$SCINA_label))
f <- ggplot(x, aes(x=cluster, y=Freq, fill=condicion))  
plot <- f + geom_bar(stat="identity") +
  scale_fill_viridis_d(name="Condición y sexo", labels=c("AD_Mujer", "AD_Hombre",
                                                         "Control_Mujer", "Control_Hombre")) +
  ggtitle("Grupos por tipo celular") +
  xlab("Tipo celular") + # Volver a ejecutar para tener las etiquetas en el eje x
  ylab("Núm celulas")
ggsave(filename=paste0(path2, "barplot_tipocelular_grupo.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "barplot_tipocelular_grupo.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")


x <- as.data.frame(table(Tipo_celular = sce.annot$SCINA_label, Condicion_Sexo = sce.annot$Group_sex))
x <- x %>%
  group_by(Condicion_Sexo) %>%
  mutate(Freq_relativa = Freq / sum(Freq) * 100)
f <- ggplot(x, aes(x=Condicion_Sexo, y=Freq_relativa, fill=Tipo_celular))  
plot <- f + geom_bar(stat="identity") +
  scale_fill_manual(name="Tipo celular", values=personal_colors) + # He cambiado scale_color_manual por scale_color_fill (hay que volver a ejecutarlo)
  ggtitle("Tipos celulares por grupo") +
  ylab("Núm celulas") +
  scale_x_discrete(name = "Grupo", labels = c("AD_Mujer", "AD_Hombre",
                                              "Control_Mujer", "Control_Hombre"))
ggsave(filename=paste0(path2, "barplot_grupo_tipocelular.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "barplot_grupo_tipocelular.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")


conteo <- as.data.frame(colData(sce.annot)) %>% 
  group_by(SCINA_label, Group_sex)  %>% 
  dplyr::summarize(count = n())

conteo$CellTypeCondSex = paste0(conteo$SCINA_label,"_",conteo$Group_sex)

axis_spanish_labels <- c(
  "Astrocitos_AD_F" = "Astrocitos_AD_Mujer",
  "Astrocitos_AD_M" = "Astrocitos_AD_Hombre",
  "Astrocitos_Control_F" = "Astrocitos_Control_Mujer",
  "Astrocitos_Control_M" = "Astrocitos_Control_Hombre",
  "Celulas_epiteliales_AD_F" = "CélulasEpiteliales_AD_Mujer",
  "Celulas_epiteliales_AD_M" = "CélulasEpiteliales_AD_Hombre",
  "Celulas_epiteliales_Control_F" = "CélulasEpiteliales_Control_Mujer",
  "Celulas_epiteliales_Control_M" = "CélulasEpiteliales_Control_Hombre",
  "Microglia_AD_F" = "Microglía_AD_Mujer",
  "Microglia_AD_M" = "Microglía_AD_Hombre",
  "Microglia_Control_F" = "Microglía_Control_Mujer",
  "Microglia_Control_M" = "Microglía_Control_Hombre",
  "Excitadoras_AD_F" = "NeuronasExcitadoras_AD_Mujer",
  "Excitadoras_AD_M" = "NeuronasExcitadoras_AD_Hombre",
  "Excitadoras_Control_F" = "NeuronasExcitadoras_Control_Mujer",
  "Excitadoras_Control_M" = "NeuronasExcitadoras_Control_Hombre",
  "Inhibidoras_AD_F" = "NeuronasInhibidoras_AD_Mujer",
  "Inhibidoras_AD_M" = "NeuronasInhibidoras_AD_Hombre",
  "Inhibidoras_Control_F" = "NeuronasInhibidoras_Control_Mujer",
  "Inhibidoras_Control_M" = "NeuronasInhibidoras_Control_Hombre",
  "Oligodendrocitos_AD_F" = "Oligodendrocitos_AD_Mujer",
  "Oligodendrocitos_AD_M" = "Oligodendrocitos_AD_Hombre",
  "Oligodendrocitos_Control_F" = "Oligodendrocitos_Control_Mujer",
  "Oligodendrocitos_Control_M" = "Oligodendrocitos_Control_Hombre",
  "OPC_AD_F" = "OPC_AD_Mujer",
  "OPC_AD_M" = "OPC_AD_Hombre",
  "OPC_Control_F" = "OPC_Control_Mujer",
  "OPC_Control_M" = "OPC_Control_Hombre"
)

plot <- ggplot(conteo, aes(x=CellTypeCondSex, y=count, fill = SCINA_label, label = count)) + 
  geom_bar(stat = "identity")+
  theme(axis.text.x = 
          element_text(angle = 90, vjust = 1, hjust = 1)) +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  coord_flip() +
  scale_x_discrete(name="Condición y sexo por tipo celular", labels=axis_spanish_labels) +
  ylab("Número de células") +
  scale_fill_manual(name="Tipo celular", values=personal_colors)
ggsave(filename=paste0(path2, "barplot_chuli.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(path2, "barplot_chuli.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")

plot <- ggplot(conteo, aes(x=CellTypeCondSex, y=count, fill = SCINA_label, label = count)) + 
  geom_bar(stat = "identity")+
  theme(axis.text.x = 
          element_text(angle = 45, vjust = 1, hjust = 1)) +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_discrete(name="Tipo celular") +
  scale_x_discrete(name = "Condición y sexo por tipo celular", labels=axis_spanish_labels) +
  ylab("Número de células") +
  scale_fill_manual(name="Tipo celular", values=personal_colors)
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
  scale_x_discrete(name="Condición y sexo por tipo celular", labels=axis_spanish_labels) +
  ylab("Número de células") +
  scale_fill_viridis_d(name="Condición y sexo")
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

saveRDS(sce.annot, file = "processed_data/scegv_dw_clusters_filtered_annot1.rds")

x <- as.data.frame(table(sce.annot$SCINA_label))
y <- as.vector(x$Freq)
names(y) <- x$Var1
z <- (y*100)/sum(y)

proporciones <- data.frame(tipo.celular = names(z),
                           porcentaje = z,
                           total = y)

require(openxlsx)
write.xlsx(proporciones, file = paste0("processed_data/abundancia_celular_gv_dw.xlsx"),rowNames=TRUE)






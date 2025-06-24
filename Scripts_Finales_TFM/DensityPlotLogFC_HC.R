# Single cell RNA-seq data analysis
# Density plot
# To set the logFC threshold per cell type after the meta-analysis
# Author: Patricia del Carmen Vaca Rubio
# Date: 05/25

setwd("/clinicfs/projects/i63/tfm_hipocampo/")

cell_types <- c("Astrocitos",	"Oligodendrocitos",
                "OPC", "Excitadoras",	"Inhibidoras")

output_dir <- "Metaanalysis_HC/"
df_list <- list()

for (cell_type in cell_types) {
  sig.genes.df <- read.delim(paste0(output_dir, cell_type, "/FM/sig.genes.tsv"))
  sig.genes.df$CellType <- cell_type
  df_list[[cell_type]] <- sig.genes.df[,c("ENSEMBL", "logFC", "CellType")]
  
  logfc_df <- do.call(rbind, df_list)
}

plot <- ggplot(data=logfc_df, aes(x=abs(logFC), group=CellType, fill=CellType)) +
  geom_vline(xintercept = 0.5) +
  geom_density(adjust=1.5, alpha=.4) 
ggsave(filename=paste0(output_dir, "density_logfc.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(output_dir, "density_logfc.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")


plot <- ggplot(data=logfc_df, aes(x=abs(logFC), group=CellType, fill=CellType)) +
  geom_density(adjust=1.5) +
  facet_wrap(~CellType) +
  geom_vline(xintercept = 0.5) +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.ticks.x=element_blank()
  )
ggsave(filename=paste0(output_dir, "grid_density_logfc.jpeg", sep=""),
       plot=plot, width=10, height=6, units="in", device="jpeg", dpi=300)
ggsave(filename=paste0(output_dir, "grid_density_logfc.svg", sep=""),
       plot=plot, width=10, height=6, units="in", device="svg")



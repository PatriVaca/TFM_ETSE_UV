# Single cell RNA-seq data analysis
# Meta-analysis significant results filtering
# Only the significant genes (FDR < 0.05) for each cell type
# that have an abs(logFC) > 0.5 are kept for later analysis with STRING.
# Author: Patricia del Carmen Vaca Rubio
# Date: 05/25

setwd("/clinicfs/projects/i63/tfm_hipocampo/")

cell_types <- c("Microglia", "Astrocitos",	"Oligodendrocitos",
                "OPC", "Excitadoras",	"Inhibidoras")

output_dir <- "Metaanalysis_HC/"

for (cell_type in cell_types) {
   for (c in contrasts_vector) {
    res_F <- read.delim(paste0(output_dir, cell_type, "/F/sig.genes.tsv"))
    res_M <- read.delim(paste0(output_dir, cell_type, "/M/sig.genes.tsv"))
    res_FM <- read.delim(paste0(output_dir, cell_type, "/FM/sig.genes.tsv"))
    res_FM_allgenes <- read.delim(paste0(output_dir, cell_type, "/FM/all.genes.tsv"))
    
    res_F_filtered <- res_F %>% 
      filter(abs(logFC) > 0.5) %>% 
      arrange(p.adjust.fdr)  # Order by significance
    
    res_M_filtered <- res_M %>% 
      filter(abs(logFC) > 0.5) %>% 
      arrange(p.adjust.fdr)  
    
    res_FM_filtered <- res_FM %>% 
      filter(abs(logFC) > 0.5) %>% 
      arrange(p.adjust.fdr)
    
    res_FM_allgenes_filtered <- res_FM_allgenes %>% 
      filter(abs(logFC) > 0.5) %>% 
      arrange(p.adjust.fdr) 
    
    # Save results
    write.table(res_F_filtered, paste0(output_dir, cell_type, "/F/sig.genes_filtered.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(res_M_filtered, paste0(output_dir, cell_type, "/M/sig.genes_filtered.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(res_FM_filtered, paste0(output_dir, cell_type, "/FM/sig.genes_filtered.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(res_FM_allgenes_filtered, paste0(output_dir, cell_type, "/FM/all.genes_filtered.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
   }
}

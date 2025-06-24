# Single cell RNA-seq data analysis
# Functional Enrichment
# Author: Patricia del Carmen Vaca Rubio
# Date: 06/25

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
# library(biomaRt)

#---------------------- Load data and perform the enrichment ----------------------------

setwd("/clinicfs/projects/i63/tfm_hipocampo/Metaanalysis_HC/HC")
dir <- getwd()
cell_types <- list.dirs(path = dir, full.names = FALSE, recursive = FALSE)
results_GO <- list()
results_KEGG <- list()

for (cell_type in cell_types) {
  contrasts <- list.dirs(path = file.path(dir, cell_type), full.names = FALSE, recursive = FALSE)
  
  for (contrast in contrasts) {
    dir.tsv <- file.path(dir, cell_type, contrast, 'sig.genes.tsv')
    
    if (file.exists(dir.tsv)) {
      tryCatch({
        # Read data
        df <- read.delim(dir.tsv, header = TRUE, sep = '\t')
        
        # We want the log2 fold change 
        original_gene_list <- df$logFC
        
        # Name the vector
        names(original_gene_list) <- df$ENSEMBL
        
        # Omit any NA values 
        gene_list <- na.omit(original_gene_list)
        
        # Sort the list in decreasing order (required for clusterProfiler)
        gene_list <- sort(gene_list, decreasing = TRUE)
        
        #---------------------- GO ---------------------------------------------
        
        gseG <- gseGO(
          geneList = gene_list,
          ont = "BP",
          keyType = "ENSEMBL",
          OrgDb = 'org.Hs.eg.db',
          minGSSize = 3,
          maxGSSize = 800,
          pvalueCutoff = 0.05,
          pAdjustMethod = "none"
        )
        
        # Save GO enrichment
        if (nrow(gseG) > 0) {
          results_GO[[paste0(cell_type, "_", contrast)]] <- gseG
        }
        
        #---------------------- KEGG -------------------------------------------
        
        # Convert gene IDs for gseKEGG function
        # We will lose some genes here because not all IDs will be converted
        ids <- bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = 'org.Hs.eg.db')
        # keytypes(org.Hs.eg.db)  # To verify that "ENSEMBL" and "ENTREZID" are available

        # Remove duplicate IDS
        dedup_ids <- ids[!duplicated(ids$ENSEMBL), ]

        # Create a new dataframe df2 which has only the genes which were
        # successfully mapped using the bitr function above
        df2 <- df[df$ENSEMBL %in% dedup_ids$ENSEMBL, ]

        # Create a new column in df2 with the corresponding ENTREZ IDs
        df2$ENTREZID <- dedup_ids$ENTREZID

        # We want the log2 fold change
        kegg_gene_list <- df2$logFC

        # Name the vector
        names(kegg_gene_list) <- df2$ENTREZID

        # Omit any NA values and sort the list in decreasing order
        kegg_gene_list <- na.omit(sort(kegg_gene_list, decreasing = TRUE))
        
        ## Another alternative because the entrezid genes obtained via bitr are not
        # recognized by gseKEGG
        # Use last version of biomart
        # mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = NULL)
        
        # # Connect to biomart Ensembl version 108
        # ensembl <- useEnsembl(
        #   biomart = "ensembl",
        #   dataset = "hsapiens_gene_ensembl",
        #   version = 108
        # )
        # 
        # # Get Entrez IDs from Ensembl IDs
        # ids <- getBM(
        #   attributes = c("ensembl_gene_id", "entrezgene_id"),
        #   filters = "ensembl_gene_id",
        #   values = names(original_gene_list),
        #   mart = ensembl  # Replace with your desired Ensembl version
        # )
        # 
        # # Remove duplicate IDS
        # dedup_ids <- ids[!duplicated(ids$ensembl_gene_id), ]
        # 
        # # Create a new dataframe df2 which has only the genes which were
        # # successfully mapped using the useEnsembl function above
        # df2 <- df[df$ENSEMBL %in% dedup_ids$ensembl_gene_id, ]
        # 
        # # Create a new column in df2 with the corresponding ENTREZ IDs
        # df2$ENTREZID <- dedup_ids$entrezgene_id
        # 
        # # We want the log2 fold change 
        # kegg_gene_list <- df2$logFC
        # 
        # # Name the vector
        # names(kegg_gene_list) <- df2$ENTREZID
        # names(kegg_gene_list) <- as.character(names(kegg_gene_list))
        # 
        # # Omit any NA values and sort the list in decreasing order
        # kegg_gene_list <- na.omit(sort(kegg_gene_list, decreasing = TRUE))
        # 
        # ## Several trials to avoid gseKEGG lack of results
        # # Map valid genes
        # valid_ids <- bitr_kegg(names(kegg_gene_list), fromType = "ncbi-geneid", toType = "kegg", organism = "hsa")
        # 
        # # Filter kegg_gene_list to keep only the valid IDs
        # kegg_gene_list_filtered <- kegg_gene_list[names(kegg_gene_list) %in% valid_ids$`ncbi-geneid`]
        # 
        # # Secure names as characters
        # names(kegg_gene_list_filtered) <- as.character(names(kegg_gene_list_filtered))
        # In case of biomaRt way, use kegg_gene_list_filtered instead of kegg_gene_list
        # in geneList parameter of gseKEGG.
        
        gseK <- gseKEGG(
          geneList = kegg_gene_list,
          organism = "hsa",
          keyType = "ncbi-geneid",
          minGSSize = 3,
          maxGSSize = 800,
          pvalueCutoff = 0.05,
          pAdjustMethod = "none"
        )
        
        # Save KEGG enrichment
        if (nrow(gseK) > 0) {
          results_KEGG[[paste0(cell_type, "_", contrast)]] <- gseK
        }
        
      }, error = function(e) {
        message(paste("Error en", cell_type, contrast, ":", e$message))
      })
    }
  }
}

#---------------------- Create plots -----------------------

for (result_name in names(results_GO)) {
  cell_type <- strsplit(result_name, "_")[[1]][1]
  contrast <- sub(paste0(cell_type, "_"), "", result_name)
  
  output_dir <- file.path(dir, cell_type, contrast, 'FunctionalEnrichment')
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  gseG <- results_GO[[result_name]]
  
  if (!is.null(gseG) && nrow(gseG) > 0){
    tryCatch({
      # Dotplot GO
      dotplotGO <- dotplot(gseG, showCategory = 9, split = ".sign", title = cell_type) + 
        facet_grid(.~.sign) +
        theme(axis.text.y = element_text(size = 8)) +
        scale_fill_gradient(low = 'dodgerblue1', high = 'salmon1')
      
      ggsave(file.path(output_dir, 'dotplot.GO.jpeg'), dotplotGO, width = 10, height = 6, dpi = 300)
      ggsave(file.path(output_dir, 'dotplot.GO.svg'), dotplotGO, width = 10, height = 6)
      
      # Enrichmap GO
      gseG <- pairwise_termsim(gseG)
      enrichmapGO <- emapplot(gseG, showCategory = 10, cex_label_category = 0.6, title = cell_type)+
        scale_fill_gradient(low = 'dodgerblue1', high = 'salmon1')
      
      jpeg(file.path(output_dir, 'enrichmap.GO.jpeg'), 
           width = 10, height = 6, units = 'in', res = 300)
      print(enrichmapGO)
      dev.off()
      svg(file.path(output_dir, 'enrichmap.GO.svg'), 
          width = 10, height = 6)
      print(enrichmapGO)
      dev.off()
      
      # Ridge plot GO
      ridgeplotGO <- ridgeplot(gseG) + labs(x = "Perfil de enrequecimiento") +
        theme(
          axis.text.y = element_text(size = 6)) +
        scale_fill_gradient(low = 'dodgerblue1', high = 'salmon1')
      
      jpeg(file.path(output_dir, 'ridgeplot.GO.jpeg'), 
           width = 10, height = 6, units = 'in', res = 300)
      print(ridgeplotGO)
      dev.off()
      svg(file.path(output_dir, 'ridgeplot.GO.svg'), 
          width = 10, height = 6)
      print(ridgeplotGO)
      dev.off()
    }, error = function(e) message("Ha ocurrido un error para ", cell_type, ": ", e$message))
  }
  
  gseK <- results_KEGG[[result_name]]
  
  if (!is.null(gseK) && nrow(gseK) > 0){
    # Dotplot KEGG
    dotplotKEGG <- dotplot(gseK, showCategory = 10, split = ".sign", title = cell_type) + 
      facet_grid(.~.sign) +
      theme(axis.text.y = element_text(size = 8))+
      scale_fill_gradient(low = 'dodgerblue1', high = 'salmon1')
    
    ggsave(file.path(output_dir, 'dotplot.KEGG.jpeg'), dotplotKEGG, width = 10, height = 6, dpi = 300)
    ggsave(file.path(output_dir, 'dotplot.KEGG.svg'), dotplotKEGG, width = 10, height = 6)
    
    # Enrichmap KEGG
    gseK <- pairwise_termsim(gseK)
    enrichmapKEGG <- emapplot(gseK, showCategory = 10, cex_label_category = 0.6, title = cell_type)+
      scale_fill_gradient(low = 'dodgerblue1', high = 'salmon1')
    
    jpeg(file.path(output_dir, 'enrichmap.KEGG.jpeg'), 
         width = 10, height = 6, units = 'in', res = 300)
    print(enrichmapKEGG)
    dev.off()
    svg(file.path(output_dir, 'enrichmap.KEGG.svg'), 
        width = 10, height = 6)
    print(enrichmapKEGG)
    dev.off()
    
    # Ridge plot KEGG
    ridgeplotKEGG <- ridgeplot(gseK) + labs(x = "Perfil de enrequecimiento") +
      theme(
        axis.text.y = element_text(size = 6)) +
      scale_fill_gradient(low = 'dodgerblue1', high = 'salmon1')
    
    jpeg(file.path(output_dir, 'ridgeplot.KEGG.jpeg'), 
         width = 10, height = 6, units = 'in', res = 300)
    print(ridgeplotKEGG)
    dev.off()
    svg(file.path(output_dir, 'ridgeplot.KEGG.svg'), 
        width = 10, height = 6)
    print(ridgeplotKEGG)
    dev.off()
    
    #GSEA Plot
    plotverde <- gseaplot(gseK, by = "all", title = gseK$Description[1], geneSetID = 1)
    jpeg(file.path(output_dir, 'plotverde.KEGG.jpeg'), 
         width = 10, height = 6, units = 'in', res = 300)
    print(plotverde)
    dev.off()
    svg(file.path(output_dir, 'plotverde.KEGG.svg'), 
        width = 10, height = 6)
    print(plotverde)
    dev.off()
  }
}


######################################
### metagene.R
### 14/06/2020, fgarcia@cipf.es
### 22/09/2020, jfcatala@cipf.es (adapted for MS output)
### 11/12/2020, jfcatala@cipf.es (RNA-Seq with DESeq2)
### 05/05/2022, jlleraoy@gmail.com  (adapted for microRNA)
### 23/02/2023, fernandogorgonz8@gmail.com  (adapted for sc-RNAseq)
######################################

# Single cell RNA-seq data analysis
# Meta-analysis
# Author: Patricia del Carmen Vaca Rubio
# Date: 05/25

#----------------------------- Load packages -------------------------------####

library(Biobase)
library(metafor)
library(tibble)
library(limma)
library(dplyr)
library(ggplot2)

#-------------------------------- Arguments --------------------------------####

setwd("/clinicfs/projects/i63/tfm_hipocampo/")

cell_types <- c("Microglia", "Astrocitos",	"Oligodendrocitos",
                "OPC", "Excitadoras",	"Inhibidoras")
contrasts_vector <- c("F", "M", "FM")

output_dir <- "Metaanalysis_HC/"
dir.create(output_dir)

for (type in cell_types) {
  for (c in contrasts_vector) {
    cat("\nProcesando:", type, "- Contraste:", c, "\n")
    args <- list(dir = "/processed_data/8.DEG/", # MAST output folder in each study
                 contrast = c,
                 cell = type,
                 tissue = "HC",
                 metadata = "studies_HC.tsv", 
                 pvalue = 0.05, 
                 outputdir = output_dir)
    
    metastudies <- read.delim(args$metadata)
    stopifnot("Not supported contrast." = args$contrast %in% contrasts_vector,
              "Tissue not among possibles
              (Tissue column from metadata file)." = args$tissue %in% metastudies$Tissue,
              "Not all directories exist. 
              Check --dir param or metadata file." = {unlist(lapply(metastudies$Studies, FUN=function(x){dir.exists(paste0(x,args$dir))}))})
    
    # Keep just the studies with the corresponding cell type
    filtering <- metastudies
    filtering <- filtering[filtering[,args$cell],]
    metastudies <- filtering[,c(1,2)]
    
    #---------------------- Make directories & Load data -----------------------####
    
    metaanalysis_name <- paste("MA", args$tissue, args$cell, args$contrast, sep="_")
    print(metaanalysis_name)
    meta_dir <- paste0(getwd(),"/",args$outputdir ,args$tissue,"/",args$cell ,"/",args$contrast)
    dir.create(meta_dir, recursive = TRUE)
    
    # Load data
    selected_studies <- metastudies$Studies[metastudies$Tissue==args$tissue]
    MAST_results_list=list()
    for (i in selected_studies){
      print(paste("Cargando datos del estudio", i))
      MAST_results_list[[i]]=read.delim(paste0(i, args$dir, args$cell ,"/",args$cell,"res_",args$contrast,".tsv"), row.names=1)
    }
    
    setwd(meta_dir)
    
    #---------------------------- Meta-analysis --------------------------------####
    
    # STEP 0. Pre-processing previous data ####
    # ===============================================================
    ## Calculate SE
    # OPCION 1: (sum[, "ci.hi"] - sum[, "ci.lo"])/ 3.92 <--
    # OPCION 2: z/logFC
    
    # STEP 1. Preparing input for meta-analysis: LOR and SE matrix ####
    # ===============================================================
    
    # Search a list including all unique ID genes for all studies
    genes <- unlist(lapply(MAST_results_list, rownames))
    genes <- unique (genes)
    genes <- sort (genes)
    
    
    # Generate matrix with all SE (Standard Error) for all studies
    mat.SE <- matrix (NA, nrow = length (genes), ncol = length(MAST_results_list))
    rownames (mat.SE) <- genes
    colnames (mat.SE) <- names(MAST_results_list)
    
    for (i in 1:length(MAST_results_list)){
      mat.SE[, names(MAST_results_list[i])] <- MAST_results_list[[i]][rownames(mat.SE), "SE"]
    }
    mat.SE[mat.SE==0]=NA
    cat("\nDimensiones originales de la matriz SE:")
    print(dim(mat.SE))
    
    # Select genes included at least in 2 or more studies
    mat.SE.NA <- is.na(mat.SE)
    sum.NA <-  apply(mat.SE.NA, 1, sum)
    min.sum.NA <- sum.NA < ncol(mat.SE) - 1
    
    # Filter by min.sum.NA
    mat.SE <- mat.SE[min.sum.NA == T, ]
    mat.SE.NA <- mat.SE.NA[min.sum.NA == T, ]
    
    cat("\nDimensiones tras filtrar la matriz SE:")
    print(dim(mat.SE))
    
    # Generate matrix with all logFC for all studies
    mat.logFC <- matrix (NA, nrow = length(genes), ncol = length(MAST_results_list))
    rownames (mat.logFC) <- genes
    colnames (mat.logFC) <- names(MAST_results_list)
    
    for (i in 1:length(MAST_results_list)){
      mat.logFC[, names(MAST_results_list[i])] <- MAST_results_list[[i]][rownames(mat.logFC), "logFC"]
    }
    
    cat("\nDimensiones originales de la matriz logFC:")
    print(dim(mat.logFC))
    
    # Filter by min.sum.NA
    mat.logFC <- mat.logFC[min.sum.NA == T, ]
    
    cat("\nDimensiones tras filtrar de la matriz logFC:")
    print(dim(mat.logFC))
    
    # STEP 2. Meta-analysis for genes ####
    # ===============================================================
    
    # suppose between-study variance is non-zero.
    # there are different methods to estimate this variance:
    # DL (Dersimonian-Laird), REML (Restricted maximum-likelihood, default)....
    # Now we have logFC and SE  (not VARIANCE), so:
    # yi -> logFC   sei -> SE
    # result.lor <- rma(yi = mat.logFC[1, ],
    #                   sei = mat.SE[1, ],   #pay attention, not vi (variance)
    #                   method = "DL") # DerSimonian-Laird.
    
    # explore the function to do the meta-analysis
    #?rma
    
    MA <- lapply(1:length(rownames(mat.logFC)),
                 function(x){rma(yi = mat.logFC[x, ],
                                 sei = mat.SE[x, ],
                                 method = "DL")})
    
    names (MA) <- rownames(mat.logFC)
    
    print(paste0("Dimensiones del metaanalisis: ", length(MA)))
    
    #data.frame including all detailed results:
    result_meta <- as.data.frame(do.call("rbind",
                                         lapply(MA,
                                                function(x){
                                                  c(x$ci.lb, x$b, x$ci.ub,
                                                    x$pval, x$QE, x$QEp, x$se,
                                                    x$tau2, x$I2, x$H2)
                                                })))
    
    colnames(result_meta) <- c("lower_bound", "logFC", "upper_bound",
                               "pvalue", "QE", "QEp", "SE", "tau2", "I2", "H2")
    
    result_meta$p.adjust.fdr <- stats::p.adjust(result_meta[,4], method = "fdr")
    result_meta$p.adjust.BY  <- stats::p.adjust(result_meta[,4], method = "BY")
    result_meta$logFDR <- -log10(result_meta$p.adjust.fdr)
    
    # Significant genes
    print(paste("Genes significativos ( p<",args$pvalue ,"):",sum(result_meta[, "pvalue"] < args$pvalue )))
    print(paste("Genes significativos ( FDR<",args$pvalue ,"):",sum(result_meta[, "p.adjust.fdr"] < args$pvalue )))
    
    # Add number of studies where the gene is evaluated
    n.studies <-  ncol(mat.logFC) - sum.NA
    table(n.studies)
    result_meta$n.studies <- n.studies[rownames(mat.logFC)]
    
    # Annotation (biomart annotation file Gcrh38)
    annot <- read.delim("/clinicfs/projects/i63/tfm_hipocampo/MAST/transcript_length_borja.txt", header = TRUE, sep = "\t")
    annot <- annot[!duplicated(annot$Gene.stable.ID),]
    rownames(annot) <- annot$Gene.stable.ID
    
    idx <- match(rownames(result_meta),rownames(annot))
    result_meta$SYMBOL <- annot[idx,"HGNC.symbol"]
    result_meta <- relocate(result_meta,"SYMBOL",.before = "lower_bound")
    
    sig.genes.df <- result_meta[result_meta$p.adjust.fdr < args$pvalue ,]
    
    write.table(x = sig.genes.df[order(sig.genes.df$p.adjust.fdr),] %>%
                  rownames_to_column('ENSEMBL'), file = "sig.genes.tsv",
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(x = result_meta[order(result_meta$p.adjust.fdr),] %>%
                  rownames_to_column('ENSEMBL'), file = "all.genes.tsv",
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    
    # STEP 3. INFLUENCE AND SENSITIVITY ANALYSIS ####
    # ===============================================================
    
    #Add 4 new variables about influence & sensitivity analysis:
    for (i in rownames(sig.genes.df)){
      #print(i)
      # Define studies for each function (not NA)
      estudios <- colnames(mat.logFC)[!mat.SE.NA[i,]]
      
      # Influence info 1:
      # Number of studies where the sign of the logOR is the same  of the global logOR:
      sig.genes.df[i, "infl.same.sign.logFC"] <- sum(sign(MA[[i]]$yi)== rep(sign(MA[[i]]$b),length(estudios)))
      
      # Influence info 2: how many studies could be influencers?
      inf <- influence(MA[[i]])
      res <- paste(estudios[inf$is.infl], collapse = ",")
      sig.genes.df[i, "infl.nstudies"] <- ifelse(res =="", "non", res)
      
      # Sensivity analysis
      l1 <-as.data.frame(leave1out(MA[[i]]))
      rownames(l1) <- estudios
      
      #1. p.value about differences between all estimates from leave one out
      #   and global estimate)
      sig.genes.df[i, "sensi.global"] <-t.test(x= l1$estimate,
                                               mu=as.numeric(MA[[i]]$b))$p.value
      #2. number of  studies where pvalue > 0.05
      # (we hope p-values < 0.05, significant estimates)
      res2 <- paste(estudios[l1$pval > 0.05], collapse = ",")
      sig.genes.df[i, "sensi.specific"] <- ifelse(res2 =="", "all.p.values < 0.05", res2)
    }
    
    ## QUESTIONS TO ASSESS META-ANALYSIS FOR EACH FUNCTION:
    
    #1. INFLUENCE STUDIES. How many logOR have the same sign to global logOR?
    cat("\n\n LORs individuales con el mismo signo que el LOR global \n")
    print(table(sig.genes.df$infl.same.sign.logFC))
    
    #2. INFLUENCE STUDIES. How many functions including influence studies?
    cat("\n\n Genes con estudios de influencia \n")
    print(table(sig.genes.df$infl.nstudies=="non"))
    
    #3. SENSITIVITY. In global, are there many functions with differences in the estimate?
    cat("\n\n Resultados leave 1 out (sensi.global < 0.05) \n")
    print(table(sig.genes.df$sensi.global < 0.05))
    
    #4. SENSITIVITY.  How many functions including changes in the significance about
    # its new estimate  after leave1out?
    cat("\n\n leave 1 out no cambia el p valor \n")
    print(table(sig.genes.df$sensi.specific == "all.p.values < 0.05"))
    
    
    # Save final results:
    write.table(x = sig.genes.df[order(sig.genes.df$p.adjust.fdr),] %>% rownames_to_column('ENSEMBL'), 
                file = "sig.genes_def.tsv", sep ="\t", quote = F, row.names = F)
    
    
    # STEP 4. Visualization of significant genes ####
    # ===============================================================
    cat("\n\nGeneración de informes a partir de los genes significativos.\n")
    
    if (nrow(sig.genes.df) > 0) {
      sig.genes.df.5 <- head(sig.genes.df[order(sig.genes.df$p.adjust.fdr),], n=20)
      #Select significant functions to visualize:
      sig.results <- result_meta[rownames(sig.genes.df.5),]
      
      dir.create("./plots")
      setwd("plots")
      
      selMethod <- "DL"
      
      for (i in 1:nrow(sig.results)){
        mygenes <- rownames(sig.results)[i]
        res <- rma(yi= mat.logFC[mygenes,], sei =mat.SE[mygenes,], method = "DL")
        
        #FOREST PLOT
        # png (filename = paste("FOREST_", mygenes,".png", sep =""), width = 960 ,
        #      height = 960, res = 200)
        png (filename = paste(gsub("-","_",mygenes),"_FOREST",".png", sep =""), width = 960 ,
             height = 960, res = 200)
        forest(res,
               slab = toupper(colnames(mat.logFC)), #Nombre de los estudios
               xlab="logFC", cex=0.7,
               mlab=paste(selMethod, "Model for All Studies", sep = " "),
               border = "black", #Color del borde del rombo
               col = "red", #Color del rombo
               main = paste("\n", mygenes, sep=""))
        text( 9,-3, "logFC [IC 95%]", pos=2, cex = 0.7)
        dev.off()
        
        #FUNNEL PLOT
        png (filename = paste(gsub("-","_",mygenes),"_FUNNEL", ".png", sep =""), width = 960 ,
             height = 960, res = 200)
        par(mfrow=c(2,2))
        funnel(res, main="Standard Error", back ="darkslategray1",
               xlab = paste("logFC (", mygenes, ")",sep =""))
        funnel(res, yaxis="vi", main="Sampling Variance", back ="darkslategray1",
               xlab = paste("logFC (", mygenes, ")",sep =""))
        funnel(res, yaxis="seinv", main="Inverse Standard Error",
               back ="darkslategray1", xlab = paste("logFC (", mygenes, ")",sep =""))
        funnel(res, yaxis="vinv", main="Inverse Sampling Variance",
               back ="darkslategray1",  xlab = paste("logFC (", mygenes, ")",sep =""))
        par(mfrow=c(1,1))
        dev.off()
        
        #INFLUENCE PLOTS
        # That shows various diagnostic measures
        png (filename = paste(gsub("-","_",mygenes), "_INFLUENCE", ".png", sep =""), width = 960 ,
             height = 960, res = 200) 
        inf <- influence(res)
        #plot(inf, plotfb = T)#"plotfb" is not a graphical parameter
        plot(inf)
        dev.off()
      }
    
    # STEP 5. Generating report ####
    # ===============================================================

      # Summary table
      res <- function(df){
        total <- nrow(df)
        sig <- sum(df$p.adjusted < 0.05)
        res <- c(total, sig)
        names(res) <- c("genes", "sig.genes")
        return(res)
      }
      
      tabla_resumen <- function(x) {
        table <- lapply(x, res)
        # Las uno en una tabla
        table <- dplyr::bind_rows(table, .id = "Study")
        
        return(table)
      }
      
      # Summary table with the genes and significant genes (p.adj.val < 0.05)
      table_sum <- tabla_resumen(MAST_results_list)
      # table_sum$Study <- gsub("_ED","", table_sum$Study)
      # if(all(names(apply(!is.na(mat.logFC), 2, sum)) == table_sum$Study)){
      #   table_sum <- dplyr::bind_cols(table_sum, apply(!is.na(mat.logFC), 2, sum))
      #   table_sum <- table_sum %>% dplyr::relocate(...4, .after = genes)
      # }
      
      # Function to create multiple tabs
      make.tabs <- function(sig.genes.df){
        res <- NULL
        for(g in rownames(sig.genes.df)){
          file_name <- gsub("-","_", g)
          res <- c(res, '## ', g, ' {.tabset .tabset-fade .tabset-pills} \n',
                   "**Statistics of ", g, " meta-analisys** \n",
                   "```{r, fig.align='center'}", '\n',
                   "kable(sig.genes.df['",g,"',])", '\n',
                   '```', '\n\n',
                   "[Gene information](https://www.genecards.org/cgi-bin/carddisp.pl?gene=", sig.genes.df[g,"SYMBOL"], ") \n\n",
                   "### Forest plot {-} \n",
                   "```{r, fig.align='center'}", '\n',
                   'knitr::include_graphics("', meta_dir, '/plots/', file_name, '_FOREST.png")\n',
                   '```', ' \n\n',
                   "### Funnel plot {-} \n",
                   "```{r, fig.align='center'}", '\n',
                   'knitr::include_graphics("', meta_dir, '/plots/', file_name, '_FUNNEL.png")\n',
                   '```', ' \n\n',
                   "### Incluence plot {-} \n",
                   "```{r, fig.align='center'}", '\n',
                   'knitr::include_graphics("', meta_dir, '/plots/', file_name, '_INFLUENCE.png")\n',
                   '```', '\n\n')
        }
        return(res)
      }
      
      
      # Generating volcano plotly
      sig_limit <- 0.05
      lfc_limit <- 1.5
      
      data <- result_meta
      data$gene_name <- data$SYMBOL
      rownames(data) <- NULL
      data <- data[,c("gene_name", "logFC", "p.adjust.fdr", "logFDR")]
      
      # Add a grouping column; default value is "not significant"
      data["Group"] <- "NotSignificant"
      
      # Change the grouping for the entries with significance but not a large enough Fold change
      data[which(data['p.adjust.fdr'] < sig_limit & abs(data['logFC']) < lfc_limit ),"Group"] <- "Significant"
      
      # Change the grouping for the entries a large enough Fold change but not a low enough p value
      data[which(data['p.adjust.fdr'] > sig_limit & abs(data['logFC']) > lfc_limit ),"Group"] <- paste0("FoldChange > ",lfc_limit," | < -",lfc_limit)
      
      # Change the grouping for the entries with both significance and large enough fold change
      data[which(data['p.adjust.fdr'] < sig_limit & abs(data['logFC']) > lfc_limit ),"Group"] <- paste0("Significant & FoldChange > ",lfc_limit," | < -",lfc_limit)
      
      data$texto <- paste0("Gene: ", data$gene_name, "\np.val.fdr: ", data$p.adjust.fdr)
      
      volcano_plot <- ggplot(data, aes(x = logFC, y = logFDR)) +
        geom_point(aes(color = Group, text = texto)) +
        scale_color_manual(values = c("#E69F00", "#DB5E00", "#0072B2", "#009E73")) + #https://stackoverflow.com/q/57153428
        xlab("logFC") + ylab("-log10(FDR)")
      
      # Create the Rmd to knit
      ns <- nrow(sig.genes.df)
      cat(
        '---
      title: "Global Meta-analysis"
      output:
        html_document:
          toc: false
          toc_float: false
          code_folding: hide
          number_sections: true
          theme: spacelab
      ---
      ## ', metaanalysis_name, ' {.tabset .tabset-pills -}
      
      ```{r, warning=F, message=F}
      library(dplyr)
      library(knitr)
      library(DT)
      library(ggplot2)
      library(plotly)
      load("', meta_dir, '/', metaanalysis_name, '.RData")
      ```  \n
      ### MA results (significant genes) {-}  \n',
            "```{r, fig.align='center'}", '\n',
            "datatable(sig.genes.df[,c(1:5,12,15)], caption='Statistics of ", ns, " significant genes',rownames=FALSE, escape = FALSE, filter = 'top')", '\n',
            '```', '\n\n',
            '### MA results (volcano plot) {-} \n',
            "```{r, fig.align='center', out.width = '99%'}", '\n',
            "ggplotly(volcano_plot)", '\n',
            '```', '\n\n',
            '### Individual results (summary) {-} \n',
            "```{r, fig.align='center'}", '\n',
            "kable(table_sum, col.names = c('Study','All genes','Significant genes'))", '\n',
            '```',
            '\n\n',
            '### sessionInfo {-}  \n',
            "```{r, fig.align='center'}", '\n',
            "date()", "\n",
            "sessionInfo()", '\n',
            '```', '\n\n',
            sep = "",
            file = paste0(meta_dir, "/", metaanalysis_name, ".Rmd"))
          
          cat(
            '---
      title: "Meta-analysis of genes"
      subtitle: "Draft, ', format(Sys.time(), "%d-%m-%Y"),'v6c"
      output:
        html_document:
          toc: true
          toc_float:
            collapsed: false
            smooth_scroll: false
          code_folding: hide
          number_sections: false
          theme: spacelab
      ---
      # ', metaanalysis_name, ' {- .tabset .tabset-fade .tabset-pills}
      
      ```{r, warning=F, message=F}
      library(dplyr)
      library(knitr)
      library(DT)
      load("', meta_dir, '/', metaanalysis_name, '.RData")
      ```  \n',
            make.tabs(sig.genes.df.5[order(row.names(sig.genes.df.5)),]), "\n\n",
            '## sessionInfo {-}  \n',
            "```{r, fig.align='center'}", '\n',
            "date()", "\n",
            "sessionInfo()", '\n',
            '```', '\n\n',
            sep = "",
            file = paste0(meta_dir, "/", metaanalysis_name, "_genes.Rmd"))
          
          save(sig.genes.df, result_meta, table_sum, volcano_plot, file = paste0(meta_dir, "/", metaanalysis_name, ".RData"))
          # Render the Rmd created into html here ####
          rmarkdown::render(paste0(meta_dir, "/", metaanalysis_name, ".Rmd"))
          rmarkdown::render(paste0(meta_dir, "/", metaanalysis_name, "_genes.Rmd"))
    } else {
      cat("\nNo hay genes significativos. Se omiten los gráficos y no se generarán reportes HTML ni Rmd.\n")
    }
  # Reestablish original dir to continue with another cell type-contrast combination
  setwd("/clinicfs/projects/i63/tfm_hipocampo/")
  }
}

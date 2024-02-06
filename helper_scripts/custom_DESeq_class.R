my_results_transform <- function(result, wb){
  setClass(
    ## custom class for storing results data
    "myDESRclass",
    contains="DESeqResults",
    slots=c("visualizations", "sig.genes", 'gsea')
  ) -> myDESRclass
  ## Convert results object to new class and generate all outputs
  result <- as(result, "myDESRclass")
  result@visualizations$summary <- capture.output(summary(result))
  result@visualizations$volplot <- generate_volplot(result)
  result@visualizations$volplot_fixed <- generate_volplot(result) + ylim(c(0,25)) + xlim(c(-5,5))
  result@visualizations$greg_volplot <- greg_volplot(result)
  # result@visualizations$finalfig <- arrangeGrob(textGrob(paste(result@visualizations$summary[1:8],
  #                                                              collapse="\n")),
  #                       ggarrange(result@visualizations$volplot, result@visualizations$volplot_fixed,
  #                                 common.legend = T),
  #                       result@visualizations$greg_volplot,
  #                       ncol=2, nrow=2, layout_matrix=layout,
  #                       top=textGrob(result@metadata$contrast, gp=gpar(fontsize=25)))#, bottom="padj < 0.1, | FoldChange | > 1.3")#, fig.lab=unlist(out[[i]]$pair))
  result@gsea <- run_fgsea(result, gmt.file, nperm=10000)
  result@gsea$source <- unlist(lapply(result@gsea$pathway, function(x){str_split(x, '_')[[1]][1]}))
  volData <- result[!is.na(result$padj),]
  volData <- volData[order(volData$padj),]
  result@sig.genes <-  volData[volData$padj<analysis$config$alpha,]@rownames
  return(result)
}



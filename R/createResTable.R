createResTable <- function (fit, contr.matrix, ebayes_fun = "eBayes") {
  resultsTables_list <- list()
  
  for (contrast in colnames(fit$coefficients)) {
    if (ebayes_fun == "eBayes") {
      tab <- limma::topTable(fit, coef = contrast, n = Inf) 
    } else if (ebayes_fun == "treat") {
      tab <- limma::topTreat(fit, coef = contrast, n = Inf)
    }
    
    resultsTables_list[[contrast]] <- tab %>% 
      dplyr::rename(log2FoldChange = .data$logFC, 
                    pvalue = .data$P.Value, 
                    padj = .data$adj.P.Val)
    }
    
  resultsTable <- lapply(resultsTables_list, function(one_tbl) {
    one_tbl %>% tibble::rownames_to_column(var = "gene")
  }) %>% 
    dplyr::bind_rows(.id = "contrast") %>% 
    dplyr::as_tibble() %>% 
    dplyr::mutate(contrast = forcats::fct(contrast, levels = colnames(contr.matrix)))
  return(resultsTable)
}

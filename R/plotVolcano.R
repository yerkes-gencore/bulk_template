#' Make Volcano Plots For Differentially Expressed Genes
#'
#' @param bulk bulk list object. Differential expression analysis results must be present. 
#' @param method Differential expression analysis method (deseq or limma). No default.
#' @param contrast Constrast name for visualization. No default.
#' 
#' @returns Volcano Plot with interactive capabilities
#' @export
#'
#' @examples
plotVolcano = function(bulk = NULL, 
                       method = NULL,
                       contrast = NULL){
  ### util
  `%nin%` = Negate(`%in%`)
  
  ### Checkpoints
  if (is.null(bulk)) {
    stop("`bulk` list object must be provided with degTab populated.")
  }
  if (is.null(method) | method %nin% c('limma', 'deseq')) {
    stop("`method` must be either `limma` or `deseq`.")
  }
  if (is.null(contrast)) {
    stop("`contrast` is NULL. Please provide a contrast.")
  }else(con = contrast)
  
  
  ### Set y-axis scale
  y_lim <- range(
    -log10(bulk$degTab[[method]]$pvalue),
    finite = TRUE
  )
  y_lim[2] <- ceiling(y_lim[2])
  
  ### Find maximum unadjusted pvalue where padj is significant
  p_threshold <- bulk$degTab[[method]] %>%
    filter(contrast == con, padj < 0.05) %>%
    summarise(
      threshold = if (n() > 0) max(pvalue, na.rm = TRUE) else NA_real_
    ) %>%
    pull(threshold)
  
  ### Generate Volcano Plot
  plt =  bulk$degTab[[method]] %>%
      filter(contrast == con) %>%
      {
        ggplot(data = .) +
          geom_point(
            aes(
              x = log2FoldChange,
              y = -log(pvalue, base = 10),
              color = signif_dir,
              text = gene
            ),
            size = 0.2
          ) +
          theme_bw() +
          xlim(
            -max(abs(.$log2FoldChange)),
            max(abs(.$log2FoldChange))
          ) +
          ylim(y_lim[1], y_lim[2]) +
          scale_color_manual(
            values = c(
              "NotSig" = "black",
              "Up" = "tomato",
              "Down" = "dodgerblue"
            ),
            drop = FALSE,
            guide = guide_legend(
              override.aes = list(size = 3)
            )
          ) +
          labs(
            color = "Significance (padj < 0.05)",
            title = con,
            x = "Log2 Fold Change",
            y = "p-value"
          ) +
          geom_vline(xintercept = 0, linewidth = 0.5)
      }
  
  ### Add horizontal significance line
  if(!is.na(p_threshold)){
    plt = plt +
      geom_hline(
      yintercept = -log10(p_threshold),
      linetype = "dashed",
      color = 'red'
    )
  }
  
  ### Return Volcano Plot
  return(plt)
    
}

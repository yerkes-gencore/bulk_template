#' Plot Heatmaps
#'
#' @param bulk bulk list object. Differential expression analysis results must be present. 
#' @param genes Vector of genes to include in the heatmap. No Default.
#' @param method Differential expression analysis method (deseq or limma). No default.
#' @param contrast Constrast name for visualization. No default.
#' @param pairing_bar Boolean indicating whether to include a color bar depicting which samples are paired and unpaired across time points. Default is TRUE.
#' @param comparison_borders 
#' @param fdr False discovery rate threshold. Default is 0.05.
#' @param gene_order Method used for ordering genes in the heatmap. Default is `rank` where genes are ordered based on avglog2FC. Alternatives are 'cluster' enabling hierarchical clustering and `none` preventing any further gene ordering.
#' @param ... Pass other params for getTopDegTbl. (ex: filter_pattern = '^ENSG')
#``
#' @returns
#' @export
#'
#' @examples
plotHeatmap = function(bulk = NULL,
                       genes = NULL,
                       method = NULL,
                       contrast = NULL,
                       pairing_bar = TRUE,
                       comparison_borders = TRUE,
                       fdr = 0.05,
                       gene_order = 'rank',
                       ...){
  ### util
  `%nin%` = Negate(`%in%`)
  
  ### Checkpoints
  if (is.null(bulk)) {
    stop("`bulk` list object must be provided with degTab populated.")
  }
  if (is.null(genes)) {
    warning("`genes` vector is NULL. Using DEGs instead.")
  }
  if (is.null(method) | method %nin% c('limma', 'deseq')) {
    stop("`method` must be either `limma` or `deseq`.")
  }
  if (is.null(contrast)) {
    stop("`contrast` is NULL. Please provide a contrast.")
  }else(con = contrast)
  
  
  # defining contrast and pulling out from bulk$md$contr.matrix
  cmSubset <- bulk$md$contr.matrix %>%
    as.data.frame() %>%
    rownames_to_column("group") %>%
    dplyr::select(group, "comparison" = all_of(con))
  
  # Make Genes vector
  if(is.null(genes)){
    
    # subsetting gene table
    gene_tbl <- bulk$degTab[[method]] %>%
      getTopDegTbl(., contrast_i = con, groupID = "group",
                   arrange_by = "pvalue", direction = "unequal",
                   padj_cutoff = NULL, slice_n = 50,
                   ...) %>%
      dplyr::select(gene, padj, l2fc = log2FoldChange) %>%
      dplyr::mutate(color = dplyr::case_when(
        is.na(padj) ~ "grey",
        padj >= fdr ~ "black",
        l2fc > 0 ~ "tomato",
        l2fc < 0 ~ "dodgerblue",
        TRUE ~ NA_character_))
    
    genes = gene_tbl$gene
    
  }else{
    
    # subsetting gene table
    gene_tbl <- bulk$degTab[[method]] %>%
      getTopDegTbl(., contrast_i = con, groupID = "group",
                   arrange_by = "pvalue", direction = "unequal",
                   padj_cutoff = NULL, slice_n = Inf, ...) %>%
      filter(gene %in% genes) %>%
      dplyr::select(gene, padj, l2fc = log2FoldChange) %>%
      dplyr::mutate(color = dplyr::case_when(
        is.na(padj) ~ "grey",
        padj >= fdr ~ "black",
        l2fc > 0 ~ "tomato",
        l2fc < 0 ~ "dodgerblue",
        TRUE ~ NA_character_))
  }
  

  
  # Main plot data frame 
  plot_df <- bulk$tables[[method]]$scaleGenes[c(genes),] %>%
    rownames_to_column("gene") %>%
    tidyr::pivot_longer(
      !gene,
      names_to = "sampleID",
      values_to = "zscore"
    ) %>%
    dplyr::inner_join(
      bulk$dge$samples,
      by = "sampleID"
    ) %>%
    dplyr::inner_join(
      gene_tbl,
      by = "gene"
    ) %>%
    dplyr::inner_join(
      cmSubset,
      by = "group"
    ) %>%
    dplyr::group_by(subjectID) %>%
    dplyr::mutate(
      pairing = if_else(
        n_distinct(group) > 1,
        "Paired",
        "Unpaired"
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(group, pairing, subjectID) %>%
    dplyr::mutate(
      group = factor(group, levels = unique(group)),
      subjectID = factor(subjectID, levels = unique(subjectID)),
      pairing = factor(pairing, levels = c("Paired", "Unpaired"))
    )
  
  
  # getting gene order for heatmap
  if(gene_order == 'cluster'){
    ord <- hclust(dist(bulk$tables[[method]]$scaleGenes[c(genes),]))$order
  }else if(gene_order == 'rank'){
    ord = match(gene_tbl %>% arrange(desc(l2fc)) %>% .$gene,
                rownames(bulk$tables[[method]]$scaleGenes[c(genes),]))
  }else if(gene_order == 'none'){
    ord = seq(1:nrow(bulk$tables[[method]]$scaleGenes[c(genes),]))
  }
  
  # order plot
  plot_df = plot_df %>% 
    mutate(gene = factor(gene, levels = unique(gene)[ord]))
  
  
  # ---------------------------------------
  # Initialize plot
  # ---------------------------------------
  plt =ggplot(
    data = plot_df,
    aes(x = subjectID, y = gene, fill = zscore)
  ) +
    
    # ---------------------------------------
  # Main heatmap
  # ---------------------------------------
  geom_tile() +
    
    scale_fill_gradient2(
      low = "dodgerblue3",
      mid = "white",
      high = "tomato2",
      limits = c(
        min(plot_df$zscore),
        max(plot_df$zscore)
      )
    ) +
    scale_y_discrete(
      limits = c(
        if (pairing_bar) "Pairing",
        rev(levels(plot_df$gene))
      )
    )
  
  # ---------------------------------------
  # Pairing annotation
  # ---------------------------------------
  if(pairing_bar){
    
    plt = plot_df %>%  
      dplyr::distinct(group, subjectID, pairing) %>%
      dplyr::mutate(
        annotation = "Pairing"
      ) %>%
      {
        plt+ggnewscale::new_scale_fill()+
          geom_tile(
            data = .,
            aes(
              x = subjectID,
              y = annotation,
              fill = pairing
            ),
            inherit.aes = FALSE,
            height = 0.8
          ) +
          scale_fill_manual(
            name = "Pairing",
            values = c(
              "Paired" = "black",
              "Unpaired" = "grey80"
            ),
            guide = guide_legend(
              order = 2
            )
          ) 
      }
    
  }
  
  # ---------------------------------------
  # Comparison borders
  # ---------------------------------------
  if(comparison_borders){
    
    plt = plt+geom_rect(
      data = subset(
        plot_df,
        comparison %in% c(0.5, 1, -1, -0.5)
      ),
      aes(colour = factor(comparison)),
      linewidth = 1.5,
      fill = NA,
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf,
      inherit.aes = FALSE
    ) +
      scale_colour_manual(
        values = c(
          "-1"   = "dodgerblue3",
          "-0.5" = "dodgerblue3",
          "0.5"  = "tomato2",
          "1"    = "tomato2"
        ),
        breaks = c("-1", "-0.5", "0.5", "1"),
        labels = c(
          "-1"   = "Reference",
          "-0.5" = "Reference",
          "0.5"  = "Test",
          "1"    = "Test"
        ),
        name = "Comparison direction"
      )
  }
  
  
  
  # ---------------------------------------
  # Themes
  # ---------------------------------------
  plt = plt+facet_grid(
    ~ group,
    labeller = labeller(group = newLineLabel),
    scales = "free_x",
    space = "free_x"
  ) +
    labs(
      title = paste0(con),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold"
      ),
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        color = "black"
      ),
      axis.text.y = element_text(
        color = c(
          setNames(rev(plot_df$color[ord]), rev(levels(plot_df$gene))),
          "Pairing" = "black"
        )[c(if (pairing_bar) "Pairing", rev(levels(plot_df$gene)))],
        size = 10
      ),
      strip.text = element_text(
        size = 10,
        face = "bold"
      )
    )
  
  ### Return Plot
  return(plt)

}

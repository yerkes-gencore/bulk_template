## Custom helper functions for voomEdgeR workflow
## Will move this to gencoreBulk package

## QC

getVoomByGroupData <- function(bulk, group, ...) {
  vbg <- voomByGroup(counts = bulk$dge$counts, 
                     design = bulk$md$design, 
                     group = bulk$dge$samples[[group]], 
                     plot = "combine", save.plot = TRUE)
  
  vbg_plot_data <- 
    lapply(names(vbg$voom.line), function(grp) {
      as_tibble(vbg$voom.line[[grp]]) %>%
        mutate(group = grp) %>%
        distinct()
    }) %>% bind_rows()
  
  return(vbg_plot_data)
}

plotVoomByGroupData <- function(vbg_data, ...) {
  vbg_data %>%
    ggplot(data = ., aes(x = x, y = y, color = group)) +
    geom_line() +
    scale_color_brewer(palette = "Paired") +
    xlab("log2( count size + 0.5 )") +
    ylab("Sqrt( standard deviation )") +
    theme_classic()
}

ggplotMDS <- function(dge, group = "group", sampleID = "sampleID", 
                      gene.selection = "common", 
                      color, path, ...) {
  # Get mds data from edgeR::plotMDS()
  mds_data <- plotMDS(dge, top = 500,
                      labels = dge$samples[[group]], 
                      col = col.group, plot = FALSE, 
                      gene.selection = gene.selection, ...)
  
  mds_xy <- mds_data[c("x","y")] %>% as_tibble()
  mds_xy[[sampleID]] <- colnames(dge)
  # mds_xy[[sampleID]] <- str_extract(colnames(dge), "(^[A-Za-z0-9]*)")
  mds_xy <- mds_xy %>% full_join(dge$samples, by = sampleID)
  
  x_varex <- round(mds_data$var.explained[1]*100, digits = 0)
  y_varex <- round(mds_data$var.explained[2]*100, digits = 0)
  
  mds_xy %>%
    ggplot(aes(x = .data$x, y = .data$y, 
               color = .data[[color]], 
               label = .data[[sampleID]])) +
    geom_text() + 
    geom_path(aes(linetype = .data[[path]])) +
    xlab(paste0(mds_data$axislabel, " 1 (", x_varex, "%)")) +
    ylab(paste0(mds_data$axislabel, " 2 (", y_varex, "%)")) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    coord_cartesian(clip = "off") +
    ggtitle(ifelse(gene.selection == "common", "PCA", "MDS"))
}

## DE

runEdgeRVoomLmFit <- function(bulkObj, contr.matrix, block = NULL, sample.weights = TRUE, var.design = NULL, var.group = NULL, plotVoom = TRUE) {
  ## This function adapts the limma voom method (Law et al, 2014) to allow for loss of residual degrees of freedom due to exact zero counts (Lun and Smyth, 2017). 
  ## The function is analogous to calling voom followed by duplicateCorrelation and lmFit except for the modified residual df values and residual standard deviation sigma values.
  ## If block is specified, then the intra-block correlation is estimated using duplicateCorrelation. In that case, the voom weights and the intra-block correlation are each estimated twice to achieve effective convergence.
  ## Empirical sample quality weights will be estimated if sample.weights=TRUE or if var.design or var.group are non-NULL (Liu et al 2015). In that case, voomLmFit is analogous to running voomWithQualityWeights followed by lmFit.
  bulkObj$fit <- edgeR::voomLmFit(counts = bulkObj$dge$counts, 
                                  design = bulkObj$md$design,
                                  block = block,
                                  sample.weights = sample.weights,
                                  var.design = var.design,
                                  var.group = var.group,
                                  plot = plotVoom)
  ## Re-orient the fitted model object from the coefficients of the original design matrix to the set of contrasts defined above in contr.matrix
  bulkObj$fit <- contrasts.fit(bulkObj$fit, 
                               contrasts = contr.matrix)
  ## Run empirical Bayes moderation; borrowing information across all the genes to obtain more precise estimates of gene-wise variability
  bulkObj$fit <- eBayes(bulkObj$fit, robust=TRUE)
  ## Plot the model's residual variances against average expression values; demonstrating that the variance is no longer dependednt on the mean expression level
  plotSA(bulkObj$fit, main="Final model: Mean-variance trend")
  ## Identify which genes are significantly differentially expressed for each contrast from a fit object containing p-values and test statistics
  bulkObj$de <- decideTests(bulkObj$fit, lfc = 0)
  return(bulkObj)
}

runGlmQLFit <- function(bulkObj, contr.matrix) {
  bulkObj$dge <- estimateDisp(bulkObj$dge, bulkObj$md$design)
  plotBCV(bulkObj$dge)
  bulkObj$fit <- glmQLFit(bulkObj$dge, bulkObj$md$design, robust = TRUE)
  bulkObj$res <- list()
  for (i in 1:ncol(contr.matrix)){
    bulkObj$res[[i]]<- glmQLFTest(bulkObj$fit, contrast = contr.matrix[,i])
  }
  names(bulkObj$res) <- colnames(contr.matrix)
  return(bulkObj)
}

createResTable <- function(clustObj, contr.matrix) {
  resultsTables_list <- list()
  for (contrast in colnames(clustObj$fit$coefficients)) {
    resultsTables_list[[contrast]] <- topTable(clustObj$fit, coef = contrast, n = Inf) %>%
      # as_tibble(rownames = "gene") %>%
      dplyr::rename(log2FoldChange = logFC, pvalue = P.Value, padj = adj.P.Val)
  }
  
  resultsTable <- lapply(resultsTables_list, function(one_tbl) {
    one_tbl %>% rownames_to_column(var = "gene")
  }) %>% bind_rows(., .id = "contrast") %>%
    as_tibble() %>%
    mutate(contrast = fct(contrast, levels = colnames(contr.matrix)))
  
  return(resultsTable)
}

## Check sample sizes per contrast (probably not necessary for bulk studies)

makeSamplesPerContrastTable <- function(contr.matrix, bulkObj) {
  ## Make table of the numerators and denominators of each contrast
  contr.levels <- tibble(contrast = colnames(contr.matrix))
  for (icol in 1:ncol(contr.matrix)) {
    contr.levels[icol, 2] <- names(which(contr.matrix[,icol]==1))
    contr.levels[icol, 3] <- names(which(contr.matrix[,icol]==-1))
  }
  colnames(contr.levels)[2] <- "numerator"
  colnames(contr.levels)[3] <- "denominator"
  
  # Column of sample metadata table to pull grouping factor from
  fct_name <- names(attr(bulkObj$md$design, "contrasts"))
  
  contr.levels_tbl <- contr.levels %>%
    pivot_longer(cols = -contrast, names_to = "type", values_to = "level") %>%
    mutate(level = str_remove(level, paste0("^", fct_name))) %>%
    full_join(count(bulkObj$dge$samples, .data[[fct_name]], name = "nSamples"), 
              by = c("level" = fct_name)) %>%
    select(-level) %>%
    pivot_wider(names_from = type, values_from = nSamples)
  
  return(contr.levels_tbl)
}

## Heatmaps

# Modified from gencoreBulk
plotHeatmapFromGenelist <- function(
    data = NULL, geneList = ., 
    column_split = NULL, 
    column_names_gp_col = NULL,
    row_names_gp_col = NULL,
    slice_labels = NULL, 
    colors = c("blue", "white", "red"), 
    column_labels = colnames(data), 
    slice_labels_rot = 90, 
    slice_labels_col = "black", 
    box_width = unit(3.5, "mm"), 
    box_height = unit(3.5, "mm"), 
    width_buffer = unit(5, "mm"), 
    height_buffer = unit(10, "mm"), 
    column_title = " ", 
    cluster_rows = FALSE, 
    cluster_columns = FALSE, 
    column_gap = unit(2, "mm"), 
    scale_min = -2, 
    scale_max = 2, 
    heatmap_legend_param = list(at = c(scale_min, 0, scale_max), 
                                labels = c(scale_min, 0, scale_max), 
                                title = "log2 fold\ndifference\nfrom\nmedian\nexpression")) {
  
  duds <- geneList[!geneList %in% rownames(data)]
  if (length(duds) > 0) {
    geneList <- geneList[geneList %in% rownames(data)]
    if (length(geneList) == 0) {
      stop("No data for requested genes")
    }
    else {
      message(paste0("Genes ", paste0(duds, collapse = ", "), 
                     " not found in data"))
    }
  }
  
  hmap <- data[geneList, ]
  baseline <- matrixStats::rowMedians(hmap)
  hmap <- hmap - baseline
  ComplexHeatmap::Heatmap(
    hmap, 
    heatmap_legend_param = heatmap_legend_param, 
    width = ncol(hmap) * box_width + width_buffer, 
    height = nrow(hmap) * box_height + height_buffer, 
    column_title = column_title, 
    cluster_rows = cluster_rows, 
    cluster_columns = cluster_columns, 
    column_split = column_split, 
    column_names_gp = grid::gpar(col = column_names_gp_col),
    row_names_gp = grid::gpar(col = row_names_gp_col),
    top_annotation = (if (!is.null(slice_labels)) {
      if (is.null(column_split)) {
        warning("Setting labels requires slices to also be set")
      }
      ComplexHeatmap::HeatmapAnnotation(
        foo = ComplexHeatmap::anno_block(gp = grid::gpar(col = NA), 
                                         labels = slice_labels, 
                                         labels_gp = grid::gpar(col = slice_labels_col, fontsize = 10), 
                                         labels_rot = slice_labels_rot, 
                                         height = unit(2, "cm")))
    }
    else {
      NULL
    }),
    column_gap = column_gap, 
    col = circlize::colorRamp2(c(scale_min, 0, scale_max), colors))
}

plot_heatmap_from_resTable <- function(
    bulk, resultsTable, n_genes = 50, contrast_id, 
    sampleID = "sampleID", groupID = "group", 
    design_md_colors = NULL) {
  hm_tbl <- resultsTable %>%
    filter(contrast == contrast_id) %>%
    filter(!grepl("^ENSMMUG", gene)) %>%
    slice_head(n = n_genes) %>%
    mutate(row_name_color = ifelse(padj >= 0.05, "black",
                                   ifelse(log2FoldChange > 0, "red", 
                                          "blue")))
  
  
  if (is.null(design_md_colors)) {
    design_md_colors <- bulk$dge$samples %>%
      dplyr::select({{ sampleID }}, {{ groupID }}) %>%
      mutate(color = "black")
  }
  # # for label colors
  # design_md_colors_tmp <- design_md_colors %>% 
  #   mutate({{ groupID }} := as.character(.data[[groupID]]),
  #          {{ sampleID }} := as.character(.data[[sampleID]]))
  
  if (is.null(bulk$dge$lcpm)) {
    bulk$dge$lcpm <- cpm(bulk$dge, log = TRUE)
  }
  
  gex_data <- bulk$dge$lcpm[, design_md_colors_tmp[[sampleID]]] 
  gex_data[hm_tbl$gene,]
  
  hm_tbl$gene %>%
    plotHeatmapFromGenelist(geneList = .,
                            data = gex_data,
                            cluster_rows = TRUE,
                            slice_labels = unique(design_md_colors_tmp[[groupID]]),
                            slice_labels_col = unique(design_md_colors_tmp$color),
                            column_split = design_md_colors_tmp[[groupID]],
                            column_names_gp_col = design_md_colors_tmp$color,
                            row_names_gp_col = hm_tbl$row_name_color
    )
}
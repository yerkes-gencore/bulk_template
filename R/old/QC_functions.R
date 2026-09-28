## Part 1:
ggplotMDS <- function(dge, sampleID = "sampleID", gene.selection = "common", 
                      dims = c(1,2),
                      color = NULL, shape = NULL, size = 4,
                      ellipse = NULL, path = NULL,
                      show.labels = TRUE, label.size = 4, custom.labels = NULL, 
                      alpha = 1, ...) {
  # Get mds data from edgeR::plotMDS()
  mds_data <- limma::plotMDS(dge, top = 500, plot = FALSE, 
                             gene.selection = gene.selection, dim.plot = dims, ...)
  
  mds_xy <- mds_data[c("x","y")]
  mds_xy[[sampleID]] <- colnames(dge)
  # mds_xy[[sampleID]] <- str_extract(colnames(dge), "(^[A-Za-z0-9]*)")
  mds_xy <- dplyr::as_tibble(mds_xy) %>% dplyr::full_join(dge$samples, by = sampleID)
  
  x_varex <- round(mds_data$var.explained[dims[1]]*100, digits = 0)
  y_varex <- round(mds_data$var.explained[dims[2]]*100, digits = 0)
  
  if (!is.null(custom.labels)) {
    mds_xy <- mds_xy %>%
      mutate(sampleID = ifelse(sampleID %in% custom.labels, sampleID, NA))
  }
  
  mds_xy %>%
    ggplot(aes(x = .data$x, y = .data$y)) +
    geom_point(
      aes(
        color = (
          if (!is.null(color)) { .data[[color]] } else { NULL }
        ),
        shape = (
          if (!is.null(shape)) { .data[[shape]] } else { NULL }
        ),
        # text = (if (is.null(analysis$qc_config$pcaMapping$hover)) {
        #   NULL
        # } else {
        #   .data[[analysis$qc_config$pcaMapping$hover]]
        # })
      ),
      size = size,
      alpha = alpha
    ) +
    (if (show.labels) { ggrepel::geom_text_repel(aes(label = .data[["sampleID"]]), box.padding = 0.5, na.rm = TRUE) } else { NULL }) +
    (if (!is.null(path)) { geom_path(aes(linetype = .data[[path]]), alpha = 0.5) } else { NULL }) +
    (if (!is.null(ellipse)) { stat_ellipse(aes(color = .data[[ellipse]]), type = "norm", level = 0.67)} else { NULL }) +
    labs(color = color, shape = shape) +
    xlab(paste0(mds_data$axislabel, " ", dims[1], " (", x_varex, "%)")) +
    ylab(paste0(mds_data$axislabel, " ", dims[2]," (", y_varex, "%)")) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    coord_cartesian(clip = "off") +
    ggtitle(ifelse(gene.selection == "common", "PCA", "MDS"))
}

## Part 2:
fitVoomLm <- function(bulkObj, design = bulk$md$design, contr.matrix, 
                      block = NULL, sample.weights = TRUE, 
                      var.design = NULL, var.group = NULL, plotVoom = TRUE,
                      ebayes_fun = "eBayes", fc = 1.2) {
  bulkObj$fit <- edgeR::voomLmFit(counts = bulkObj$dge, # Defaults to the normalized (effective) library sizes in counts if counts is a DGEList or to the columnwise count totals if counts is a matrix.
                                  design = design,
                                  block = block,
                                  sample.weights = sample.weights,
                                  var.design = var.design,
                                  var.group = var.group,
                                  plot = plotVoom)
  bulkObj$fit.contr <- limma::contrasts.fit(bulkObj$fit, contrasts = contr.matrix)
  if (ebayes_fun == "eBayes") {
    bulkObj$fit.contr <- limma::eBayes(bulkObj$fit.contr, robust=TRUE)
  } else if (ebayes_fun == "treat") {
    bulkObj$fit.contr <- limma::treat(bulkObj$fit.contr, fc = fc, robust=TRUE)
  }
  limma::plotSA(bulkObj$fit.contr, main="Final model: Mean-variance trend")
  return(bulkObj)
}

createResTable <- function (fit, contr.matrix, ebayes_fun = "eBayes") {
  resultsTables_list <- list()
  for (contrast in colnames(fit$coefficients)) {
    if (ebayes_fun == "eBayes") {
      resultsTables_list[[contrast]] <- limma::topTable(fit, 
                                                        coef = contrast, n = Inf) %>% dplyr::rename(log2FoldChange = .data$logFC, 
                                                                                                    pvalue = .data$P.Value, padj = .data$adj.P.Val)
    } else if (ebayes_fun == "treat") {
      resultsTables_list[[contrast]] <- limma::topTreat(fit, 
                                                        coef = contrast, n = Inf) %>% dplyr::rename(log2FoldChange = .data$logFC, 
                                                                                                    pvalue = .data$P.Value, padj = .data$adj.P.Val)
    }
  }
  resultsTable <- lapply(resultsTables_list, function(one_tbl) {
    one_tbl %>% tibble::rownames_to_column(var = "gene")
  }) %>% dplyr::bind_rows(.id = "contrast") %>% dplyr::as_tibble() %>% 
    dplyr::mutate(contrast = forcats::fct(contrast, levels = colnames(contr.matrix)))
  return(resultsTable)
}

checkRLEv2 <- function(bulk, DEmethod = "DESeq2") {
  raw_counts <- bulk$dge$counts
  rawLogCounts <- log(raw_counts[rowMins(raw_counts) > 0, ])
  rawMedianLogs <- matrixStats::rowMedians(rawLogCounts)
  rawLogRatios <- rawLogCounts - rawMedianLogs
  
  if(DEmethod == "DESeq2"){
    normCounts <- bulk$deseq$normCounts
    normTitle <- "DESeq2 RLE Normalized"
  } else if(DEmethod == "limma") { 
    normCounts <- bulk$dge$cpm
    normTitle <- "Limma RLE Normalized"
  }
  
  normLogCounts <- log(normCounts)
  
  normMedianLogs <- matrixStats::rowMedians(normLogCounts)
  normLogRatios <- normLogCounts - normMedianLogs
  RLE_raw <- .plotRLE(rawLogRatios, "RLE Raw")
  RLE_norm <- .plotRLE(normLogRatios, normTitle)
  return(list(RLE_raw = RLE_raw, RLE_norm = RLE_norm))
}
.plotRLE <- function(data, title) {
  tibble::as_tibble(data, rownames = NA) %>%
    tidyr::pivot_longer(everything(), names_to = "Sample", values_to = "RLE") %>%
    ggplot(aes(x = .data$Sample, y = .data$RLE)) +
    geom_hline(yintercept = 0, color = "red") +
    geom_violin(draw_quantiles = c(0.25, 0.75), trim = TRUE, color = "lightgreen", alpha = 0.1) +
    geom_boxplot(alpha = 0) +
    theme_bw() +
    theme(axis.text.x = element_text(vjust = 0.5, angle = 90), axis.title.x = element_blank(), aspect.ratio = 0.55) +
    ggtitle(title) +
    aes(forcats::fct_inorder(.data$Sample)) #+ scale_y_continuous(limits = c(-9,3.25), expand = c(0,0))
}

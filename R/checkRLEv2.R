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

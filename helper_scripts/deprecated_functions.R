plotCountsGKT <- function (dds, gene, intgroup = "condition", normalized = TRUE, 
                           transform = TRUE, main, xlab = "group", returnData = FALSE, 
                           replaced = FALSE, pc, color, ...) 
{
  stopifnot(length(gene) == 1 & (is.character(gene) | (is.numeric(gene) & 
                                                         (gene >= 1 & gene <= nrow(dds)))))
  if (!all(intgroup %in% names(colData(dds)))) 
    stop("all variables in 'intgroup' must be columns of colData")
  if (!returnData) {
    if (!all(sapply(intgroup, function(v) is(colData(dds)[[v]], 
                                             "factor")))) {
      stop("all variables in 'intgroup' should be factors, or choose returnData=TRUE and plot manually")
    }
  }
  if (missing(pc)) {
    pc <- if (transform) 
      0.5
    else 0
  }
  if (is.null(sizeFactors(dds)) & is.null(normalizationFactors(dds))) {
    dds <- estimateSizeFactors(dds)
  }
  cnts <- counts(dds, normalized = normalized, replaced = replaced)[gene, 
  ]
  group <- if (length(intgroup) == 1) {
    colData(dds)[[intgroup]]
  }
  else if (length(intgroup) == 2) {
    lvls <- as.vector(t(outer(levels(colData(dds)[[intgroup[1]]]), 
                              levels(colData(dds)[[intgroup[2]]]), function(x, 
                                                                            y) paste(x, y, sep = ":"))))
    droplevels(factor(apply(as.data.frame(colData(dds)[, 
                                                       intgroup, drop = FALSE]), 1, paste, collapse = ":"), 
                      levels = lvls))
  }
  else {
    factor(apply(as.data.frame(colData(dds)[, intgroup, drop = FALSE]), 
                 1, paste, collapse = ":"))
  }
  #  data <- data.frame(count = cnts + pc, group = as.integer(group)) # orig
  #  dataGKT <- data.frame(count = cnts, group = as.integer(group))
  #  dataGKT <- data.frame(count = counts(dds, normalized = normalized, replaced = replaced)[gene,], colData(dds)[intgroup])
  #  dataGKT <- data.frame(count = counts(dds, normalized = normalized, replaced = replaced)[gene,], group = as.factor(group))
  dataGKT <- data.frame(count = counts(dds, normalized = normalized, replaced = replaced)[gene,], group = group)
  #  str(dataGKT) # diagnostic
  #  str(group) # diagnostic
  #  str(xgroup) #diagnostic
  logxy <- if (transform) 
    "y"
  else ""
  if (missing(main)) {
    main <- if (is.numeric(gene)) {
      rownames(dds)[gene]
    }
    else {
      gene
    }
  }
  ylab <- ifelse(normalized, "normalized count", "count")
  if (returnData) 
    #    return(data.frame(count = data$count, colData(dds)[intgroup])) # orig
    return(data.frame(count = dataGKT$count, colData(dds)[intgroup]))
  #    plot(data$group + runif(ncol(dds), -0.05, 0.05), data$count, 
  #       xlim = c(0.5, max(data$group) + 0.5), log = logxy, xaxt = "n", 
  #       xlab = xlab, ylab = ylab, main = main, ...)
  #  axis(1, at = seq_along(levels(group)), levels(group))
  (ggplot(dataGKT, aes(group, count))
    + geom_violin()
    + geom_boxplot(width = 0.2,outlier.shape = NA)
    #    + geom_sina()
    + xlab(paste(intgroup,collapse = ":"))
    + theme(axis.text.x = element_text(angle = 90))
    + labs(title = main)
  )
}


as.mutate.data.frame <- function(object, var = "gene_id", subset = "all", relabel = "", contrastlabel = "", select = c("gene","baseMean","log2FoldChange","pvalue","padj","significant","filterThreshold","subset","contrast"))
{
  df <- rownames_to_column(as.data.frame(object), var = var) %>%
    mutate(significant = factor(if_else(is.na(padj),"removed",if_else(padj<object@metadata$alpha,"significant","not_significant")),
                                levels = c("significant","not_significant","removed")),
           subset = subset,
           contrast = contrastlabel,
           filterThreshold = as.numeric(object@metadata$filterThreshold)) %>%
    rename_at(vars(2:7), function(x) paste0(names(object)[1:6], " (", object@elementMetadata$description[1:6], ")"))
  
  if(subset != "all") {
    subset
    names(df)[2:7]
    #    df <- df %>% rename_at(vars(2:7), function(x) gsub(": ", paste0(": ",subset," "), gsub("all",subset, names(df)[2:7])))
    df <- df %>% rename_at(vars(2:7), function(x) gsub("all", subset, names(df)[2:7]))
  }
  if(relabel != "") {
    df <- df %>% rename_at(vars(2:7), function(x) gsub(relabel, "", names(df)[2:7], fixed = TRUE))
    relabel_nous <- gsub("_"," ",relabel)
    df <- df %>% rename_at(vars(2:7), function(x) gsub(relabel_nous, "", names(df)[2:7], fixed = TRUE))
    # df <- df %>% rename_at(vars(2:7), function(x) gsub("p-values", paste0("p-values: ",newlabel), names(df)[2:7], fixed = TRUE))
  }
  # if(relabel != "" & newlabel != "") {
  #   df <- df %>% rename_at(vars(2:7), function(x) gsub(relabel, newlabel, names(df)[2:7], fixed = TRUE))
  #   relabel_nous <- gsub("_"," ",relabel)
  #   df <- df %>% rename_at(vars(2:7), function(x) gsub(relabel_nous, newlabel, names(df)[2:7], fixed = TRUE))
  #   df <- df %>% rename_at(vars(2:7), function(x) gsub("p-values", paste0("p-values: ",newlabel), names(df)[2:7], fixed = TRUE))
  # }
  if(length(select) > 1) {
    df <- df %>% select(contains(select))
  }
  df
}

greg_volplot <- function(result){
  ## Greg's volcano plot for his filter
  gVolData <- as.mutate.data.frame(result)
  gVolPlot <- gVolData %>%
    dplyr::arrange(!is.na(.data[[colnames(gVolData)[5]]]),desc(.data[[colnames(gVolData)[5]]])) %>%
    ggplot(aes(x=.data[[colnames(gVolData)[2]]],
               y=.data[[colnames(gVolData)[3]]])) +
    geom_vline(aes(xintercept = filterThreshold), color = "goldenrod") +
    geom_text(aes(x = filterThreshold, y=-5,
                  label = format(filterThreshold, digits=3)),
              color = "goldenrod", hjust = 0.5, vjust = 1) +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    geom_point(aes(color = significant), size = 0.1) +
    scale_x_log10() +
    theme_bw() +
    scale_color_manual(values = c("significant" = "red",
                                  "not_significant" = "green",
                                  "removed" = "blue",
                                  "filterThreshold" = "goldenrod"), drop = FALSE) + labs(y="Log2FoldChange") + labs(y="Log2FoldChange")
  return(gVolPlot)
}

##### Derrik functions

### makes a heatmap of the given list of genes, separating samples slightly by group variable





#### DESeq2 results object
### Generate volcano plot from DESeq2 results table
### also references config GeneList to label genes of interest


gen_text_summary <- function(result){
  ## Text summary
  volData <- result[!is.na(result$padj),]
  volData <- volData[order(volData$padj),]
  summary_out <- capture.output(summary(result))
  
}

make_worksheet <- function(result, wb){
  ## Write results to excel sheeet in an open workbook
  pair <- result@metadata$contrast
  if (nchar(pair)>31){
    pair <- substr(pair,1,31)
  }
  volData <- result[!is.na(result$padj),]
  volData <- volData[order(volData$padj),]
  addWorksheet(wb, pair)
  writeData(wb, sheet=pair, x=as.data.frame(volData), rowNames=TRUE)
}

## Deprecated
# write_output_workbook <- function(wb, results, analysis){
#   change_colnames <- function(x){
#     a <- rownames_to_column(as.data.frame(x), var="Gene")
#     colnames(a) <- c("Gene", "baseMean (mean of normalized counts)", "log2FoldChange (Log2 fold change MLE)", "lfcSE (Log fold change standard error)", "stat","pvalue (Wald test p-value)", "padj (BH adjusted p-values)")
#     a
#   }
#   merge_tables <- lapply(results, change_colnames)
#   suffixes <- paste0(".", names(merge_tables))
#   res <- merge_tables[[1]]
#   for (i in head(seq_along(merge_tables),-1)){
#     res <- merge(res, merge_tables[[i+1]], all=TRUE,
#                  suffixes = suffixes[i:(i++1)], by="Gene")
#   }
#   res <- res %>% select(-contains("stat"))
#   addWorksheet(wb, "all")
#   writeData(wb, sheet="all", x=res, rowNames=FALSE)
#   
#   tmp <- analysis$rldDrop
#   baseline <- rowMeans(assay(tmp[,tmp@colData$Group %in% analysis$config$designBaseline]))
#   tmp <- assay(tmp) - baseline
#   addWorksheet(wb, "rle")
#   writeData(wb, sheet="rle", x= tmp, rowNames=TRUE)
# }

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

##### GSEA

### Read Cls File for fGSEA
GSEA.ReadClsFile <- function(file = "NULL") {
  ### Read Cls File for fGSEA
  cls.cont <- readLines(file)
  num.lines <- length(cls.cont)
  cls.cont[[3]] <- gsub("\\t", " ", cls.cont[[3]])  #Converts any tabs to spaces
  class.list <- unlist(strsplit(cls.cont[[3]], " "))  #Splits CLS on spaces
  s <- length(class.list)
  t <- table(class.list)[c(unique(class.list))]
  l <- length(t)
  phen <- vector(length = l, mode = "character")
  phen.label <- vector(length = l, mode = "numeric")
  class.v <- vector(length = s, mode = "numeric")
  for (i in 1:l) {
    phen[i] <- noquote(names(t)[i])
    phen.label[i] <- i - 1
  }
  for (i in 1:s) {
    for (j in 1:l) {
      if (class.list[i] == phen[j]) {
        class.v[i] <- phen.label[j]
      }
    }
  }
  return(list(phen = phen, class.v = class.v))
}


GSEA.GeneRanking <- function(A, class.labels, gene.labels, nperm, permutation.type = 0, 
                             sigma.correction = "GeneCluster", fraction = 1, replace = F, reverse.sign = F, 
                             rank.metric) {
  
  A <- A + 1e-08
  B <- A
  B[is.na(B)] <- 0
  
  N <- length(A[, 1])
  Ns <- length(A[1, ])
  
  subset.mask <- matrix(0, nrow = Ns, ncol = nperm)
  reshuffled.class.labels1 <- matrix(0, nrow = Ns, ncol = nperm)
  reshuffled.class.labels2 <- matrix(0, nrow = Ns, ncol = nperm)
  class.labels1 <- matrix(0, nrow = Ns, ncol = nperm)
  class.labels2 <- matrix(0, nrow = Ns, ncol = nperm)
  
  order.matrix <- matrix(0, nrow = N, ncol = nperm)
  obs.order.matrix <- matrix(0, nrow = N, ncol = nperm)
  rnk.matrix <- matrix(0, nrow = N, ncol = nperm)
  obs.rnk.matrix <- matrix(0, nrow = N, ncol = nperm)
  
  obs.gene.labels <- vector(length = N, mode = "character")
  obs.gene.descs <- vector(length = N, mode = "character")
  obs.gene.symbols <- vector(length = N, mode = "character")
  
  M1 <- matrix(0, nrow = N, ncol = nperm)
  M2 <- matrix(0, nrow = N, ncol = nperm)
  S1 <- matrix(0, nrow = N, ncol = nperm)
  S2 <- matrix(0, nrow = N, ncol = nperm)
  Gn1 <- matrix(0, nrow = N, ncol = nperm)
  Gn2 <- matrix(0, nrow = N, ncol = nperm)
  
  gc()
  
  C <- split(class.labels, class.labels)
  class1.size <- length(C[[1]])
  class2.size <- length(C[[2]])
  class1.index <- seq(1, class1.size, 1)
  class2.index <- seq(class1.size + 1, class1.size + class2.size, 1)
  
  for (r in 1:nperm) {
    class1.subset <- sample(class1.index, size = ceiling(class1.size * fraction), 
                            replace = replace)
    class2.subset <- sample(class2.index, size = ceiling(class2.size * fraction), 
                            replace = replace)
    class1.subset.size <- length(class1.subset)
    class2.subset.size <- length(class2.subset)
    subset.class1 <- rep(0, class1.size)
    for (i in 1:class1.size) {
      if (is.element(class1.index[i], class1.subset)) {
        subset.class1[i] <- 1
      }
    }
    subset.class2 <- rep(0, class2.size)
    for (i in 1:class2.size) {
      if (is.element(class2.index[i], class2.subset)) {
        subset.class2[i] <- 1
      }
    }
    subset.mask[, r] <- as.numeric(c(subset.class1, subset.class2))
    fraction.class1 <- class1.size/Ns
    fraction.class2 <- class2.size/Ns
    
    if (permutation.type == 0) {
      # random (unbalanced) permutation
      full.subset <- c(class1.subset, class2.subset)
      label1.subset <- sample(full.subset, size = Ns * fraction.class1)
      reshuffled.class.labels1[, r] <- rep(0, Ns)
      reshuffled.class.labels2[, r] <- rep(0, Ns)
      class.labels1[, r] <- rep(0, Ns)
      class.labels2[, r] <- rep(0, Ns)
      for (i in 1:Ns) {
        m1 <- sum(!is.na(match(label1.subset, i)))
        m2 <- sum(!is.na(match(full.subset, i)))
        reshuffled.class.labels1[i, r] <- m1
        reshuffled.class.labels2[i, r] <- m2 - m1
        if (i <= class1.size) {
          class.labels1[i, r] <- m2
          class.labels2[i, r] <- 0
        } else {
          class.labels1[i, r] <- 0
          class.labels2[i, r] <- m2
        }
      }
    } else if (permutation.type == 1) {
      # proportional (balanced) permutation
      
      class1.label1.subset <- sample(class1.subset, size = ceiling(class1.subset.size * 
                                                                     fraction.class1))
      class2.label1.subset <- sample(class2.subset, size = floor(class2.subset.size * 
                                                                   fraction.class1))
      reshuffled.class.labels1[, r] <- rep(0, Ns)
      reshuffled.class.labels2[, r] <- rep(0, Ns)
      class.labels1[, r] <- rep(0, Ns)
      class.labels2[, r] <- rep(0, Ns)
      for (i in 1:Ns) {
        if (i <= class1.size) {
          m1 <- sum(!is.na(match(class1.label1.subset, i)))
          m2 <- sum(!is.na(match(class1.subset, i)))
          reshuffled.class.labels1[i, r] <- m1
          reshuffled.class.labels2[i, r] <- m2 - m1
          class.labels1[i, r] <- m2
          class.labels2[i, r] <- 0
        } else {
          m1 <- sum(!is.na(match(class2.label1.subset, i)))
          m2 <- sum(!is.na(match(class2.subset, i)))
          reshuffled.class.labels1[i, r] <- m1
          reshuffled.class.labels2[i, r] <- m2 - m1
          class.labels1[i, r] <- 0
          class.labels2[i, r] <- m2
        }
      }
    }
  }
  
  if (rank.metric == "S2N") {
    
    # compute S2N for the random permutation matrix
    P <- reshuffled.class.labels1 * subset.mask
    for (m in 1:nperm) {
      P2 <- do.call("rbind", replicate(nrow(A), P[,m], simplify = FALSE))
      P2[is.na(A)] <- NA
      Gn1[,m] <- rowSums(P2, na.rm=TRUE)}
    M1 <- B %*% P
    M1 <- M1/Gn1
    gc()
    B2 <- B * B
    S1 <- B2 %*% P
    S1 <- S1/Gn1 - M1 * M1
    S1 <- sqrt(abs((Gn1/(Gn1 - 1)) * S1))
    gc()
    P <- reshuffled.class.labels2 * subset.mask
    for (m in 1:nperm) {
      P2 <- do.call("rbind", replicate(nrow(A), P[,m], simplify = FALSE))
      P2[is.na(A)] <- NA
      Gn2[,m] <- rowSums(P2, na.rm=TRUE)}
    M2 <- B %*% P
    M2 <- M2/Gn2
    gc()
    B2 <- B * B
    S2 <- B2 %*% P
    S2 <- S2/Gn2 - M2 * M2
    S2 <- sqrt(abs((Gn2/(Gn2 - 1)) * S2))
    rm(P)
    rm(B2)
    gc()
    
    if (sigma.correction == "GeneCluster") {
      # small sigma 'fix' as used in GeneCluster
      S2 <- ifelse(0.2 * abs(M2) < S2, S2, 0.2 * abs(M2))
      S2 <- ifelse(S2 == 0, 0.2, S2)
      S1 <- ifelse(0.2 * abs(M1) < S1, S1, 0.2 * abs(M1))
      S1 <- ifelse(S1 == 0, 0.2, S1)
      gc()
    }
    
    M1 <- M1 - M2
    rm(M2)
    gc()
    S1 <- S1 + S2
    rm(S2)
    gc()
    
    rnk.matrix <- M1/S1
    
    if (reverse.sign == T) {
      rnk.matrix <- -rnk.matrix
    }
    gc()
    
    for (r in 1:nperm) {
      order.matrix[, r] <- order(rnk.matrix[, r], decreasing = T)
    }
    
    # compute S2N for the 'observed' permutation matrix  
    P <- class.labels1 * subset.mask
    for (m in 1:nperm) {
      P2 <- do.call("rbind", replicate(nrow(A), P[,m], simplify = FALSE))
      P2[is.na(A)] <- NA
      Gn1[,m] <- rowSums(P2, na.rm=TRUE)}
    M1 <- B %*% P
    M1 <- M1/Gn1
    gc()
    B2 <- B * B
    S1 <- B2 %*% P
    S1 <- S1/Gn1 - M1 * M1
    S1 <- sqrt(abs((Gn1/(Gn1 - 1)) * S1))
    gc()
    P <- class.labels2 * subset.mask
    for (m in 1:nperm) {
      P2 <- do.call("rbind", replicate(nrow(A), P[,m], simplify = FALSE))
      P2[is.na(A)] <- NA
      Gn2[,m] <- rowSums(P2, na.rm=TRUE)}
    M2 <- B %*% P
    M2 <- M2/Gn2
    gc()
    B2 <- B * B
    S2 <- B2 %*% P
    S2 <- S2/Gn2 - M2 * M2
    S2 <- sqrt(abs((Gn2/(Gn2 - 1)) * S2))
    rm(P)
    rm(B2)
    gc()
    
    if (sigma.correction == "GeneCluster") {
      # small sigma 'fix' as used in GeneCluster
      S2 <- ifelse(0.2 * abs(M2) < S2, S2, 0.2 * abs(M2))
      S2 <- ifelse(S2 == 0, 0.2, S2)
      S1 <- ifelse(0.2 * abs(M1) < S1, S1, 0.2 * abs(M1))
      S1 <- ifelse(S1 == 0, 0.2, S1)
      gc()
    }
    
    M1 <- M1 - M2
    rm(M2)
    gc()
    S1 <- S1 + S2
    rm(S2)
    gc()
    
    obs.rnk.matrix <- M1/S1
    gc()
  }
  if (rank.metric == "ttest") {
    
    # compute TTest for the random permutation matrix
    P <- reshuffled.class.labels1 * subset.mask
    for (m in 1:nperm) {
      P2 <- do.call("rbind", replicate(nrow(A), P[,m], simplify = FALSE))
      P2[is.na(A)] <- NA
      Gn1[,m] <- rowSums(P2, na.rm=TRUE)}
    M1 <- B %*% P
    M1 <- M1/Gn1
    gc()
    B2 <- B * B
    S1 <- B2 %*% P
    S1 <- S1/Gn1 - M1 * M1
    S1 <- sqrt(abs((Gn1/(Gn1 - 1)) * S1))
    gc()
    P <- reshuffled.class.labels2 * subset.mask
    for (m in 1:nperm) {
      P2 <- do.call("rbind", replicate(nrow(A), P[,m], simplify = FALSE))
      P2[is.na(A)] <- NA
      Gn2[,m] <- rowSums(P2, na.rm=TRUE)}
    M2 <- B %*% P
    M2 <- M2/Gn2
    gc()
    B2 <- B * B
    S2 <- B2 %*% P
    S2 <- S2/Gn2 - M2 * M2
    S2 <- sqrt(abs((Gn2/(Gn2 - 1)) * S2))
    rm(P)
    rm(B2)
    gc()
    
    if (sigma.correction == "GeneCluster") {
      # small sigma 'fix' as used in GeneCluster
      S2 <- ifelse(0.2 * abs(M2) < S2, S2, 0.2 * abs(M2))
      S2 <- ifelse(S2 == 0, 0.2, S2)
      S1 <- ifelse(0.2 * abs(M1) < S1, S1, 0.2 * abs(M1))
      S1 <- ifelse(S1 == 0, 0.2, S1)
      gc()
    }
    
    M1 <- M1 - M2
    rm(M2)
    gc()
    S1 <- (S1^2)/Gn1
    S2 <- (S2^2)/Gn2
    S1 <- S1 + S2
    S1 <- sqrt(S1)
    rm(S2)
    gc()
    
    rnk.matrix <- M1/S1
    
    if (reverse.sign == T) {
      rnk.matrix <- -rnk.matrix
    }
    gc()
    
    for (r in 1:nperm) {
      order.matrix[, r] <- order(rnk.matrix[, r], decreasing = T)
    }
    
    # compute TTest for the 'observed' permutation matrix  
    P <- class.labels1 * subset.mask
    for (m in 1:nperm) {
      P2 <- do.call("rbind", replicate(nrow(A), P[,m], simplify = FALSE))
      P2[is.na(A)] <- NA
      Gn1[,m] <- rowSums(P2, na.rm=TRUE)}
    M1 <- B %*% P
    M1 <- M1/Gn1
    gc()
    B2 <- B * B
    S1 <- B2 %*% P
    S1 <- S1/Gn1 - M1 * M1
    S1 <- sqrt(abs((Gn1/(Gn1 - 1)) * S1))
    gc()
    P <- class.labels2 * subset.mask
    for (m in 1:nperm) {
      P2 <- do.call("rbind", replicate(nrow(A), P[,m], simplify = FALSE))
      P2[is.na(A)] <- NA
      Gn2[,m] <- rowSums(P2, na.rm=TRUE)}
    M2 <- B %*% P
    M2 <- M2/Gn2
    gc()
    B2 <- B * B
    S2 <- B2 %*% P
    S2 <- S2/Gn2 - M2 * M2
    S2 <- sqrt(abs((Gn2/(Gn2 - 1)) * S2))
    rm(P)
    rm(B2)
    gc()
    
    if (sigma.correction == "GeneCluster") {
      # small sigma 'fix' as used in GeneCluster
      S2 <- ifelse(0.2 * abs(M2) < S2, S2, 0.2 * abs(M2))
      S2 <- ifelse(S2 == 0, 0.2, S2)
      S1 <- ifelse(0.2 * abs(M1) < S1, S1, 0.2 * abs(M1))
      S1 <- ifelse(S1 == 0, 0.2, S1)
      gc()
    }
    
    M1 <- M1 - M2
    rm(M2)
    gc()
    S1 <- (S1^2)/Gn1
    S2 <- (S2^2)/Gn2
    S1 <- S1 + S2
    S1 <- sqrt(S1)
    rm(S2)
    gc()
    
    obs.rnk.matrix <- M1/S1
    gc()
  }
  
  if (reverse.sign == T) {
    obs.rnk.matrix <- -obs.rnk.matrix
  }
  for (r in 1:nperm) {
    obs.order.matrix[, r] <- order(obs.rnk.matrix[, r], decreasing = T)
  }
  
  return(list(rnk.matrix = rnk.matrix, obs.rnk.matrix = obs.rnk.matrix, order.matrix = order.matrix, 
              obs.order.matrix = obs.order.matrix))
}

GSEA.ReadClsFile <-
  function(file = "NULL") {
    
    cls.cont <- readLines(file)
    num.lines <- length(cls.cont)
    cls.cont[[3]] <- gsub("\\t", " ", cls.cont[[3]])  #Converts any tabs to spaces
    class.list <- unlist(strsplit(cls.cont[[3]], " "))  #Splits CLS on spaces
    s <- length(class.list)
    t <- table(class.list)[c(unique(class.list))]
    l <- length(t)
    phen <- vector(length = l, mode = "character")
    phen.label <- vector(length = l, mode = "numeric")
    class.v <- vector(length = s, mode = "numeric")
    for (i in 1:l) {
      phen[i] <- noquote(names(t)[i])
      phen.label[i] <- i - 1
    }
    for (i in 1:s) {
      for (j in 1:l) {
        if (class.list[i] == phen[j]) {
          class.v[i] <- phen.label[j]
        }
      }
    }
    return(list(phen = phen, class.v = class.v))
  }


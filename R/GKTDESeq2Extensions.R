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

plotCountsVsVarGKT <- function (dds, gene, var, intgroup = "condition", normalized = TRUE, 
                                transform = TRUE, main, xlab = "group", returnData = FALSE, 
                                replaced = FALSE, pc, color, ...) 
{
  stopifnot(length(gene) == 1 & (is.character(gene) | (is.numeric(gene) & 
                                                         (gene >= 1 & gene <= nrow(dds)))))
  if (!all(var %in% names(colData(dds))))
    stop("var must be column of colData")
  
  # if (!all(intgroup %in% names(colData(dds)))) 
  #   stop("all variables in 'intgroup' must be columns of colData")
  # if (!returnData) {
  #   if (!all(sapply(intgroup, function(v) is(colData(dds)[[v]], 
  #                                            "factor")))) {
  #     stop("all variables in 'intgroup' should be factors, or choose returnData=TRUE and plot manually")
  #   }
  # }
  # if (missing(pc)) {
  #   pc <- if (transform) 
  #     0.5
  #   else 0
  # }
  if (is.null(sizeFactors(dds)) & is.null(normalizationFactors(dds))) {
    dds <- estimateSizeFactors(dds)
  }
  cnts <- counts(dds, normalized = normalized, replaced = replaced)[gene, 
  ]
  # group <- if (length(intgroup) == 1) {
  #   colData(dds)[[intgroup]]
  # }
  # else if (length(intgroup) == 2) {
  #   lvls <- as.vector(t(outer(levels(colData(dds)[[intgroup[1]]]), 
  #                             levels(colData(dds)[[intgroup[2]]]), function(x, 
  #                                                                           y) paste(x, y, sep = ":"))))
  #   droplevels(factor(apply(as.data.frame(colData(dds)[, 
  #                                                      intgroup, drop = FALSE]), 1, paste, collapse = ":"), 
  #                     levels = lvls))
  # }
  # else {
  #   factor(apply(as.data.frame(colData(dds)[, intgroup, drop = FALSE]), 
  #                1, paste, collapse = ":"))
  # }
  #  data <- data.frame(count = cnts + pc, group = as.integer(group)) # orig
  #  dataGKT <- data.frame(count = cnts, group = as.integer(group))
  #  dataGKT <- data.frame(count = counts(dds, normalized = normalized, replaced = replaced)[gene,], colData(dds)[intgroup])
  #  dataGKT <- data.frame(count = counts(dds, normalized = normalized, replaced = replaced)[gene,], group = as.factor(group))
  dataGKT <- data.frame(count = counts(dds, normalized = normalized, replaced = replaced)[gene,],  colData(dds)[var])
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
    #    return(data.frame(count = dataGKT$count, colData(dds)[var]))
    return(dataGKT)
  #    plot(data$group + runif(ncol(dds), -0.05, 0.05), data$count, 
  #       xlim = c(0.5, max(data$group) + 0.5), log = logxy, xaxt = "n", 
  #       xlab = xlab, ylab = ylab, main = main, ...)
  #  axis(1, at = seq_along(levels(group)), levels(group))
  (ggplot(dataGKT, aes(colData(dds)[[var]], count))
    + geom_point()
    #    + geom_violin()
    #    + geom_boxplot(width = 0.2,outlier.shape = NA)
    #    + geom_sina()
    + xlab(paste(intgroup,collapse = ":"))
    + theme(axis.text.x = element_text(angle = 90))
    + labs(title = main)
  )
}

linePlotsGKT <- function (dds, gene, subject, connect, group, main, normalized = TRUE, replaced = FALSE)
{
  dataGKT <- data.frame(count = counts(dds, normalized = normalized, replaced = replaced)[gene,], colData(dds))
  if (missing(main)) {
    main <- if (is.numeric(gene)) {
      rownames(dds)[gene]
    }
    else {
      gene
    }
  }
  ylab <- ifelse(normalized, "normalized count", "count")
  (gglineplotf <- ggplot(dataGKT, aes(x=colData(dds)[[connect]], y=count))
    + scale_color_manual(values = c("red","purple4"))
    + guides(color=guide_legend(title=group))
    + geom_line(aes(color=colData(dds)[[group]], group = colData(dds)[[subject]]))
    #  + geom_point(aes(color=Group))
    #     + geom_text(aes(label=SubjectID), size = 2)
    #    + facet_wrap(~GeneSymbol, scales = "free_y")
    + theme_bw()
    + theme(axis.text.x = element_text(size = 6))
    + ylab(ylab)
    + xlab(connect)
    + labs(title = main)
    #    + scale_y_log10()
  )
  
}

pcaDataGKT <- function (object, intgroup = "condition", ntop = 500) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(pca$x, group = group, 
                  intgroup.df, name = colnames(object))
  attr(d, "percentVar") <- percentVar
  return(d)
}

pcaPlotGKT <- function(object, intgroup = "condition", xpc = 1, ypc = 2, ntop = 500)
{
  d <- pcaDataGKT(object, intgroup, ntop)
  percentVar <- round(100 * attr(d, "percentVar"))
  ggplot(d,aes_string(x = names(d)[xpc], y = names(d)[ypc])) + labs(x=paste0(names(d)[xpc],": ", percentVar[xpc], "% variance"), y=paste0(names(d)[ypc],": ", percentVar[ypc], "% variance"))
}

###sigGKT <- function(object, alpha = 0.05, lfcThreshold) #lfc single: abs threshold, 2vector, up and down thresholds

pForAlpha <- function(object, alpha = 0.05)
{ 
  clean <- as.data.frame(object) %>% 
    filter_all(all_vars(!is.na(.))) %>% 
    arrange(padj)
  
  if(clean[1,]$padj > alpha) { return(0.0) }
  
  highnot <- clean %>% 
    filter(padj > alpha) %>%
    head(1) %>%
    select(pvalue)
  lowsig <- clean %>% 
    filter(padj < alpha) %>% 
    tail(1) %>% 
    select(pvalue)
  mean(highnot$pvalue,lowsig$pvalue)
}

sigGenes <- function(object, alpha = 0.05, lfcThreshold = 1)
{
  as.data.frame(object) %>% 
    mutate(genesymbol = rownames(.)) %>% 
    filter_all(all_vars(!is.na(.))) %>% 
    filter(lfcSE < 1) %>% 
    filter(padj <= alpha) %>%
    filter(abs(log2FoldChange) >= lfcThreshold) %>%
    select(genesymbol)
}

volcanoPlotDataGKT <- function(object, alpha = 0.05, lfcThreshold = 1)
{
  ### should check result object for alpha and lfcThreshold values if there
  ### maybe default lfcThreshold should be 0
  d <- as.data.frame(object) %>%
    filter_all(all_vars(!is.na(.))) %>%
    filter(lfcSE < 1) %>%
    mutate(sig = if_else(padj <= alpha & abs(log2FoldChange) >= lfcThreshold, "significant", "not significant"))
  
  d$sig <- factor(d$sig, levels = c("significant","not significant"))
  return(droplevels(d))
  #  ggplot(d, aes(x=log2FoldChange, y=-log10(pvalue)))
  #  ggplot(d, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point(aes(color = sig), size = 0.3) + scale_color_manual(values = c("red","black"))
}

volcanoPlotGKT <- function(object, alpha = 0.05, lfcThreshold = 1)
{
  ### should check result object for alpha and lfcThreshold values if there
  ### maybe default lfcThreshold should be 0
  d <- as.data.frame(object) %>%
    filter_all(all_vars(!is.na(.))) %>%
    filter(lfcSE < 1) %>%
    mutate(sig = if_else(padj <= alpha & abs(log2FoldChange) >= lfcThreshold, "significant", "not significant"))
  
  d$sig <- factor(d$sig, levels = c("significant","not significant"))
  
  ggplot(d, aes(x=log2FoldChange, y=-log10(pvalue)))
  #  ggplot(d, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point(aes(color = sig), size = 0.3) + scale_color_manual(values = c("red","black"))
}

volcanoPointGKT <- geom_point(aes(color = sig), size = 0.7)
volcanoColorGKT <- scale_color_manual(values = c("red","black"))
volcanoColorNoSigGKT <- scale_color_manual(values = c("black"))

stdVolcanoPlotGKT <- function(object, alpha = 0.05, lfcThreshold = 1, main_title) 
{
  pline = round(-log10(pForAlpha(object, alpha)), 2)
  #  as.numeric(format(round(min(res[res$sigGenes == 1 | res$sigGenes == -1,]$negLogPval), 2), nsmall = 2))
  d = volcanoPlotDataGKT(object, alpha, lfcThreshold)
  p = ggplot(d, aes(x=log2FoldChange, y=-log10(pvalue))) +
    #    volcanoPlotGKT(object,alpha,lfcThreshold) +
    volcanoPointGKT +
    #    volcanoColorGKT + 
    scale_x_continuous(limits = c(-9,9), breaks = c(-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10)) +
    scale_y_continuous(limits = c(-1,50), breaks = c(0,pline,5,10,15,20,25,30,35,40,45)) +
    theme_bw() +
    theme(
      text = element_text(size = 12, family = "ArialMT", face = "bold"),
      axis.text = element_text(size = 10), axis.title = element_text(size = 12), 
      plot.title = element_text(hjust = 0.5, size = 12), 
      panel.grid.minor = element_blank(), 
      legend.title = element_blank(), legend.position = c(1,1), legend.justification = c(1,1), 
      legend.background = element_blank(), legend.key = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    geom_vline(xintercept = 1, linetype = "dotdash", color = "#636363") +
    geom_vline(xintercept = -1, linetype = "dotdash", color = "#636363") +
    geom_hline(yintercept = pline, linetype = "dotdash", color = "#636363") +
    labs(x=expression(fold~change~(log[2])), y=expression(p~value~(-log[10])), title=main_title)
  if(length(levels(d$sig)) > 1) { return(p + volcanoColorGKT) }
  else { return(p + volcanoColorNoSigGKT) }
}

maxMinFilter <- function(object, intgroup = "condition")
{
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  
  group <- if (length(intgroup) > 1) {
    interaction(lapply(intgroup, function(factorname) colData(ddsMoTime20034)[[factorname]]), drop = TRUE)
    #    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  #  if_else(rowMaxs(sapply(c("d0","d30"), function(lvl) rowMins(counts(ddsMoTime20034, normalize = TRUE)[,ddsMoTime20034$Timepoint == lvl]))) > 0,rowMeans2(counts(ddsMoTime20034, normalized = TRUE)),0)
  if_else(rowMaxs(sapply(levels(group), function(lvl) rowMins(counts(object, normalize = TRUE)[,group == lvl]))) > 0,rowMeans2(counts(object, normalized = TRUE)),0)
}

maxMinFilter <- function(object, intgroup = "condition", comp = c("ctrl","trmt"), thresh = 0)
{
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  
  group <- if (length(intgroup) > 1) {
    interaction(lapply(intgroup, function(factorname) colData(ddsMoTime20034)[[factorname]]), drop = TRUE)
    #    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  if (!all(comp %in% levels(group))) {
    stop("the argument 'comp' should specify levels of intgroup")
  }
  #  if_else(rowMaxs(sapply(c("d0","d30"), function(lvl) rowMins(counts(ddsMoTime20034, normalize = TRUE)[,ddsMoTime20034$Timepoint == lvl]))) > 0,rowMeans2(counts(ddsMoTime20034, normalized = TRUE)),0)
  if_else(rowMaxs(sapply(comp, function(lvl) rowMins(counts(object, normalize = TRUE)[,group == lvl]))) > thresh,rowMeans2(counts(object, normalized = TRUE)),0)
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

as.annot.data.frame.mutate <- function(object, var = "gene_id", subset = "all", relabel = "", newlabel = "", select = c("gene","baseMean","log2FoldChange","pvalue","padj"), mutate = "Contrast")
{
  if(mutate != "") {
    select <- c(select,"significant")
    df <- rownames_to_column(as.data.frame(object), var = var) %>% 
      mutate(significant = factor(if_else(is.na(padj),"removed",if_else(padj<object@metadata$alpha,"significant","not_significant")), levels = c("significant","not_significant","removed"))) %>%
      rename_at(vars(2:7), function(x) paste0(names(object)[1:6], " (", object@elementMetadata$description[1:6], ")"))
  } else {
    df <- rownames_to_column(as.data.frame(object), var = var) %>% rename_at(vars(2:7), function(x) paste0(names(object)[1:6], " (", object@elementMetadata$description[1:6], ")"))
  }
  
  if(subset != "all") {
    subset
    names(df)[2:7]
    #    df <- df %>% rename_at(vars(2:7), function(x) gsub(": ", paste0(": ",subset," "), gsub("all",subset, names(df)[2:7])))
    df <- df %>% rename_at(vars(2:7), function(x) gsub("all", subset, names(df)[2:7]))
  }
  if(relabel != "" & newlabel != "") {
    df <- df %>% rename_at(vars(2:7), function(x) gsub(relabel, newlabel, names(df)[2:7], fixed = TRUE))
    relabel_nous <- gsub("_"," ",relabel)
    df <- df %>% rename_at(vars(2:7), function(x) gsub(relabel_nous, newlabel, names(df)[2:7], fixed = TRUE))
    df <- df %>% rename_at(vars(2:7), function(x) gsub("p-values", paste0("p-values: ",newlabel), names(df)[2:7], fixed = TRUE))
  }
  if(length(select) > 1) {
    df <- df %>% select(contains(select))
  }
  if(mutate != "") {
    df <- df %>% mutate(subset = subset)
    mut_val <- if_else(!is.na(newlabel) && newlabel != "",newlabel, names(df)[3])
    df <- df %>% mutate(Contrast = mut_val)
    #    df <- df %>% mutate(significant = if_else(is.na(as.name(names(df)[5])),"removed",if_else(as.name(names(df)[5])<object@metadata$alpha,"significant",names(df)[5])))
    #    df$significant <- factor(df$significant, levels = c("significant","not_significant","removed"))
  }
  df
}

as.annot.data.frame <- function(object, var = "gene_id", subset = "all", relabel = "", newlabel = "", select = c("gene","baseMean","log2FoldChange","pvalue","padj"))
{
  df <- rownames_to_column(as.data.frame(object), var = var) %>% rename_at(vars(2:7), function(x) paste0(names(object)[1:6], " (", object@elementMetadata$description[1:6], ")"))
  if(subset != "all") {
    subset
    names(df)[2:7]
    #    df <- df %>% rename_at(vars(2:7), function(x) gsub(": ", paste0(": ",subset," "), gsub("all",subset, names(df)[2:7])))
    df <- df %>% rename_at(vars(2:7), function(x) gsub("all", subset, names(df)[2:7]))
  }
  if(relabel != "" & newlabel != "") {
    df <- df %>% rename_at(vars(2:7), function(x) gsub(relabel, newlabel, names(df)[2:7], fixed = TRUE))
    relabel_nous <- gsub("_"," ",relabel)
    df <- df %>% rename_at(vars(2:7), function(x) gsub(relabel_nous, newlabel, names(df)[2:7], fixed = TRUE))
    df <- df %>% rename_at(vars(2:7), function(x) gsub("p-values", paste0("p-values: ",newlabel), names(df)[2:7], fixed = TRUE))
  }
  if(length(select) > 1) {
    df <- df %>% select(contains(select))
  }
  df
}

annotCount <- function(object, normalized = FALSE)
{
  origCols <- colnames(object)
  valueLabel <- 
    rownames_to_column(as.data.frame(colData(object)),var = "SampleID") %>% left_join(gather_(rownames_to_column(as.data.frame(counts(object,normalized = normalized)), var = "EnsemblID"), "SampleID", ifelse(normalized,"NormalizedCount","Count"), origCols))
}


alphaPvalue <- function(object, alpha = NULL)
{
  #  object@metadata$alpha
  #  hasArg(alpha)
  #  alphaval <- if_else(hasArg(alpha),as.double(alpha),as.double(object@metatdata$alpha))
  if(hasArg(alpha)) {
    alphaval <- alpha
  } else {
    alphaval <- object@metadata$alpha
  }
  #  wupper <- {x<-resB2_filt_g1456_STIPosVsNeg_inHIVPos$padj;x[x<alphaval]<-max(x);which.min(x)}
  wupper <- {x<-object$padj;x[x<alphaval]<-max(x);which.min(x)}
  alphaupper <- object[wupper,]$padj
  pupper <- object[wupper,]$pvalue
  #  wlower <- {x<-resB2_filt_g1456_STIPosVsNeg_inHIVPos$padj;x[x>alphaval]<-min(x);which.max(x)}
  wlower <- {x<-object$padj;x[x>alphaval]<-min(x);which.max(x)}
  alphalower <- object[wlower,]$padj
  plower <- object[wlower,]$pvalue
  
  plower + (alphaval-alphalower) * (pupper-plower) / (alphaupper - alphalower)
}

# volcano plots with EnhancedVolcano

DESeq2ResultsEnhancedVolcano <- function(object, fcThresh = 1, padjThresh = 0.05, baseMeanThresh = 20, title = "volcano plot", subtitle = "comparison") 
{
  volcanoData <- Mmul10Annot %>% right_join(as.data.frame(object) %>% rownames_to_column(var = "gene_id"))  %>% filter(gene_id %in% Mmul10ProteinCodingByENSID$gene_id & baseMean > baseMeanThresh & lfcSE < 1) %>% arrange(-log2FoldChange)
  pThresh <- mean((as.data.frame(object) %>% rownames_to_column(var = "gene_id") %>% filter(padj<0.05) %>% arrange(-padj) %>% select(pvalue))$pvalue[1],(as.data.frame(object) %>% rownames_to_column(var = "gene_id") %>% filter(padj>0.05) %>% arrange(padj) %>% select(pvalue))$pvalue[1])
  volcanoDataFiltered <- volcanoData %>% filter(padj < padjThresh & abs(log2FoldChange) >= fcThresh)
  upCount <- count(volcanoDataFiltered %>% filter(log2FoldChange > 0))
  downCount <- count(volcanoDataFiltered %>% filter(log2FoldChange < 0))
  keyvals <- ifelse(volcanoData$padj < padjThresh,ifelse(volcanoData$log2FoldChange > fcThresh,'#d95f02',ifelse(volcanoData$log2FoldChange < -fcThresh,'#1b9e77','grey')),'grey')
  keyvals[is.na(keyvals)] <- 'grey'
  names(keyvals)[keyvals == '#d95f02'] <- "Up"
  names(keyvals)[keyvals == '#1b9e77'] <- "Down"
  names(keyvals)[keyvals == 'grey'] <- "NS"
  EnhancedVolcano(volcanoData, lab = volcanoData$gene_symbol, x = 'log2FoldChange', y = 'pvalue', pCutoff = pThresh, FCcutoff = fcThresh, selectLab = volcanoDataFiltered$gene_symbol[c(1:31, length(volcanoDataFiltered$gene_symbol) - (20:0))], xlim = c(-4,8), ylim = c(0,35), axisLabSize = 7, title = title, subtitle = str_remove(object@elementMetadata$description[2],".*: "), titleLabSize = 12, subtitleLabSize = 8, caption = "", captionLabSize = 0, legendPosition = 'none', labSize = 3, labhjust = 0.5, pointSize = 0.5, colCustom = keyvals, labCol = "black",  labFace = 'bold', hlineWidth = 0.05, vlineWidth = 0.05, borderWidth = 0.1, boxedLabels = FALSE, drawConnectors = TRUE, widthConnectors = 0.1, lengthConnectors = unit(0.001, "npc"), colConnectors = "grey") + annotate("text", label = paste0("Down=",downCount), size = 3, x = -1.5, y=0, hjust = 0.5) + annotate("text", label = paste0("Up=",upCount), size = 3, x = 1.5, y=0, hjust = 0.5) + theme(plot.title = element_text(hjust = 0.5, face = "plain", margin = margin(b = 0)), plot.caption = element_blank(), plot.margin = margin(2,2,2,2), axis.ticks = element_line(size = 0.1), axis.ticks.length = unit(2, "pt"), axis.title.x = element_text(margin = margin(0,0,0,0)), axis.text.x = element_text(hjust = 0.5, margin = margin(1,0,0,0)), axis.title.y = element_text(margin = margin(0,0,0,0)), axis.text.y = element_text(vjust = 0.5, margin = margin(0,0,0,1)), panel.grid = element_line(size = 0.1))
}

## GSEA stuff

phenorderGSEArun <- function (gsearundir)
{
  #  phenorder <- unlist(str_split(dir(gsearundir, glob2rx("ranked_gene_list*.tsv")),"_"))[c(4,6)]
  phencomp <- unlist(str_split(dir(gsearundir, glob2rx("ranked_gene_list*.tsv")),"_"))
  vsindex <- which(phencomp == "versus")
  phenNumerator <- paste(phencomp[4:(vsindex-1)],collapse = "_")
  phenDenominator <- paste(phencomp[(vsindex+1):(length(phencomp)-1)],collapse = "_")
  phenorder <- c(phenNumerator,phenDenominator)
  phenorder
}

resultsGSEArun <- function (gsearundir, runlabel = "")
{
  #  phenorder <- unlist(str_split(dir(gsearundir, glob2rx("ranked_gene_list*.tsv")),"_"))[c(4,6)]
  phenorder <- phenorderGSEArun(gsearundir)
  runTables <- dir(gsearundir,glob2rx("gsea_report_for*.tsv"))
  #  phenotypes <- unlist(lapply(str_split(runTables,"_"),"[",4))
  phenotypes <- str_remove(str_remove(dir(gsearundir,glob2rx("gsea_report_for*.tsv")),"_[0-9]*.tsv"), "gsea_report_for_")
  grep(phenorder[1],runTables,value = TRUE)
  compPath <- paste0(gsearundir,"/",grep(paste0("_",phenorder[1],"_"),runTables,value = TRUE))
  refPath <- paste0(gsearundir,"/",grep(paste0("_",phenorder[2],"_"),runTables,value = TRUE))
  compTable <- read_tsv(file = compPath, col_types = "c--dddddddc-")
  refTable <- read_tsv(file = refPath, col_types = "c--dddddddc-")
  resTable <- compTable %>% bind_rows(refTable) %>% dplyr::rename(GENESET = NAME)
  resTable$`NOM p-val` <- ifelse(resTable$`NOM p-val` == 0,0.001,resTable$`NOM p-val`)
  resTable <- resTable %>% mutate(topPhenotype = phenorder[1], bottomPhenotype = phenorder[2], .before = 1)
  if(length(runlabel)>0) resTable <- resTable %>% mutate(run = runlabel, .before = 1)
  # if(length(runlabel)>0) compTable %>% bind_rows(refTable) %>% mutate(run = runlabel, topPhenotype = phenorder[1], bottomPhenotype = phenorder[2], .before = 1)
  # else compTable %>% bind_rows(refTable) %>% mutate(topPhenotype = phenorder[1], bottomPhenotype = phenorder[2], .before = 1)
  # compTable %>% bind_rows(refTable) %>% mutate()
  # if(length(runid>0)) label = paste0(runid," ")
  # label <- paste0("(",label,phenorder[1],"vs",phenorder[2],")")
  # compTable %>% bind_rows(refTable) %>% rename_with( ~ paste0(.x,label), .cols = 3:9)
  #    compTable %>% bind_rows(refTable) %>% rename_with( ~ paste0(.x," (",phenorder[1],"vs",phenorder[2],")"), .cols = 3:9)
  #%>% bind_rows(read_tsv(file = paste0(gseaPath,"CD4Drop2_EFMvsMONO_H.Gsea.1598643059317/gsea_report_for_MpONO_1598643059317.tsv"), col_types = "c--dddddddc-")) %>% rename_with( ~ paste0(.x," (EFMvsMONO)"), .cols = 3:9)
}
# resultsGSEArun <- function (gsearundir, runlabel = "")
# {
#     phenorder <- unlist(str_split(dir(gsearundir, glob2rx("ranked_gene_list*.tsv")),"_"))[c(4,6)]
#     runTables <- dir(gsearundir,glob2rx("gsea_report_for*.tsv"))
#     phenotypes <- unlist(lapply(str_split(runTables,"_"),"[",4))
#     grep(phenorder[1],runTables,value = TRUE)
#     compPath <- paste0(gsearundir,"/",grep(phenorder[1],runTables,value = TRUE))
#     refPath <- paste0(gsearundir,"/",grep(phenorder[2],runTables,value = TRUE))
#     compTable <- read_tsv(file = compPath, col_types = "c--dddddddc-")
#     refTable <- read_tsv(file = refPath, col_types = "c--dddddddc-")
#     resTable <- compTable %>% bind_rows(refTable) %>% rename(GENESET = NAME)
#     resTable$`NOM p-val` <- ifelse(resTable$`NOM p-val` == 0,0.001,resTable$`NOM p-val`)
#     if(length(runlabel)>0) resTable <- resTable %>% mutate(run = runlabel, .before = 1)
#     resTable %>% mutate(topPhenotype = phenorder[1], bottomPhenotype = phenorder[2], .before = 1)
#     # if(length(runlabel)>0) compTable %>% bind_rows(refTable) %>% mutate(run = runlabel, topPhenotype = phenorder[1], bottomPhenotype = phenorder[2], .before = 1)
#     # else compTable %>% bind_rows(refTable) %>% mutate(topPhenotype = phenorder[1], bottomPhenotype = phenorder[2], .before = 1)
#     # compTable %>% bind_rows(refTable) %>% mutate()
#     # if(length(runid>0)) label = paste0(runid," ")
#     # label <- paste0("(",label,phenorder[1],"vs",phenorder[2],")")
#     # compTable %>% bind_rows(refTable) %>% rename_with( ~ paste0(.x,label), .cols = 3:9)
# #    compTable %>% bind_rows(refTable) %>% rename_with( ~ paste0(.x," (",phenorder[1],"vs",phenorder[2],")"), .cols = 3:9)
#     #%>% bind_rows(read_tsv(file = paste0(gseaPath,"CD4Drop2_EFMvsMONO_H.Gsea.1598643059317/gsea_report_for_MpONO_1598643059317.tsv"), col_types = "c--dddddddc-")) %>% rename_with( ~ paste0(.x," (EFMvsMONO)"), .cols = 3:9)
# }

setListGSEArun <- function (gsearundir)
{
  tsvs <- dir(gsearundir, glob2rx("*.tsv"))
  tsvs <- tsvs[-grep("gene_set_sizes.tsv",tsvs)]
  tsvs <- tsvs[-grep(glob2rx("gsea_report_for*.tsv"),tsvs)]
  tsvs <- tsvs[-grep(glob2rx("ranked_gene_list*.tsv"),tsvs)]
  tsvs <- tsvs[-grep("Symbol_to_probe_set_mapping_details.tsv",tsvs)]
  tsvs
}

atomresGSEAset <- function(geneset, gsearundir)
{
  settsv <- paste0(gsearundir,"/",geneset,".tsv")
  if(file.exists(settsv)) {
    read_tsv(settsv, col_types = "-cciddf-")
  }
}

atomgenelistGSEAset <- function(resGSEAset, leadingEdge = FALSE)
{
  if(leadingEdge) filter(resGSEAset, `CORE ENRICHMENT` == "Yes") %>%
    pull(SYMBOL)
  else
    pull(resGSEAset, SYMBOL)
}

resultsGSEAsetraw <- function (geneset, gsearundir)
{
  sapply(geneset, function(set)
  {
    if(file.exists(paste0(gsearundir,"/",set,".tsv"))) {
      read_tsv(paste0(gsearundir,"/",set,".tsv"), col_types = "-cciddf-")
    }
  },
  simplify = FALSE,
  USE.NAMES = TRUE
  )
}

resultsGSEAPrerankedsetraw <- function (geneset, gsearundir)
{
  sapply(geneset, function(set)
  {
    if(file.exists(paste0(gsearundir,"/",set,".tsv"))) {
      read_tsv(paste0(gsearundir,"/",set,".tsv"), col_types = "-ciddf-")
    }
  },
  simplify = FALSE,
  USE.NAMES = TRUE
  )
}

genelistGSEAset <- function(geneset, gsearundir)
{
  resultsGSEAsetraw(geneset, gsearundir)
}

resultsGSEAset <- function (gsearundir, geneset)
{
  phenorder <- phenorderGSEArun(gsearundir)
  if(file.exists(paste0(gsearundir,"/",geneset,".tsv"))) { read_tsv(paste0(gsearundir,"/",geneset,".tsv"), col_types = "-cciddf-") %>% rename_with(~ paste0(.x," (",phenorder[1],"vs",phenorder[2],")"), .cols = 3:6)
  }
}

LE_GSEA_RMS <- function (geneset, gsearundir)
{
  sets <- resultsGSEAsetraw(geneset, gsearundir)
  sets <- sets[!sapply(sets,is.null)]
  sapply(names(sets), function(name) 
  {
    #name = names(sets)[[i]]
    filter(sets[[name]], `CORE ENRICHMENT` == "Yes") %>%
      select(c(SYMBOL, `RANK METRIC SCORE`)) %>%
      transmute(SYMBOL = SYMBOL,
                "{name}" := `RANK METRIC SCORE`)
  },
  simplify = FALSE,
  USE.NAMES = TRUE
  )
  # if(file.exists(paste0(gsearundir,"/",geneset,".tsv"))) { 
  #   read_tsv(paste0(gsearundir,"/",geneset,".tsv"),
  #            col_types = "-cciddf-") %>% 
  #     filter(`CORE ENRICHMENT` == "Yes") %>%
  #     select(c(SYMBOL,`RANK METRIC SCORE`)) %>%
  #     transmute(SYMBOL = SYMBOL,
  #               "{geneset}" := `RANK METRIC SCORE`)
  # }
}

LE_GSEA_geneList <- function (geneset, gsearundir)
{
  sets <- resultsGSEAsetraw(geneset, gsearundir)
  sets <- sets[!sapply(sets,is.null)]
  sapply(names(sets), function(name) 
  {
    #name = names(sets)[[i]]
    filter(sets[[name]], `CORE ENRICHMENT` == "Yes") %>%
      select(c(SYMBOL)) %>%
      transmute(SYMBOL = SYMBOL)
  },
  simplify = FALSE,
  USE.NAMES = TRUE
  )
  # if(file.exists(paste0(gsearundir,"/",geneset,".tsv"))) { 
  #   read_tsv(paste0(gsearundir,"/",geneset,".tsv"),
  #            col_types = "-cciddf-") %>% 
  #     filter(`CORE ENRICHMENT` == "Yes") %>%
  #     select(c(SYMBOL,`RANK METRIC SCORE`)) %>%
  #     transmute(SYMBOL = SYMBOL,
  #               "{geneset}" := `RANK METRIC SCORE`)
  # }
}

LE_GSEAPreranked_RMS <- function (geneset, gsearundir)
{
  sets <- resultsGSEAPrerankedsetraw(geneset, gsearundir)
  sets <- sets[!sapply(sets,is.null)]
  sapply(names(sets), function(name) 
  {
    #name = names(sets)[[i]]
    filter(sets[[name]], `CORE ENRICHMENT` == "Yes") %>%
      select(c(SYMBOL, `RANK METRIC SCORE`)) %>%
      transmute(SYMBOL = SYMBOL,
                "{name}" := `RANK METRIC SCORE`)
  },
  simplify = FALSE,
  USE.NAMES = TRUE
  )
}

LE_GSEAPreranked_geneList <- function (geneset, gsearundir)
{
  sets <- resultsGSEAPrerankedsetraw(geneset, gsearundir)
  sets <- sets[!sapply(sets,is.null)]
  sapply(names(sets), function(name) 
  {
    #name = names(sets)[[i]]
    filter(sets[[name]], `CORE ENRICHMENT` == "Yes") %>%
      select(c(SYMBOL)) %>%
      transmute(SYMBOL = SYMBOL)
  },
  simplify = FALSE,
  USE.NAMES = TRUE
  )
}

saveGSEADotPlot <- function(object, name)
{
  gseaDotPlot <- object %>% filter(GENESET %in% unlist(c(GSEA_MET[,name]))) %>% ggplot(aes(x=topPhenotype, y=GENESET, color=NES, size=1-log(`NOM p-val`))) + geom_point(na.rm = TRUE) + scale_color_gradientn(colors = c("blue","grey75","grey75","red"), limits = c(-2.8,2.8), breaks = c(-2,-1,0,1,2)) + theme_bw() + theme(axis.text.y = element_text(angle = 30, vjust = 0.5, hjust=1, size = 4, face = "bold"), axis.text.x = element_text(angle = 30)) + scale_radius("1-log(NOM p-val)", range = c(1,(12-0.3*length(unlist(c(na.omit(GSEA_MET[,name])))))), breaks =  c(2,4,6,8), limits = c(1,8.2), labels = c("0.37","0.05","0.007","<0.001")) + facet_grid(cols = vars(run)) + ggtitle(paste(name," Genesets of Interest"))
  ggsave(gseaDotPlot, filename = paste0("resGSEA_",name,".pdf"), units = "in", width = 9, height = 6, scale = 1)
}

getGeneLists <- function(genesetList, source, contrast) {
  genelists <- lapply(unique(genesetList), function(GENESET) {
    genelist <- resultsGSEAset(paste0(gseaPath,filter(studyrunsGSEA$runTable, Source == source & Contrast == contrast) %>% select(RunDir) %>% unlist()),GENESET) %>% .$SYMBOL
    if(length(genelist) > 50) {
      genelist <- resultsGSEAset(paste0(gseaPath,filter(studyrunsGSEA$runTable, Source == source & Contrast == contrast) %>% select(RunDir) %>% unlist()),GENESET) %>% filter(`CORE ENRICHMENT (2dpivsminus5dpi)` == "Yes") %>% .$SYMBOL
    }
    genelist
  })
  names(genelists) <- genesetList
  genelists
}

saveGSEAHeatmap <- function(plotset, source, contrast) {
  pdf(paste0(plotset,"_Heatmap.pdf"), width = 9, height = 6)
  genelists <- getGeneLists(plotset, source, contrast)
  ht = Heatmap(rlogFromMedBaselineBALPathP20097[genelists[[plotset]],], col = heatscale1, column_title = plotset,name = "log2 foldChange\nfrom median Baseline", heatmap_legend_param = list(title_gp = gpar(fontsize = 6)), column_labels = rldBALPathP20097forPlots$SubjectID, column_names_rot = 0, column_names_centered = TRUE, row_order = genelists[[plotset]], column_order = order(rldBALPathP20097forPlots$Timepoint), row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 5), top_annotation = HeatmapAnnotation(foo = anno_block(labels = c("Baseline","2 dpi","4 dpi","7 dpi","Necropsy"))), column_split = rldBALPathP20097forPlots$Timepoint)
  #  heatmap_GeneLists[[plotset]]
  draw(ht)
  dev.off()
  
}

## functions to compute RLE plots and get stats for filtering

relLogExp <- function(counts) {
  y = log(counts)
  y - rowMedians(y)
  #  log(counts) - log(rowMedians(counts))
}

relLogExpTibble <- function(relLogExp) {
  as_tibble(relLogExp) %>% pivot_longer(everything(), names_to = "SampleID", values_to = "relLogExp")
}

relLogExpPlot <- function(relLogExpTib) {
  ggplot(relLogExpTib, aes(x=SampleID, y=relLogExp))
}

stat_relLogExpGKT <- function(relLogExp) {
  relLogExpPlot(relLogExp) +
    geom_hline(yintercept = 0, color = "red") +
    geom_violin(draw_quantiles = c(0.25,0.75), trim = TRUE, color = "lightgreen", alpha = 0.1) +
    geom_boxplot(alpha = 0) +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() +
    theme(axis.text.x = element_text(vjust = 0.5, angle = 90))
}

str_wrap_balance <- function(string, width = 80, indent = 0, exdent = 0,USE.NAMES = FALSE) {
  out <- str_wrap(string, width, indent, exdent)
  vapply(out, function(string,width,indent,exdent) {
    wraps <- str_count(string,"\n")
    if(wraps > 0 && width > 1) {
      bwidth <- width
      repeat {
        bwidth <- bwidth - 1
        bstring <- str_wrap(string, bwidth, indent, exdent)
        bwraps <- str_count(bstring,"\n")
        if(bwraps > wraps || bwidth <= 1) break
        string <- bstring
      }
    }
    string
  },character(1),width,indent,exdent,USE.NAMES = USE.NAMES)
}

wrap_underscore_strings <- function(string, width) {
  str_replace_all(str_wrap(str_replace_all(string,"_"," "),width)," ","_")
}

wrap_underscore_strings_balance <- function(string, width) {
  str_replace_all(str_wrap_balance(str_replace_all(string,"_"," "),width)," ","_")
}
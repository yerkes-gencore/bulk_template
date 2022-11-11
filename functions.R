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


wrap_underscore_strings_balance <- function(string, width) {
  str_replace_all(str_wrap_balance(str_replace_all(string,"_"," "),width)," ","_")
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

maxMinFilter <- function(object, intgroup = "condition", comp = c("ctrl","trmt"), thresh = 0)
{
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  
  # group <- if (length(intgroup) > 1) {
  #   interaction(lapply(intgroup, function(factorname) colData(ddsMoTime20034)[[factorname]]), drop = TRUE)
  #   #    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  # }
  # else {
  group <- colData(object)[[intgroup]]
  # }
  if (!all(comp %in% levels(group))) {
    stop("the argument 'comp' should specify levels of intgroup")
  }
  #  if_else(rowMaxs(sapply(c("d0","d30"), function(lvl) rowMins(counts(ddsMoTime20034, normalize = TRUE)[,ddsMoTime20034$Timepoint == lvl]))) > 0,rowMeans2(counts(ddsMoTime20034, normalized = TRUE)),0)
  if_else(rowMaxs(sapply(comp, function(lvl) rowMins(counts(object, normalize = TRUE)[,group == lvl]))) > thresh,rowMeans2(counts(object, normalized = TRUE)),0)
}


##### Derrik functions

### makes a heatmap of the given list of genes, separating samples slightly by group variable
heatmap_from_genelist <- function(geneList, baseline_grouping, baseline, data=analysis$rldDrop){
  ### makes a heatmap of the given list of genes, separating samples slightly by group variable
  ## data should be of type DESeqTransform
  hmap <- data[rownames(data) %in% geneList,
               !colnames(data) %in% analysis$config$dropSamples]
  slices <- as.numeric(as.factor(hmap@colData$Group))
  labels <- levels(as.factor(hmap@colData$Group))[unique(slices)]
  baseline <- rowMedians(assay(hmap[,as.character(hmap@colData[[baseline_grouping]]) %in% baseline]))
  hmap <- assay(hmap) - baseline
  Heatmap(hmap, show_row_names = TRUE, heatmap_legend_param = list(title="log2 fold\ndifference\nfrom\nmedian\nweek 0\nexpression"),border="black",column_order=colnames(hmap),
          #  width = ncol(hmap)*unit(5, "mm"), 
          # height = nrow(hmap)*unit(5, "mm"),
          rect_gp = gpar(color="black"),
          column_title=" ",
          cluster_rows=FALSE,
          column_split=slices,
          top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(col=NA),
                                                              #,fill = 2:length(unique(exp.design$Group))+1),
                                                              labels = labels, 
                                                              labels_gp = gpar(col = "black", fontsize = 10), labels_rot=70, height=unit(2,"cm"))),
          column_gap=unit(2, "mm"),
          col=colorRamp2(c(-4, 0, 4), c("blue", "white", "red")))
}


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


### Mapping results bins plot
mapping_plot <- function(mapBins=analysis$mapBins, title=FALSE){
  data <- rownames_to_column(as_tibble(mapBins, rownames = NA), var = "map_result") %>%
  pivot_longer(!map_result, names_to = "SampleID", values_to = "count") 
  data$map_result <- unlist(lapply(data$map_result, function(x)toTitleCase(str_split(x, '_')[[1]][2])))
  data$map_result <- factor(data$map_result, levels = c("Unmapped","Multimapping","noFeature","Ambiguous","Identified"))
  #mutate(SampleID = str_sub(SampleID,end = -18)) %>%
  ggplot(data, aes(x = SampleID, y = count, fill = map_result)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c('Unmapped' = 'grey', 'Multimapping' = 'orange', 'noFeature' = 'red', 'Ambiguous' = 'yellow', 'Identified' = '#00CC33' )) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  guides(fill=guide_legend(title="Map Result"))  + aes(fct_inorder(SampleID))+ labs(x="Sample", y='Number of reads') +
    (if (title){ ggtitle(paste0(analysis$config$analysis," Mapping")) } else NULL)
}

### RLE plots
RLE_plot <- function(data, title){
  as_tibble(data, rownames = NA) %>% pivot_longer(everything(), names_to = "Sample", values_to = "RLE") %>% 
    ggplot(aes(x=Sample,y=RLE)) +
    geom_hline(yintercept = 0, color = "red") + geom_violin(draw_quantiles = c(0.25,0.75), trim = TRUE, color = "lightgreen", alpha = 0.1) +
    geom_boxplot(alpha = 0) + theme_bw() + theme(axis.text.x = element_text(vjust = 0.5,angle = 90),axis.title.x = element_blank(),aspect.ratio = 0.55) + 
    ggtitle(title) + aes(fct_inorder(Sample))#+ scale_y_continuous(limits = c(-9,3.25), expand = c(0,0))
}

### PCA plot from config, dependant on pcaPlotGKT data, just adds some overlays
PCA_plot_from_config <- function(analysis=analysis){
  pcaPlotGKT(analysis$vst,
             intgroup = names(colData(analysis$vst)),
             xpc = 1, ypc = 2) +
    #  geom_path(aes(group = SubjectID)) 
    geom_point(aes(color = (if (is.null(analysis$config$pcaMapping$color)) NULL else .data[[analysis$config$pcaMapping$color]]),
                   shape = (if (is.null(analysis$config$pcaMapping$shape)) NULL else .data[[analysis$config$pcaMapping$shape]])),
               size = 5) +
    labs(color=analysis$config$pcaMapping$color, shape=analysis$config$pcaMapping$shape) +
    (if (is.null(analysis$config$pcaMapping$label)) NULL else geom_text_repel(aes(label = .data[[analysis$config$pcaMapping$label]]),
                    size = 4, hjust = 0.5, vjust = -0.5, alpha=0.5)) +
    scale_x_continuous(expand = c(0.5,0)) +
    theme_bw() + ggtitle(paste0(analysis$config$analysis," PCA")) +
    (if (!is.null(analysis$config$pcaMapping$path)){
      geom_path(aes(linetype=.data[[analysis$config$pcaMapping$path]]))
    }) +
    # geom_path(aes(linetype=(if (is.null(analysis$config$pcaMapping$path)) NULL else .data[[analysis$config$pcaMapping$path]]))) +
    theme(legend.key.width = unit(1.2, "cm")) +
    #labs(color="Weeks") +
    (if (!is.null(analysis$config$pcaMapping$ellipse)){
      stat_ellipse(aes(color=.data[[analysis$config$pcaMapping$ellipse]]), type="norm", level=0.67)
    })+
    theme(text = element_text(size=10)) #, arrow=arrow(ends="last", type="closed", length=unit(0.1, "inches")))
  
  
}

#### DESeq2 results object
### Generate volcano plot from DESeq2 results table
### also references config GeneList to label genes of interest
generate_volplot <- function(result){
  ## Volcano plot
  volData <- result[!is.na(result$padj),]
  volData <- volData[order(volData$padj),]
  labs <- if (!is.null(analysis$config$geneList[[1]])) analysis$config$geneList else rownames(volData[1:20,])
  volplot <- (EnhancedVolcano(volData,
                              x = 'log2FoldChange',
                              y = 'padj',
                              lab = rownames(volData),
                              selectLab=labs,
                              drawConnectors = TRUE,
                              colConnectors = "lightgrey",
                              pCutoff = analysis$config$alpha,
                              FCcutoff = log2(1.3),
                              title =NULL,# paste0(numer_contrast, " vs ", denom_contrast),
                              caption=NULL,
                              subtitle=NULL
                              
  )) 
  return(volplot)
}

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

run_fgsea <- function(result, gmt.file=gmt.file, nperm=1000){
  result <- result %>% na.omit()
  tmp1 <- result$stat
  names(tmp1) <- rownames(result)
  tmp1 <- tmp1[!is.na(tmp1)]
  res <- fgseaSimple(pathways=gmt.file,
                     stats=tmp1,
                     nperm=nperm)
  return(res)
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

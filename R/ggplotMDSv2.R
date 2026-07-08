ggplotMDSv2 <- function(dge, sampleTable = dge$samples,
                        mode = "limma", deseqData = NULL,
                        sampleID = "sampleID", gene.selection = "common", 
                        dims = c(1,2), nDEGs = 20,
                        color = NULL, shape = NULL, size = 4,
                        ellipse = NULL, path = NULL,
                        show.labels = TRUE, label.size = 4, custom.labels = NULL, 
                        alpha = 1,...) {
  # mode = "limma" or "deseq"
  ## just changes where data comes from 
  # deseqData: from assay(bulk$deseq$vst)
  
  # Get mds data from edgeR::plotMDS()
  if(mode == "deseq"){
    mds_data <- limma::plotMDS(deseqData, top = 500, plot = FALSE, 
                               gene.selection = gene.selection, dim.plot = dims)
  } else{
    mds_data <- limma::plotMDS(dge, top = 500, plot = FALSE, 
                               gene.selection = gene.selection, dim.plot = dims)
  }
  
  # get degs
  pca_fit <- edgeR::voomLmFit(dge, design =
                                cbind(Intercept=1, 
                                      xaxis=mds_data$x, 
                                      yaxis=mds_data$y)) %>% limma::eBayes()
  xaxis_degs <- limma::topTable(pca_fit, coef="xaxis", number = nDEGs)
  yaxis_degs <- limma::topTable(pca_fit, coef="yaxis", number = nDEGs)
  xaxis_name <- paste0("PC",dims[1],"_DEGs")
  yaxis_name <- paste0("PC",dims[2],"_DEGs")
  
  # format coordinate table for ggplot
  mds_xy <- mds_data[c("x","y")]
  mds_xy[[sampleID]] <- colnames(dge)
  mds_xy <- dplyr::as_tibble(mds_xy) %>% dplyr::full_join(sampleTable, by = sampleID)

  # get percent variation
  x_varex <- round(mds_data$var.explained[dims[1]]*100, digits = 0)
  y_varex <- round(mds_data$var.explained[dims[2]]*100, digits = 0)
  
  if (!is.null(custom.labels)) {
    mds_xy <- mds_xy %>%
      dplyr::mutate(sampleID = ifelse(sampleID %in% custom.labels, sampleID, NA))
  }
  
  # fit environmental variables to the ordination, return table, plot, and DEGs as list
    fit <- vegan::envfit(mds_xy[,c("x","y")] ~ .,
                  data = dplyr::select(mds_xy, !c(x,y)), permutations = 999)
    fitTab <- rbind(
      cbind("variable" = names(fit$factors$r),
            "type" = "factor",
            "r" = round(fit$factors$r,3),
            "pval" = fit$factors$pvals),
      cbind("variable" = names(fit$vectors$r),
            "type" = "vector",
            "r" = round(fit$vectors$r,3),
            "pval" = fit$vectors$pvals)
    ) %>% as.data.frame() %>%
      dplyr::arrange(pval,desc(r)) %>%
      dplyr::filter(r !=1)
    
    ordPlot <- mds_xy %>%
      ggplot2::ggplot(aes(x = .data$x, y = .data$y)) +
      ggplot2::geom_point(
        aes(
          color = (
            if (!is.null(color)) { .data[[color]] } else { NULL }
          ),
          shape = (
            if (!is.null(shape)) { .data[[shape]] } else { NULL }
          ),),
        size = size,
        alpha = alpha) +
       (if (show.labels) { ggrepel::geom_text_repel(aes(label = .data[["sampleID"]]), box.padding = 0.5, na.rm = TRUE) } else { NULL }) +
       (if (!is.null(path)) { ggplot2::geom_path(aes(linetype = .data[[path]]), alpha = 0.5) } else { NULL }) +
       (if (!is.null(ellipse)) { ggplot2::stat_ellipse(aes(color = .data[[ellipse]]), type = "norm", level = 0.67)} else { NULL }) +
      ggplot2::labs(color = color, shape = shape) +
      xlab(paste0(mds_data$axislabel, " ", dims[1], " (", x_varex, "%)")) +
      ylab(paste0(mds_data$axislabel, " ", dims[2]," (", y_varex, "%)")) +
      ggplot2::theme_classic() +
      ggplot2::theme(aspect.ratio = 1) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::ggtitle(ifelse(gene.selection == "common", "PCA", "MDS"))
    
    allObj <- list("ord" = ordPlot,
                   "envfit" = fitTab)
    allObj[[xaxis_name]] <- xaxis_degs
    allObj[[yaxis_name]] <- yaxis_degs
    
    return(allObj)
  }

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

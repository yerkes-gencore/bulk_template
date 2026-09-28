#' Flag redundant pathways with significant leading edge overlap
#'
#' @param df Data frame containing results from GSEA
#' @param threshold Ratio of allowable overlap. Default is 1 where pathways 
#' who's leading edge genes are found to have 100% overlap with another pathway are flagged. 
#'
#' @returns data frame with flags column, ratio, and the parent pathway
#' @export
#'
#' @examples
flag_redundant_leading_edge <- function(df, threshold = 1) {
  
  n <- nrow(df)
  leading_edges <- lapply(df$leadingEdge, unique)
  
  overlap_mat <- matrix(
    0,
    nrow = n,
    ncol = n,
    dimnames = list(df$pathway, df$pathway)
  )
  
  for (x in seq_len(n)) {
    for (y in seq_len(n)) {
      
      if (x == y) next
      
      if(is.na(sign(df$NES[x])) | is.na(sign(df$NES[y]))) next
      
      # Only compare pathways with the same NES direction
      if (sign(df$NES[x]) != sign(df$NES[y])) next
      
      overlap_mat[x, y] <- sum(
        leading_edges[[x]] %in% leading_edges[[y]]
      ) / length(leading_edges[[x]])
    }
  }
  
  max_overlap <- apply(overlap_mat, 1, max)
  
  containing_pathway <- apply(overlap_mat, 1, function(x) {
    if (max(x) >= threshold) {
      names(which.max(x))
    } else {
      NA_character_
    }
  })
  
  df %>%
    mutate(
      leading_edge_overlap = max_overlap,
      leading_edge_contained_in = containing_pathway,
      redundant_leading_edge = max_overlap >= threshold
    )
}

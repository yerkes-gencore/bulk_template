getTopDegTbl <- function(
    resultsTable,            # generally bulk$degTab$deseq or bulk$degTab$limma
    contrast_id,             # single contrast, or loop through multiple
    groupID = "group",       # isn't actually used
    arrange_by = "pvalue",   # 'pvalue' or 'log2FoldChange'
    direction = "unequal",   # 'unequal', 'equal', 'up', or 'down'
    padj_cutoff = 0.05,      # adjusted p-value filter, only matters if arranging by log2FoldChange 
    slice_n = 50,            # how many top genes to pull
    filter_pattern = NULL)   # "^LOC//d+" to remove NCBI LOC genes, "^ENS[A-Z]{1,3}\\d+" for ensembl IDs
  {
  
  if (! arrange_by %in% c("pvalue", "log2FoldChange")) {
    stop("'arrange_by' parameter must be set to 'pvalue' (default) or 'log2FoldChange'")
  }
  
  if (arrange_by == "log2FoldChange" | !is.null(padj_cutoff)) {
    if (nrow(filter(resultsTable,
                    contrast == contrast_id,
                    padj < padj_cutoff)) == 0) {
      print("No significant DEGs.")
      return(NULL)
    } else if (nrow(filter(resultsTable,
                           contrast == contrast_id,
                           padj < padj_cutoff)) == 1) {
      print("Only one significant DEG.")
      return(NULL)
    }
  }
  
  deg_tbl <- resultsTable %>% filter(contrast == contrast_id)
  if (!is.null(filter_pattern)) {
    deg_tbl <- deg_tbl %>% filter(!grepl(filter_pattern, gene))
  }
  if (direction == "equal") {
    deg_tbl <- deg_tbl %>%
      group_by(-sign(log2FoldChange))
    if (arrange_by == "pvalue") {
      deg_tbl <- deg_tbl %>% arrange(pvalue) %>% slice_head(n = round(slice_n/2)) %>%
        ungroup()
    } else if (arrange_by == "log2FoldChange") {
      deg_tbl <- deg_tbl %>% filter(padj < padj_cutoff) %>%
        arrange(desc(abs(log2FoldChange))) %>%
        slice_head(n = round(slice_n/2)) %>%
        ungroup() %>%
        arrange(desc(log2FoldChange))
    }
  } else if (direction == "up") {
    deg_tbl <- deg_tbl %>% filter(log2FoldChange > 0)
    if (arrange_by == "pvalue") {
      deg_tbl <- deg_tbl %>% arrange(pvalue) %>%
        slice_head(n = slice_n)
    } else if (arrange_by == "log2FoldChange") {
      deg_tbl <- deg_tbl %>% arrange(desc(log2FoldChange)) %>%
        filter(padj < padj_cutoff) %>% slice_head(n = slice_n)
    }
  } else if (direction == "down") {
    deg_tbl <- deg_tbl %>% filter(log2FoldChange < 0)
    if (arrange_by == "pvalue") {
      deg_tbl <- deg_tbl %>% 
        arrange(pvalue) %>% 
        slice_head(n = slice_n)
    } else if (arrange_by == "log2FoldChange") {
      deg_tbl <- deg_tbl %>% 
        arrange(log2FoldChange) %>%
        filter(padj < padj_cutoff) %>% 
        slice_head(n = slice_n)
    }
  } else if (direction == "unequal") {
    if (arrange_by == "pvalue") {
      deg_tbl <- deg_tbl %>% 
        arrange(pvalue) %>% 
        slice_head(n = slice_n)
    } else if (arrange_by == "log2FoldChange") {
      deg_tbl <- deg_tbl %>% 
        arrange(desc(abs(log2FoldChange))) %>%
        filter(padj < padj_cutoff) %>% slice_head(n = slice_n) %>%
        arrange(desc(log2FoldChange))
    }
  } else {
    stop("'direction' parameter must be set to 'unequal' (default), 'equal', 'up', or 'down'")
  }
  return(deg_tbl)
}

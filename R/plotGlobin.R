plotGlobinVln <- function(countTab = bulk$prefilt_dge$counts,
                          globinGenes = globinGenes){
countTab %>%
  as_tibble(rownames = "gene") %>%
  dplyr::mutate(globin = factor(gene %in% globinGenes, 
                                levels = c("TRUE","FALSE"))) %>%
  dplyr::mutate(globin = recode(globin, 
                                "TRUE" = "globin", "FALSE" = "not_globin")) %>%
  pivot_longer(cols = -c(gene, globin), 
               names_to = "sampleID", 
               values_to = "counts") %>%
  mutate(sampleID = fct(sampleID, levels = levels(bulk$dge$samples$sampleID))) %>%
  group_by(sampleID, globin) %>%
  dplyr::summarize(lib.size = sum(counts)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(sampleID),
              names_from = globin, 
              values_from = lib.size) %>%
  mutate(globinPct = globin/(globin+not_globin)) %>% 
  ggplot(aes(x=1,y=globinPct)) +
  geom_violin(fill = NA) + 
  geom_sina() + 
  theme_bw()
}

plotGlobinBars <- function(countTab = bulk$prefilt_dge$counts,
                           globinGenes = globinGenes){
  countTab %>%
    as_tibble(rownames = "gene") %>%
    dplyr::mutate(globin = factor(gene %in% globinGenes, 
                                  levels = c("TRUE","FALSE"))) %>%
    dplyr::mutate(globin = recode(globin, 
                                  "TRUE" = "globin", "FALSE" = "not_globin")) %>%
    pivot_longer(cols = -c(gene, globin), 
                 names_to = "sampleID", 
                 values_to = "counts") %>%
    mutate(sampleID = fct(sampleID, levels = levels(bulk$dge$samples$sampleID))) %>%
    group_by(sampleID, globin) %>%
    dplyr::summarize(lib.size = sum(counts)) %>%
    ungroup() %>%
    ggplot(aes(x = sampleID, y = lib.size, fill = globin)) +
    geom_bar(stat="identity") +
    ylab("Library size (sum of counts)") + xlab("Sample ID") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=6))
}

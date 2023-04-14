library(msigdbr)

m_t2g_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% 
  dplyr::select(gs_name, gene_symbol)
m_t2g_biocarta <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:BIOCARTA") %>% 
  dplyr::select(gs_name, gene_symbol)
m_t2g_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_name, gene_symbol)
m_t2g_h <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)
m_t2g <- bind_rows(m_t2g_reactome,m_t2g_h,m_t2g_biocarta,m_t2g_kegg)
## Script to assign colors for various experimental designs

## 1 - Repeated measures with 6 observations per 6 animals, animals associated with two factors, one with 2 levels, the other with 3 levels

bulk$md$sampleTable

n_factors = 2 # trt, grp
n_levels = c(2, 4)

pal_tbl <- tibble(
  pal = c("Reds", "Oranges", "YlOrBr", "Greens", "Blues", "Purples", "RdPu", "Greys"),
  side = c("warm", "warm", "warm", "cool", "cool", "cool", "cool", "warm"),
  priority = c(1, 5, 3, 4, 2, 6, 8, 7)
) %>%
  arrange(priority) %>%
  rowwise() %>%
  mutate(paired1 = RColorBrewer::brewer.pal(n = 9, name = pal)[4],
         paired2 = RColorBrewer::brewer.pal(n = 9, name = pal)[7]) %>%
  ungroup()


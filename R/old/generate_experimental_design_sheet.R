library(yaml)
library(here)
library(stringr)

config <- yaml.load_file(here("config/QC_config.yml"))
data_dir <- config$alignmentDir
files <- dir(data_dir, recursive = FALSE,
               include.dirs = TRUE, pattern = 'p[0-9]+.*')
samples <- str_extract(files, 's[0-9]+')

exp_des <- data.frame(FileID = files, SampleID = samples, Group = '?')
write.table(exp_des, file = here('config/exp_design.txt'),
            sep = '\t', quote = FALSE,
            col.names = TRUE, row.names = FALSE)

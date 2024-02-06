# generate_experimental_design_sheet.R

Quick script to generate a starter experimental design file. Parses the names of
files in the directory specified in the config to populate a table which you 
can fill in with metadata.

# results_chunk_template.Rmd

Example script of how to use dynamic chunk rendering for lenghty RMD reports.
You can call the script in a function like `knitr::knit_expand()`.

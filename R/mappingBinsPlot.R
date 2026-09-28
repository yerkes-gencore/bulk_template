#' Parse STAR readsPerGene.out files
#'
#' Returns the counts per gene and read count bins for a list of STAR output files
#'
#' @rdname parseReadPerGeneFiles
#' @param file.paths Named vector of file paths to STAR output files, where names
#'    correspond to sample names
#' @param library.type String of column to pull from in the output file. Options
#'  should be 'unstranded', 'sense', or 'antisense'
#'
#' @returns List of 2 data frames, one containing read mapping bins, the other
#'  read per gene counts
#'
#' @examples
#' \dontrun{
#' readfiles <- sapply(
#'   analysis$samplefileIDs,
#'   function(sid) {
#'     paste0(
#'       dir("my_STAR_output_dir",
#'         pattern = sid, full.names = TRUE
#'       ),
#'       "/", sid, "readsPerGene.out.tab"
#'     )
#'   }
#' )
#'
#' outs <- readCountFiles(readfiles, "unstranded")
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom stringr str_to_title
#' @importFrom readr read_tsv
#' @importFrom matrixStats colSums2
#' @importFrom rlang .data
#'
#' @export
parseReadPerGeneFiles <- function(file.paths, library.type = "unstranded") {
  raw_out <- sapply(file.paths,
                    function(x) {
                      readr::read_tsv(x,
                                      col_names = c(
                                        "gene_id",
                                        "unstranded_count",
                                        "sense_count",
                                        "antisense_count"
                                      ),
                                      col_types = list(
                                        "gene_id" = "c",
                                        "unstranded_count" = "d",
                                        "sense_count" = "d",
                                        "antisense_count" = "d"
                                      )
                      ) %>%
                        select(.data$gene_id, starts_with(library.type))
                    },
                    simplify = FALSE,
                    USE.NAMES = TRUE
  )
  
  map_bins <- sapply(raw_out,
                     function(x) {
                       x[c(1:4), ][[paste0(library.type, "_count")]]
                     },
                     USE.NAMES = TRUE
  )
  rownames(map_bins) <- unlist(raw_out[[1]][1:4, 1])
  
  read_counts <- sapply(raw_out,
                        function(x) {
                          x[-c(1:4), ][[paste0(library.type, "_count")]]
                        },
                        USE.NAMES = TRUE
  )
  rownames(read_counts) <- unlist(genes <- raw_out[[1]][-c(1:4), 1])
  
  map_bins <- rbind(map_bins, "N_identified" = matrixStats::colSums2(read_counts))
  
  return(list(map_bins = map_bins, read_counts = read_counts))
}


#' Mapping bins plot
#'
#' Plot values of mapping bins for each sample. Takes the output of `parseReadPerGeneFiles()`
#'
#' @rdname mappingBinsPlot
#' @param mapBins Data frame of 5 mapping bin categories for each sample
#' @param title   Optional title for the plot
#'
#' @examples 
#' \dontrun{
#' mappingBinsPlot(analysis$mapBins)
#' }
#'
#' @returns ggplot object
#'
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom forcats fct_inorder
#' @importFrom stringr str_split str_to_title
#' @importFrom rlang .data
#'
#' @export
mappingBinsPlot <- function(mapBins, title = "") {
  data <- tibble::rownames_to_column(tibble::as_tibble(mapBins, rownames = NA),
                                     var = "map_result"
  ) %>%
    tidyr::pivot_longer(!.data$map_result, names_to = "SampleID", values_to = "count")
  data$map_result <- unlist(lapply(
    data$map_result,
    function(x) stringr::str_to_title(stringr::str_split(x, "_")[[1]][2])
  ))
  data$map_result <- factor(data$map_result,
                            levels = c(
                              "Unmapped", "Multimapping",
                              "Nofeature", "Ambiguous", "Identified"
                            )
  )
  
  ggplot(data, aes(x = .data$SampleID, y = .data$count, fill = .data$map_result)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c(
      "Unmapped" = "grey",
      "Multimapping" = "#DB890E",
      "Nofeature" = "#C60A19",
      "Ambiguous" = "#C4B31C",
      "Identified" = "#209964"
    )) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    guides(fill = guide_legend(title = "Map Result")) +
    aes(fct_inorder(.data$SampleID)) +
    labs(x = "Sample", y = "Number of reads") +
    ggtitle(title)
}
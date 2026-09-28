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
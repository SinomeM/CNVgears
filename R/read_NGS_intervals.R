#' Read genomic intervals
#'
#' \code{read_intervals} handles inputs of data used as the markers in CNVs
#' calling/segmentation using NGS data (WES or WGS)
#'
#' This function is used to load data in interval list or BED like formats into
#' a \code{data.table} that integrates with the other functions of the package.
#' This is usually done at the beginning of a project involving CNVs
#' calling/segmentation on NGS data (WES or WGS) pipelines' results analysis.
#' The function should automatically skip any eventual header.
#' The parameters default values are for file in GATK interval list like format.
#'
#' @param DT_path, path to the input file.
#' @param chr_col, name of the column containing the chromosome information in
#'   the input file.
#' @param start_col, name of the column containing the start information in the
#'   input file.
#' @param end_col, name of the column containing the end information in the input
#'   file.
#'
#' @return a \code{data.table}, will be of \code{Markers} class in future
#'   versions.
#'
#' @export
#'
#' @import data.table
#'
#' @examples
#' read_NGS_intervals(DT_path = system.file("extdata", "markers_WES_example.tsv",
#' package = "CNVgears"), chr_col = "chr", start_col = "start", end_col = "end")

# tested, OK!


read_NGS_intervals <- function(DT_path, chr_col = "CONTIG",
                               start_col = "START", end_col = "END") {
  # check inputs
  if (missing(DT_path))
    stop("Missing parameter(s)!\n")

  # uniformation subfunction
  DT_uniform_internal <- function(DT_in) {
    colnames(DT_in) <- c("chr", "start", "end")
    DT_in <- chr_uniform(DT_in)
    DT_in[, `:=` (start = as.integer(start), end = as.integer(end))]
    setorder(DT_in, chr, start)
    DT_in$P_ID <- 1:nrow(DT_in)
    return(DT_in)
  }

  # read file, excluding eventual header using grep whithin fread()
  columns <- c(chr_col, start_col, end_col)
  DT <- fread(DT_path, skip = chr_col, select = columns)

  # process in DT_uniform()
  DT <- DT_uniform_internal(DT_in = DT)

  return(DT)
}

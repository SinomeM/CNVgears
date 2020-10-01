#' Read Illumina array raw data
#'
#' \code{read_finalreport} handles inputs of data in FinalReport-like format
#'
#' This function is used to load data in FinalReport like format into a
#' \code{data.table} containing the SNPs markers information (i.e. chromosome
#' and position). The function expect a single file where each markers is
#' present one single time. Similar files are often required or produced by the
#' calling algorithm/pipeline, e.g. the PFB file in PennCNV can be used here.
#'
#' @param DT_path, path to the input file.
#' @param mark_ID_col, name of the column containing the SNP ID information in
#'   the input file. Default is \code{"SNP Name"}.
#' @param chr_col, name of the column containing the chromosome information in
#'   the input file. Default is \code{"Chr"}.
#' @param pos_col, name of the column containing the SNPs position information
#'   in the input file. Default is \code{"Position"}.
#'
#' @return a \code{data.table}, will be of \code{Markers} class in future
#'   versions.
#'
#' @export

read_finalreport_snps <- function(DT_path, mark_ID_col = "SNP Name",
                                  chr_col = "Chr", pos_col = "Position") {
  # check inputs
  if (!file.exists(DT_path)) stop("File do not exist, typo(s)?\n")

  # read file and depending on DT_type drop unneeded columns and rows
  columns <- c(mark_ID_col, chr_col, pos_col)

  DT <- fread(DT_path, skip = chr_col, select = columns)

  # standardize columns name and content using DT_uniform()
  DT <- DT_uniform_snps(DT_in = DT)

  return(DT)
}

DT_uniform_snps <- function(DT_in) {
  colnames(DT_in) <- c("SNP_ID", "chr", "pos")
  DT_in[, `:=` (start = as.integer(pos), end = as.integer(pos))][, pos := NULL]
  DT_in <- chr_uniform(DT_in)
  setorder(DT_in, chr, start)
  DT_in$P_ID <- 1:nrow(DT_in)
  return(DT_in)
}

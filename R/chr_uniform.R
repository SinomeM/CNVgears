#' Uniform chromosome notation
#'
#' This is a function for internal use in the package, it handles
#' the standardize process of chromosome notation within the other functions.
#'
#' @param DT_in, a \code{data.table} with a columns named "chr"
#'
#' @return the same \code{data.table} in input with the "chr" uniformed to the
#'   notation "1", "2", ... , "23", "24" for the chromosomes 1:22, X and Y
#'
#' @export
#'
#' @import data.table
#'
#' @examples
#' DT <- data.table::data.table(chr= c("chr1", "chrX", "chr20"))
#' DT <- chr_uniform(DT)
#' DT


chr_uniform <- function(DT_in) {
  if (!"chr" %in% colnames(DT_in))
    stop("No 'chr' columns found!\n")
  if (!is.data.table(DT_in))
    stop("'DT_in' must be a dat.table1\n")

  DT_in[, chr := tolower(gsub(" ", "", chr))]
  # this won't work if the the notation in the column is not coherent, like the
  # results of biomaRt::getBM()
  if (substr(DT_in$chr[1], 1, 3) == "chr")
    DT_in[, chr := substring(chr, 4)]
  # or sub(".+(\\d+|x|y)", "\\1", chr)
  DT_in[chr == "x", chr := "23"][chr == "y", chr := "24"]
  # drop calls not in chrs 1:22, X, Y
  DT_in <- DT_in[chr %in% as.character(1:24), ]

  return(DT_in)
}

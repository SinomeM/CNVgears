#' Rapid genomic locus annotator for CNV calls
#'
#' \code{genomic_locus} add the locus information to a \code{data.table}
#' contaning CNV calls-like data
#'
#' This function takes a \code{data.table} in the format of
#' \code{\link{read_results}} output and add a columns containing the genomic
#' locus information. By default the file containing the CytoBands location for
#' the specified assebmly are downloaded from UCSC website, but also a local
#' object (as \code{data.table} or \code{data.frame}) can be used via the
#' argument \code{bands}.
#'
#' @param DT_in, a \code{data.table}, usually \code{\link{read_results}} output.
#' @param remote_cytobands, logical, indicates whether the function should
#'   download the CytoBands file or use a local object.
#' @param bands, the local object if \code{remote_cytobands} is set to FALSE.
#' @param assembly, character, specify the genomic assembly. Can be either
#'   "hg18", "hg19" or "hg38".
#' @param keep_str_end, logical, specify if intermediate columns (locus_start
#'   and locus_end) must be kept or discarded.
#'
#' @return a \code{data.table} with one or three additional columns containig
#'   genomic locus notation.
#'
#' @export
#'
#' @import data.table

genomic_locus <- function(DT_in, remote_cytobands = T, bands, assembly = "hg19",
                          keep_str_end = T) {
  # check inputs
  if (!remote_cytobands %in% c(T, F))
    stop("Wrong 'remote_cytobands' format!\n")
  if (!assembly %in% c("hg18", "hg19", "hg38"))
    stop("Wrong 'assembly' format!\n")
  if (!is.data.table(DT_in))
    stop("'DT_in' must be a data.table (usually the output of read_results())!\n")
  if (remote_cytobands == F & missing(bands))
    stop("'remote_cytobands' is set to FALSE but no local cytobands object provided!\n")
  # check also colnames?

  # data.table "set" and ":=" functions act by reference, I create a copy to
  # avoid modifying the original object
  DT <- copy(DT_in)
  rm(DT_in)

  # load cytobands
  if  (remote_cytobands == T)
    bands <- fread(paste0("https://hgdownload.cse.ucsc.edu/goldenPath/",
                          assembly, "/database/cytoBand.txt.gz"))
  else {
    cat("Using local cytoBands file!\n")
    if (!is.data.table(bands)) setDT(bands)
  }
  # check colnames and sort
  colnames(bands) <- c("chr", "start", "end", "locus", "buh")
  setorder(bands, chr, start)
  # uniform chromosome notation
  bands <- chr_uniform(bands)

  # match band subfunction
  match_band <- function(chrom, pos) {
    locus <- bands[chr == chrom & start <= pos & end > pos, locus]
    return(locus)
  }

  # apply the function and crate the columns locus_start & locus_end
  DT[, `:=` (locus_start = mapply(match_band, chr, start),
                locus_end = mapply(match_band, chr, end))][
                  locus_start == locus_end, locus := paste0(chr, locus_start)][
                    locus_start != locus_end, locus := paste0(chr, locus_start,
                                                              "-", locus_end)]
  # deleted intermediate columns if needed
  if (keep_str_end == F)
    DT[, `:=` (locus_start = NULL, locus_end = NULL)]

  return(DT)
}

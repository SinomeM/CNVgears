
#' CNVresults to GRanges
#'
#' \code{CNVresults_to_GRanges} convert \code{CNVresults} objects into \code{GRanges}
#'
#' A simple wrapper for the function \code{GenomicRanges::makeGRangesFromDataFrame}.
#' Retained metadata columns are: sample_ID, GT and meth_ID.
#'
#' @param DT a \code{CNVresults} object.
#'
#' @export
#'
#' @import data.table

CNVresults_to_GRanges <- function(DT) {
  GR <- GenomicRanges::makeGRangesFromDataFrame(DT)
  GenomicRanges::mcols(GR) <- DT[, .(sample_ID, GT, meth_ID)]
  return(GR)
}


## #' Title
## #'
## #' @return
## #' @export
## #'
## #' @import data.table
##
## GRanges_to_CNVresults <- function(GR) {
##
##   return(DT)
## }

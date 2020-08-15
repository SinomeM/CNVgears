
#' Title
#'
#' @return
#' @export
#'
#' @import data.table

CNVresults_to_GRanges <- function(DT) {
  GR <- GenomicRanges::makeGRangesFromDataFrame()
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

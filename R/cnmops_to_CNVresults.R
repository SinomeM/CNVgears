


#' Convert cn.mops results into \code{CNVgears} format
#'
#' @param cnvs
#'
#' @import data.table
#'
#' @export


cnmops_to_CNVresults <- function(cnvs, samples_list) {

  res <- GenomicRanges::as.data.frame(cnvs)
  setDT(res)
  setnames(res, c("seqnames", "sampleName"), c("chr", "sample_ID"))
  res[, .(chr, start, end, sample_ID, CN)]

  # compute GT from  CN and sex
  for (s in unique(res$sample_ID)) {
    sex <- samples_list[sample_ID == s, sex]
    # ...
  }

  class(res) <- c("CNVresults", class(res))
  return(res)
}

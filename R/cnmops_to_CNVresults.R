


#' Convert cn.mops results into \code{CNVgears} format
#'
#' @param cnRes cn.mops results (after integer CN calling)
#' @param sample_list minimal cohort metadata, a \code{data.table} produced by the
#'   function \code{\link{read_metadt}}.
#' @param markers a \code{data.table} containing the marker list, the output
#'   \code{\link{read_finalreport_snps}} with \code{DT_type} set to "markers"
#'   or \code{\link{read_NGS_intervals}}.
#'
#' @return the input object \code{cnvs} converted into \code{CNVresults}
#'
#' @import data.table
#'
#' @export
#'
#' @examples
#' library(cn.mops)
#' data(cn.mops)
#' resCNMOPS <- cn.mops(XRanges)
#' resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)
#' resCNMOPS_cnvs <- cnvs(resCNMOPS)
#' # cnmops_calls <- cnmops_to_CNVresults(resCNMOPS_cnvs, sample_list, markers)


cnmops_to_CNVresults <- function(cnRes, sample_list, markers) {

  tmp <- GenomicRanges::as.data.frame(cn.mops::cnvs(cnRes))
  setDT(tmp)
  setnames(tmp, c("seqnames", "sampleName"), c("chr", "sample_ID"))
  tmp <- tmp[, .(chr, start, end, sample_ID, CN)]
  tmp[, CN := as.integer(substr(CN, 3,3))]

  res <- data.table()
  # compute GT from  CN and sex
  for (s in unique(tmp$sample_ID)) {
    sex <- sample_list[sample_ID == s, sex]
    res <- rbind(res, DT_uniform_internal(tmp[sample_ID == s, ], markers, sex))
  }

  class(res) <- c("CNVresults", class(res))
  return(res)
}

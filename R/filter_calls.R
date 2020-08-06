#' Filter CNVs calls based on several parameter
#'
#' \code{cleaning_filter} "clean" a CNV call dataset based on measures such as
#' length number of calls per sample and more.
#'
#' This function can be used together with \code{\link{summary_stats}} in order
#' to clean the dataset from possible noise and unwanted calls. It is generally
#' recommended to briefly explore the data using \code{\link{summary_stats}} and
#' then proceeding to filter out any unwanted group of events. Mandatory
#' arguments of the function are "results" and "min_len"/"min_NP", default
#' values are the authors suggested minimal filtering step, however its quite
#' common to filter anything shorter than 10 or even 50 kb, and/or any call made
#' by less than 10 points. The use of blacklist of any kind is optional and
#' should be done with caution, as it can filter potential biologically relevant
#' events. Over-segmented samples the user wishes to exclude can be specified via
#' \code{blacklist_samples}. Immunoglobulin regions can be generated with the function
#' \code{\link{immuno_regions}}, while telomeric and/or centromeric regions can be obtained with
#' \code{\link{telom_centrom}}.
#'
#' @param results, a \code{data.table}, the output of
#'   \code{\link{read_results}}.
#' @param min_len, minimum CNVs length, any shorter event will be filtered out.
#'   Default is 5000.
#' @param min_NP, minimum CNVs points, any shorter event will be filtered out.
#'   Default is 5.
#' @param blacklists, blacklist, a list of \code{data.table} containing at least
#'   the following columns: "chr", "start", "end". Any event in a blacklists'
#'   region will be filtered out.
#' @param blacklist_samples, character vector containing samples ID to filter
#'   out.
#' @param blacklist_chrs, character vector containing chromosomes names in the
#'   package format, 1:22 for autosomes and 23 24 for chr X and Y.
#'
#' @export

# add also the possibility to filter based on the number of calling algorithms


cleaning_filter <- function(results, min_len = 10000, min_NP = 10,
                            blacklists = NULL, blacklist_samples = NULL,
                            blacklist_chrs = NULL) {
  if (!is.data.table(results))
    stop("Input must be a 'data.table'!\n")
  if (!is.list(blacklists))
    stop("blacklists must be a list object!\n")

  # create a local copy
  DT <- results

  # blacklist_samples and blacklist_chrs
  if (length(blacklist_samples) != 0)
    DT <- DT[!sample_ID %in% blacklist_samples, ]
  if (length(blacklist_chrs) != 0)
    DT <- DT[!chr %in% blacklist_chrs, ]

  # min_len min_NP
  if (!is.na(min_len)) DT <- DT[len >= min_len, ]
  if (!is.na(min_NP)) DT <- DT[NP >= min_NP, ]

  filter_region <- function(DT, bl) {
    bl[, `:=` (start = as.integer(start), end = as.integer(end))]
    # there shouldn't be a lot of regions, un-optimized for loop OK for now
    for (i in 1:nrow(bl)) {
      st <- bl$start[i]
      en <- bl$end[i]
      cc <- bl$chr[i]
      DT <- DT[!(chr == cc & (between(start, st, en) | between(end, st, en) |
                              (start < st & end > en))), ]
    }
    return(DT)
  }

  # user should be careful using these region-based filters, they can remove
  # biologically relevant CNVs

  # blacklists
  if (length(blacklists) == 0) {
    for (i in 1:length(blacklists))
      DT <- filter_region(DT, blacklists[[i]])
  }

  return(DT)
}

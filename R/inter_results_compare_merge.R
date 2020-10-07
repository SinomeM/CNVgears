
#' Combine the results from multiple methods in a single object
#'
#' \code{inter_res_merge} combines the results of multiple CNV calling methods
#' imported via \code{\link{read_results}} into a single \code{CNVresults} object.
#'
#' Multiple \code{CNVresults} objects are combined into a single object and replicated
#' calls are merged. The amount of reciprocal overlap necessary to define two calls
#' as "replicated" is controlled with the \code{prop} parameter.
#' When merging two or more events there will be two borders, inner and outer. By
#' default the outer border are kept as "start" and "end" of the final table. If
#' the user want to keep all the information \code{inner_outer} can be set to NA.
#' In this case "start" and "end" will also be the outer borders. This is done because
#' a \code{CNVresults} object need explicit "start" and "end" columns.
#'
#' @param res_list a \code{list} of \code{CNVresults}.
#' @param sample_list cohort information
#' @param chr_arms a \code{data.table} containing chromosomal arms locations. They
#'   are bundled in the package for hg18, hg19 and hg38 (\code{hgXX_chr_arms}).
#' @param prop proportion of reciprocal overlap to define two calls as "replicated".
#' @param inner_outer keep "inner" or "outer" borders? If NA all columns will be kept.
#'
#' @return a \code{CNVresults} containing the merge of the one provided via
#'   \code{res_list}.
#'
#' @export
#'
#' @examples
#'
#' DT <- inter_res_merge(res_list = list(penn_22, quanti_22), sample_list= cohort_examples, chr_arms= hg19_chr_arms)

# The biggest problem here is the "classic" big CNV detected in method "A"
# and splitted or only partially called in method "B" and the opposite one
# i.e. two separate CNVs detected by method "A" and erroneously merged in
# method "B". From our point of view the two cases are more or less identical
# because in fact we can not determine the actual "true call".
# To solve the situation there are two approaches in my opinion:
# 1. Bigger "wins", meaning that the function should screen the bigger
# calls first and let them incorporate the smaller one. In this way
# several information are lost, in particular the big call is treated
# as replicated even if in fact only a part of it is actually replicated.
# 2. Only CNVs of compatible sizes can be merged.
# In this case there is a second choice, i.e. only merged calls should be
# treated as replicated? This is quite an important question given the fact that
# an heavy filtering step is based on the number of calling methods.

# A combination of the two is possible, as it is possible to let the user
# decide after reading the manual. In fact it also depend on how much the
# user wants to rely on the filter based on number of calling methods and
# whether it is preferred to lower the false positives or the false negatives.

inter_res_merge <- function(res_list, sample_list, chr_arms, prop = 0.3,
                            inner_outer = "outer") {
  if (!is.list(res_list))
    stop("res_list is not a list and it should!\n")
  if (!is.data.table(chr_arms))
    stop("Inputs must be data.table!\n")
  if (!is.na(inner_outer)) {
    if (!inner_outer %in% c("inner", "outer"))
      stop("Invalid 'inner_outer' format!\n")
  }

  # check if the meth_IDs are unique and if the list elements are of the correct
  # class
  mids <- c()
  for (i in 1:length(res_list)) {
    if (!"CNVresults" %in% class(res_list[[i]]))
      stop("All results in the list should be of CNVresults class!\n")
    mids <- c(mids, res_list[[i]][, meth_ID][1]) # ugly but should work (?)
  }
  if (length(mids) != length(unique(mids)))
    stop("Non-unique methods IDs provided in the results!\n")

  # the number of algorithms is used to reduce the number of comparison
  # downstream, in particular only n_algs events on each side of an event will
  # be compared.
  n_alg <- length(mids)*5

  # combine results in a single data.table, keep only CNVs
  res <- data.table()
  for (i in 1:length(res_list)) {
    # compute length if for some reason is not present
    if (!"len" %in% colnames(res_list[[i]])) res_list[[i]][, len := end-start+1]

    res <- rbind(res, res_list[[i]][GT != 0, .(chr, start, end, sample_ID,
                                               GT, CN, meth_ID, len)])
  }
  setorder(res, len, chr, start)

  res_merge <- data.table()
  for (arm in chr_arms$arm_ID) {
    cat("\nMerging calls in chromosomal arm:", arm, "...\n")

    a_chr <- chr_arms[arm_ID == arm, chr]
    a_start <- chr_arms[arm_ID == arm, start]
    a_end <- chr_arms[arm_ID == arm, end]

    if (a_chr %in% 23:24) {
      cat("chr X e Y skipped for now")
      next
    }

    for (samp in sample_list$sample_ID) {

      # subset compatible cnvs
      tmp <- res[sample_ID == samp & (chr == a_chr &
                                      between(start, a_start, a_end) &
                                      between(end, a_start, a_end)), ]
      setorder(tmp, len, chr, start)
      tmp$ix <- 1:nrow(tmp)
      tmp$used <- FALSE

      for (i in tmp$ix) {

        tmp_line <- tmp[ix == i, ]
        # if the line has already been used, skip
        if (tmp_line$used) next

        reg_ref <- get_region(tmp_line, prop)
        tmp_GT <- tmp_line$GT

        merge_ixs <- c(i)

        # questo for puo' diventare una funzione (compare_neighbors())
        for (k in -n_alg:n_alg) {
          if (k == 0) next
          ii <- i + k

          # check if successive calls have desired reciprocal overlap, if GT is correct
          if (nrow(tmp[ix == ii, ]) == 0) next
          if (tmp[ix == ii, used]) next
          if (tmp[ix == ii, GT] != tmp_GT) next
          if (tmp[ix == ii, meth_ID] == tmp_line$meth_ID) next

          reg <- get_region(tmp[ix == ii, ])

          overlap <- min(reg[3], reg_ref[3]) - max(reg[2], reg_ref[2]) + 1
          if (overlap >= reg_ref[4] & overlap >= reg[4] * prop)
            merge_ixs <- c(merge_ixs, ii)
        }
        # now merge_ixs contains the ix values of "mergable" cnv

        # questo passaggio anche puo' diventare una funzione separata, puo' uscire
        # anche semplicemente la data.table da fare rbind con res_merge (esternamente)
        if (length(merge_ixs) == 1) {
          # no hits
          res_merge <-
            rbind(res_merge,
                  data.table("chr" = a_chr, "inner_start" = reg_ref[2],
                             "inner_end" = reg_ref[3], "outer_start" = reg_ref[2],
                             "outer_end" = reg_ref[3], "sample_ID" = samp,
                             "GT" = tmp_GT, "CN" = tmp_line$CN,
                             "meth_ID" = tmp_line$meth_ID, "n_meth" = 1))
        }
        else {
          # do merge
          sel_lines <- tmp[ix %in% merge_ixs, ]
          merged_line <- do_merge(sel_lines)
          meths <- paste(sort(sel_lines$meth_ID), collapse = "-")
          res_merge <-
            rbind(res_merge,
                  data.table("chr" = merged_line[1], "inner_start" = merged_line[2],
                             "inner_end" = merged_line[3], "outer_start" = merged_line[4],
                             "outer_end" = merged_line[5], "sample_ID" = samp,
                             "GT" = tmp_GT, "CN" = merged_line[6], "meth_ID" = meths,
                             "n_meth" = length(sel_lines$meth_ID)))
        }
        # mark "merged" lines, in this way an event won't be used more than once time
        tmp[ix %in% merge_ixs, used := TRUE]
      }
    }
  }
  # re-create seg_ID, unique for each sample
  setorder(res_merge, chr, outer_end, outer_start)
  res_merge[, seg_ID := 1:.N, by = sample_ID]
  # re-create start end
  res_merge <- start_end(res_merge, inner_outer)
  # re-set the class
  class(res_merge) <- c("CNVresults", class(res_merge))
  return(res_merge)
}

## Sub-Functions

compare_neighbors <- function() {

}


do_merge <- function(my_lines) {
  chr <- my_lines$chr[1]
  inner_st  <- max(my_lines$start)
  outer_st  <- min(my_lines$start)
  inner_en  <- min(my_lines$end)
  outer_en  <- max(my_lines$end)
  CN <- round(mean(my_lines$CN))
  res <- c(chr, inner_st, inner_en, outer_st, outer_en, CN)
  return(as.integer(res))
}

start_end <- function(DT, in_out) {
  if (!is.na(in_out)) {
    if (in_out == "outer") {
      setnames(DT, c("outer_start", "outer_end"), c("start", "end"))
      DT[, `:=` (inner_start = NULL, inner_end = NULL)]
    }
    if (in_out == "inner") {
      setnames(DT, c("inner_start", "inner_end"), c("start", "end"))
      DT[, `:=` (outer_start = NULL, outer_end = NULL)]
    }
    DT[, len := end - start + 1]
  }
  else
    DT[, `:=` (start = outer_start, end := outer_start, len = end - start + 1)]

  return(DT)
}

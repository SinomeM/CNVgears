
#' Combine the results from multiple methods in a single object
#'
#' @param res_list aaa
#' @param sample_list aaa
#' @param g_arms aaa
#'
#' @return
#' @export
#'

# Best option should be to load all the results in a list, rbind together
# all the elements of the list and then, iterating per each chromosomes
# arm (per sample) do something pretty similar to the CNVRs construction
# given that the problem is more or less exactly the same

# here the markers objects are not relevant anymore because we could be
# combining results from different data types.

# The biggest problem here is the "classic" big CNV detected in method "A"
# and splitted or only partially called in method "B" and the opposite one
# i.e. two separate CNVs detected by method "A" and erroneously merged in
# method "B". From our point of view the two cases are more or less identical
# because in fact we can not determine the actual "true call".
# To solve the situation there are two approaches in my opinion:
# 1. Bigger "wins", meaning that the function should screen the bigger
# calls first and let them incorporate the smaller one. In this way
# several informations are lost, in particular the big call is treated
# as replicated even if in fact only a part of it is actually replicated.
# 2. Only CNVs of compatible sizes can be merged.
# In this case there is a second choice, i.e. only merged calls should be
# treated as replicated? This is quite an important question given the fact that
# an heavy filtering step is based on the number of calling methods.

# A combination of the two is possible, as it is possible to let the user
# decide after reading the manual. In fact it also depend on how much the
# user wants to rely on the filter based on number of calling methods and
# whether it is preferred to lower the false positives or the false negatives.

inter_res_merge <- function(res_list, sample_list, g_arms, prop = 0.3,
                            inner_outer = "outer") {
  if (!is.list(res_list))
    stop("res_list is not a list and it should!\n")
  if (!is.data.table(g_arms))
    stop("Inputs must be data.table!\n")
  if (!is.na(inner_outer)) {
    if (!inner_outer %in% c("inner", "outer"))
      stop("Invalid 'inner_outer' format!\n")
  }

  # check if the meth_IDs are unique and if the list elements are of the correct
  # class
  mids <- c()
  for (i in 1:length(res_list)) {
    if (!is.data.table(res_list[[i]]))
      stop("All results in the list should be data.table!\n")
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
    # compute lenght if for some reson is not present
    if (!"len" %in% colnames(res_list[[i]])) res_list[[i]][, len := end-start+1]

    res <- rbind(res, res_list[[i]][GT != 0, .(chr, start, end, sample_ID,
                                               GT, CN, meth_ID, len)])
  }
  setorder(res, len, chr, start)

  res_merge <- data.table()
  for (arm in g_arms$arm_ID) {
    cat("\nMerging calls in chromosomal arm:", arm, "...\n")

    a_chr <- g_arms[arm_ID == arm, chr]
    a_start <- g_arms[arm_ID == arm, start]
    a_end <- g_arms[arm_ID == arm, end]

    if (a_chr %in% 23:24) {
      cat("chr X e Y skipped for now")
      next
    }

    for (samp in sample_list$sample_ID) {
      # cat("Sample :", samp, "...\n")

      # subset compatible cnvs
      tmp <- res[sample_ID == samp & (chr == a_chr &
                                      between(start, a_start, a_end) &
                                      between(end, a_start, a_end)), ]
      setorder(tmp, len, chr, start)
      tmp$ix <- 1:nrow(tmp)

      for (i in tmp$ix) {
        tmp_line <- tmp[ix == i, ]
        # if the line has already been used, skip
        if (nrow(tmp_line) == 0)
          next

        reg_ref <- get_region(tmp_line, prop)
        tmp_GT <- tmp_line$GT

        merge_ixs <- c(i)

        # questo for puo' diventare una funzione (compare_neighbors())
        for (k in -n_alg:n_alg) {
          if (k == 0) next

          # check if successive calls have desired reciprocal overlap, if GT is correct
          ii <- i + k
          if (nrow(tmp[ix == ii, ]) == 0) next
          if (tmp[ix == ii, GT] != tmp_GT) next
          if (tmp[ix == ii, meth_ID] == tmp_line$meth_ID) next

          reg <- get_region(tmp[ix == ii, ])

          overlap <- min(reg[3], reg_ref[3]) - max(reg[2], reg_ref[2]) + 1
          if (overlap >= reg_ref[4] & overlap >= reg[4] * prop)
            merge_ixs <- c(merge_ixs, ii) # bad thing grow a vector but it should not be longer than 10 in any scenario
        }
        # now merge_ixs contains the ix values of mnergable cnv, must create a unique line, add it
        # to res_merge and remove the original calls from tmp
        # the output object should have the following columns:
        # "chr", "inner_start", "inner_end", "outer_start", "outer_end", "sample_ID", "GT", ("CN"), "seg_ID", "meth_ID", "n_meth"
        # ... #
        # CAREFUL HERE, not necessarily length(merge_ixs) > 1 !!!

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
          sel_lines  <- tmp[ix %in% merge_ixs, ]
          # GT_m  <- tmp_GT # un-necessary
          inner_st  <- max(sel_lines$start)
          outer_st  <- min(sel_lines$start)
          inner_en  <- min(sel_lines$end)
          outer_en  <- max(sel_lines$end)
          CN_m <- round(mean(sel_lines$CN))
          mm <- sort(sel_lines$meth_ID)
          meth_ID_m <- paste(mm, collapse = "-")

          res_merge <- rbind(res_merge, data.table("chr" = a_chr,
                                                   "inner_start" = inner_st, "inner_end" = inner_en,
                                                   "outer_start" = outer_st, "outer_end" = outer_en,
                                                   "sample_ID" = samp, "GT" = tmp_GT, "CN" = CN_m,
                                                   "meth_ID" = meth_ID_m, "n_meth" = length(mm)))
        }

        # remove the "merged" lines, in this way an event won't be used more than
        # one time
        # invece di aggiornare tmp, aggiungi la colonna "used" quando crei la colonna "ix"
        # tutti FALSI e cambia in VERO man mano.
        tmp <- tmp[!ix %in% merge_ixs, ]
      }
    }
  }

  # re-create seg_ID, unique for each sample
  setorder(res_merge, chr, outer_end, outer_start)
  res_merge[, seg_ID := 1:.N, by = sample_ID]

  if (!is.na(inner_outer)) {
    if (inner_outer == "outer") {
      setnames(res_merge, c("outer_start", "outer_end"), c("start", "end"))
      res_merge[, `:=` (inner_start = NULL, inner_end = NULL)]
    }
    if (inner_outer == "inner") {
      setnames(res_merge, c("inner_start", "inner_end"), c("start", "end"))
      res_merge[, `:=` (outer_start = NULL, outer_end = NULL)]
    }
    res_merge[, len := end - start + 1]
  }

  # re set the class
  class(res_merge) <- c("CNVresults", class(res_merge))
  return(res_merge)
}


compare_neighbors <- function() {

}

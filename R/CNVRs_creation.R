#' Compute Copy Number Variable Regions (CNVRs)
#'
#' @param cnv, a data.table like the one produced by\code{\link{read_results}},
#'  typically the results of \code{\link{inter_res_merge}}.
#' @param g_arms, a data.table containing the genomic location of the genomic
#'  arms. For the assemblies hg38 and hg19 it is provided by the package.
#' @param prop, reciprocal overlap proportion, default 0.3 (30%).
#' @param inner_outer, specify whether inner or outer start/end should be used. If
#' NA the function assume the "start" "end" columns are present.
#'
#' lorem ipsum etc etc
#'
#' @return
#'
#' @export
#'
#' @import data.table
#'

# OK for now, it would be nice to add a third step where for all CNVs in a CNVRs
# it check the adjacent CNVRs if there is a better overlap, not super easy

cnvrs_create <- function(cnvs, g_arms, prop = 0.3) {
  # check input
  if (!is.data.table(cnvs) | !is.data.table(g_arms))
    stop("Inputs must be data.table!\n")

  # data.table "set" and ":=" functions act by reference, I create a copy to
  # avoid modifying the original object (perhaps there is a better option?)
  cnvs_cp <- cnvs
  rm(cnvs)

  # check input formats, in particular for g_arms
  g_arms[, `:=` (start = as.integer(start), end = as.integer(end))]
  g_arms <- chr_uniform(g_arms)

  # sort per chr, start & end
  setorder(cnvs_cp,	chr, start, end)
  setorder(g_arms, chr, start)
  cnvs_cp[, cnvr := NA_character_]
  # this create a line with all "NA", I remove it later
  cnvrs <- data.table("r_ID" = NA_character_, "chr" = NA_character_,
                      "start" = NA_integer_, "end" = NA_integer_)
  # cnvrs: arm_ID, chr, start, end
  res <- data.table()

  for (arm in g_arms$arm_ID) {
    cat("\n-- Computing CNVRs in chromosome arm:", arm, "--\n")
    # cnvrs for this arm
    cnvrs_tmp <- data.table("r_ID" = NA_character_, "chr" = NA_character_,
                            "start" = NA_integer_, "end" = NA_integer_)
    # arm related variables
    k <- T
    a_chr <- g_arms[arm_ID == arm, chr]
    a_start <- g_arms[arm_ID == arm, start]
    a_end <- g_arms[arm_ID == arm, end]
    # subset compatible cnvs
    DT <- cnvs_cp[chr == a_chr &
                    between(start, a_start, a_end) &
                    between(end, a_start, a_end), ]
    DT[, ix := 1:nrow(DT)]
    # if no CNVs are present in the arm, just skip
    if (nrow(DT) == 0) next
    # for each arm keep running until no more cnvs are excluded
    # need to check how DT and the index "i" are created, only the cnvs
    # with cnvr == NA need to be processed
    n_l <- 1
    n <- 1
    while (k == T) {
      # create CNVRs/fill CNVRs
      cat("Creating/filling CNVRs, loop #", n_l, "...\n")
      indexes <- DT[is.na(cnvr), ix]
      cat(length(indexes), " CNVs to be assigned\n")
      for (i in indexes) {
        st <- DT$start[i]
        en <- DT$end[i]
        len <- en - st + 1
        # possible cnvrs
        if (nrow(cnvrs_tmp[!is.na(r_ID), ]) != 0) {
          cnvr_m <- cnvrs_tmp[!is.na(r_ID),][
            between(start, st, en, incbounds = TRUE) |
            between(end, st, en, incbounds = TRUE), ]
        } # this if else is needed for the initial CNVR (of each chr), no better solution ATM
        else cnvr_m <- cnvrs_tmp[!is.na(r_ID), ]
        # no match, initialize new cnvr
        if (nrow(cnvr_m) == 0) {
          cnvrs_tmp <- rbind(cnvrs_tmp, data.table("r_ID" = paste0(arm, "-", n),
                                           "chr" = a_chr, "start" = st,
                                           "end" = en))
          DT$cnvr[i] <- paste0(arm, "-", n)
          n <- n + 1
        }
        # at least one cnvr overlap with cnv under examination
        else {
          for (r in cnvr_m$r_ID) {
            cnvs_tmp <- DT[cnvr == r, ]
            overlaps <- pmin(en, cnvs_tmp$end) - pmax(st, cnvs_tmp$start) + 1
            lens <- cnvs_tmp$end - cnvs_tmp$start + 1
            # reciprocal overlaps
            if (all(overlaps >= len * prop) &
                all(sapply(overlaps, function(x) x >= lens * prop))) {
              DT$cnvr[i] <- r
              # update cnvrs boundaries
              cnvrs_tmp[r_ID == r, `:=` (start = min(start, st), end = max(end, en))]
              break
            }
          }
          # cnvr not set mean no match with candidates cnvrs
          if (is.na(DT$cnvr[i])) {
            cnvrs_tmp <- rbind(cnvrs_tmp, data.table("r_ID" = paste0(arm, "-", n),
                                             "chr" = a_chr, "start" = st,
                                             "end" = en))
            DT$cnvr[i] <- paste0(arm, "-", n)
            n <- n + 1
          }
        }
      }

      cnvrs_tmp <- cnvrs_tmp[!is.na(r_ID) & r_ID %in% DT$cnvr, ]

      # Check if some CNVRs can be merged
      cat("Re-checking CNVRs ...\n")
      if (nrow(cnvrs_tmp > 1)) {
        # this loop should run until no more CNVRs can be merged
        wh <- 1
        while (TRUE == TRUE) {
          if (wh != 1) break
          wh <- 0

          for (i in 1:nrow(cnvrs_tmp)) {
            b <- FALSE
            st <- cnvrs_tmp$start[i]
            en <- cnvrs_tmp$end[i]
            id <- cnvrs_tmp$r_ID[i]
            # search compatible cnvrs
            cnvrs_m <- cnvrs_tmp[between(start, st, en, incbounds = FALSE) |
                                   between(end, st, en, incbounds = FALSE), ]
            # if there are no hits skip
            if (nrow(cnvrs_m) == 0 ) next
            for (mm in 1:nrow(cnvrs_m)) {
              sst <- cnvrs_m$start[mm]
              een <- cnvrs_m$end[mm]
              iid <- cnvrs_m$r_ID[mm]
              overl <- min(en, een) - max(st, sst)
              lens <- c(en-st+1, een-sst+1)
              bb <- TRUE
              # if there is a match, check the relative cnvs
              if (all(overl > lens * prop)) {
                cnvs_tmp <- DT[cnvr %in% c(id, iid), ]
                # test the CNVs
                lens <- cnvs_tmp$end - cnvs_tmp$start + 1
                for (ii in 1:nrow(cnvs_tmp)) {
                  overlaps <- pmin(cnvs_tmp$end[ii], cnvs_tmp$end) -
                              pmax(cnvs_tmp$start[ii], cnvs_tmp$start) + 1
                  # reciprocal overlaps
                  if (!(all(overlaps >= (cnvs_tmp$end[ii]-cnvs_tmp$start[ii]+1) * prop) &
                      all(sapply(overlaps, function(x) x >= lens * prop)))) {
                    bb <- FALSE
                    break
                  }
                }
                if (bb) {
                  # if all cnvs are compatible, update CNVR, remove the old ones
                  # and update DT
                  cnvrs_tmp <- rbind(cnvrs_tmp,
                                       data.table("r_ID" = paste0(arm, "-", n),
                                                  "chr" = a_chr, "start" = min(st, sst),
                                                  "end" = max(en, een)))
                  cnvrs_tmp <- cnvrs_tmp[!r_ID %in% c(id, iid)]
                  DT[cnvr %in% c(id, iid), cnvr := paste0(arm, "-", n)]
                  n <- n + 1
                  wh <- 1
                  cat("cnvr updated \n")
                  # at this point all the for loops are interrupted and another
                  # while loop starts
                  b <- TRUE
                  break
                }
              }
              if (b) break
            }
            if (b) break
          }
        }
      }

      # Recheck CNVs
      cat("Re-checking CNVs ...\n")

      # # Move CVNs if a better overlap is found
      # for (i in DT[, ix]) {
      #
      # }

      # Remove CNVs if necessary
      for (i in DT[, ix]) {
        st <- DT$start[i]
        en <- DT$end[i]
        len <- en - st + 1
        # all cnvs in the cnvr
        cnvs_tmp <- DT[cnvr == DT$cnvr[i]]
        # this is copy-pasted from above, check twice
        overlaps <- pmin(en, cnvs_tmp$end) - pmax(st, cnvs_tmp$start) + 1
        lens <- cnvs_tmp$end - cnvs_tmp$start + 1
        # note that here there is a "!"
        if (!(all(overlaps >= len * prop) &
              all(sapply(overlaps, function(x) x >= lens * prop))))
          DT$cnvr[i] <- NA_character_
      }
      cat(length(is.na(DT$cnvr)[is.na(DT$cnvr) == T]),
          "CNVs removed from the assigned CNVR\n")
      # no CNVRs set to NA mean that no CNV has been excluded in the last iteration
      if (all(!is.na(DT$cnvr))) k <- F
      n_l <- n_l + 1
    }

    # Recreate the output
    res <- rbind(res, DT[, ix := NULL])
    cnvrs <- rbind(cnvrs, cnvrs_tmp)
  }

  # cnvrs can no longer have CNVs, clean those ones
  cnvrs <- cnvrs[r_ID %in% unique(res$cnvr), ]
  # count the final CNVs frequency per CNVRs
  freqs <- res[, .N, by = cnvr]

  cnvrs[, freq := freqs[match(cnvrs$r_ID, freqs$cnvr), N]]

  return(list(cnvrs[!is.na(chr), ], res))
}

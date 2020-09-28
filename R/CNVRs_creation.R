#' Compute Copy Number Variable Regions (CNVRs)
#'
#' \code{cnvrs_create} compute CNVRs on a \code{CNVResutls} object, typically the
#' output of \code{\link{inter_res_merge}}.
#'
#' Copy Number Variable Regions (CNVRs) are defined as groups of reciprocal overlapping
#' CNVs. This function first try to assign every call to a CNVR (or create a new
#' one if it is not possible), then check if adjacent CNVRs cam be merged, and
#' finally recheck all the CNVs in each CNVRs and un-assign them if the reciprocal
#' overlap is no longer fulfilled with all the member of the CNVR. If any event
#' is touched by this last step, a new cycle begins, this until no more CNVs can
#' be removed from the assigned CNVR.
#'
#' @param cnvs a \code{CNVresults} produced by \code{\link{read_results}}.
#' @param chr_arms a \code{data.table} containing chromosomal arms locations. They
#'   are bundled in the package for hg18, hg19 and hg38 (\code{hgXX_chr_arms}).
#' @param prop reciprocal overlap proportion, default 0.3 (30\%).
#'
#' @return a \code{list} of two elements. The first element is a \code{data.table}
#'   that contains the actual CNVR information, genomic location and frequency in
#'   the cohort. The second element is the \code{CNVresults}
#'
#'
#' @export
#'
#' @import data.table
#'

# OK for now, it would be nice to add a third step where for all CNVs in a CNVRs
# it check the adjacent CNVRs if there is a better overlap, not super easy

cnvrs_create <- function(cnvs, chr_arms, prop = 0.3) {
  # check input
  if (!is.data.table(cnvs) | !is.data.table(chr_arms))
    stop("Inputs must be data.table!\n")

  # data.table "set" and ":=" functions act by reference, I create a copy to
  # avoid modifying the original object (perhaps there is a better option?)
  cnvs_cp <- cnvs
  rm(cnvs)

  # check input formats, in particular for chr_arms
  chr_arms[, `:=` (start = as.integer(start), end = as.integer(end))]
  chr_arms <- chr_uniform(chr_arms)

  # sort per chr, start & end
  setorder(cnvs_cp,	chr, start, end)
  setorder(chr_arms, chr, start)
  cnvs_cp[, cnvr := NA_character_]

  # cnvrs: arm_ID, chr, start, end
  # this create a line with all "NA", I remove it later
  cnvrs <- data.table("r_ID" = NA_character_, "chr" = NA_character_,
                      "start" = NA_integer_, "end" = NA_integer_)
  res <- data.table()

  # proceed per chromosomal arms
  for (arm in chr_arms$arm_ID) {
    cat("\n-- Computing CNVRs in chromosome arm:", arm, "--\n")

    # cnvrs for this arm
    cnvrs_tmp <- data.table("r_ID" = NA_character_, "chr" = NA_character_,
                            "start" = NA_integer_, "end" = NA_integer_)
    # arm related variables
    reg_arm <- get_region(chr_arms[arm_ID == arm,])

    # subset compatible cnvs
    DT <- cnvs_cp[chr == reg_arm[1] &
                    between(start, reg_arm[2], reg_arm[3]) &
                    between(end, reg_arm[2], reg_arm[3]), ]
    DT[, ix := 1:nrow(DT)]
    # if no CNVs are present in the arm, just skip
    if (nrow(DT) == 0) next

    # keep track of CNVR and loop number
    n_loop <- 1
    n <- 1

    # for each arm keep running until no more cnv is excluded
    while (any(is.na(DT$cnvr))) {

      # create CNVRs/fill CNVRs
      cat("Creating/filling CNVRs, loop #", n_loop, "...\n")
      indexes <- DT[is.na(cnvr), ix]
      cat(length(indexes), " CNVs to be assigned\n")

      tmp <- create_fill_CNVR(cnvrs_tmp, DT, n, prop, indexes, arm, reg_arm)
      DT <- tmp[[2]]
      cnvrs_tmp <- tmp[[1]][!is.na(r_ID) & r_ID %in% DT$cnvr, ]
      n <- tmp[[3]]

      # Check if some CNVRs can be merged
      cat("Re-checking CNVRs ...\n")
      if (nrow(cnvrs_tmp > 1)) {
        # first run
        tmp <- check_cnvrs(cnvrs_tmp, DT, n, prop)

        # if CNVRs are updated tmp[[4]] is TRUE
        while(tmp[[4]]) {
          tmp <- check_cnvrs(cnvrs_tmp, DT, prop)
        }

        cnvrs_tmp <- tmp[[1]]
        DT <- tmp[[2]]
        n <- tmp[[3]]
      }

      # Recheck CNVs
      cat("Re-checking CNVs ...\n")

      # # Move CVNs if a better overlap is found
      # for (i in DT[, ix]) {
      #
      # }

      # Remove CNVs if necessary
      DT <- remove_cnvs(DT, prop)

      n_loop <- n_loop + 1
    }

    # Recreate the output
    res <- rbind(res, DT[, ix := NULL])
    cnvrs <- rbind(cnvrs, cnvrs_tmp)
  }

  # CNVRs can no longer have CNVs, clean those ones
  cnvrs <- cnvrs[r_ID %in% unique(res$cnvr), ]
  # count the final CNVs frequency per CNVRs
  freqs <- res[, .N, by = cnvr]
  cnvrs[, freq := freqs[match(cnvrs$r_ID, freqs$cnvr), N]]

  return(list(cnvrs[!is.na(chr), ], res))
}


## SUBFUNCTIONS

create_fill_CNVR <- function(cnvrs, DT, n, prop, ixs, arm, reg_arm) {

  for (i in ixs) {
    my_reg <- get_region(DT[ix == i, ], prop)
    # possible cnvrs
    if (nrow(cnvrs[!is.na(r_ID), ]) != 0) {
      cnvr_m <- cnvrs[!is.na(r_ID),][
        between(start, my_reg[2], my_reg[3], incbounds = TRUE) |
          between(end, my_reg[2], my_reg[3], incbounds = TRUE), ]
    } # this if else is needed for the initial CNVR (of each chr)
    else cnvr_m <- cnvrs[!is.na(r_ID), ]

    # no match, initialize new cnvr
    if (nrow(cnvr_m) == 0) {
      cnvrs <- rbind(cnvrs, data.table("r_ID" = paste0(arm, "-", n),
                                               "chr" = reg_arm[1], "start" = my_reg[2],
                                               "end" = my_reg[3]))
      DT$cnvr[i] <- paste0(arm, "-", n)
      n <- n + 1
    }
    # at least one cnvr overlap with cnv under examination
    else {
      for (r in cnvr_m$r_ID) {
        reg_list <- get_regions_list(DT[cnvr == r, ], prop)
        overlaps <- pmin(my_reg[3], reg_list[[3]]) - pmax(my_reg[2], reg_list[[2]]) + 1
        # reciprocal overlaps
        if (all(overlaps >= my_reg[4]) &
            all(sapply(overlaps, function(x) x >= reg_list[[4]]))) {
          DT$cnvr[i] <- r
          # update cnvrs boundaries
          cnvrs[r_ID == r, `:=` (start = min(start, my_reg[2]), end = max(end, my_reg[3]))]
          break
        }
      }

      # cnvr not set mean no match with candidates cnvrs
      if (is.na(DT$cnvr[i])) {
        cnvrs <- rbind(cnvrs, data.table("r_ID" = paste0(arm, "-", n),
                                                 "chr" = reg_arm[1], "start" = my_reg[2],
                                                 "end" = my_reg[3]))
        DT$cnvr[i] <- paste0(arm, "-", n)
        n <- n + 1
      }
    }
  }

  return(list(cnvrs, DT, n))
}

check_cnvrs <- function(cnvrs, DT, n, prop) {

  for (i in 1:nrow(cnvrs)) {
    b <- FALSE
    my_reg <- get_region_with_rID(cnvrs[i])
    # search compatible cnvrs
    cnvrs_m <-
      cnvrs[between(start, my_reg[[1]][2], my_reg[[1]][3], incbounds = FALSE) |
              between(end, my_reg[[1]][2], my_reg[[1]][3], incbounds = FALSE), ]

    # if there are no hits skip
    if (nrow(cnvrs_m) == 0 ) next

    # if there is a match, check the relative cnvs
    for (ii in 1:nrow(cnvrs_m)) {

      my_reg_m <- get_region_with_rID(cnvrs_m[ii], prop)
      overl <-
        min(my_reg[[1]][3], my_reg_m[[1]][3]) - max(my_reg[[1]][2], my_reg_m[[1]][2])
      pass_check <- TRUE

      if (overl > my_reg[[1]][4] & overl > my_reg_m[[1]][4]) {

        reg_list <-
          get_regions_list(DT[cnvr %in% c(my_reg[[2]], my_reg_m[[2]]), ], prop)
        cnvs_tmp <- DT[cnvr %in% c(my_reg[[2]], my_reg_m[[2]]), ]

        # test the CNVs one against all the other
        for (iii in 1:nrow(cnvs_tmp)) {
          overlaps <- pmin(cnvs_tmp$end[iii], cnvs_tmp$end) -
            pmax(cnvs_tmp$start[iii], cnvs_tmp$start) + 1

          # reciprocal overlaps
          if (!(all(overlaps >= (cnvs_tmp$end[iii]-cnvs_tmp$start[iii]+1) * prop) &
                all(sapply(overlaps, function(x) x >= reg_list[[4]])))) {
            pass_check <- FALSE
            break
          }
        }
        if (pass_check) {
          # if all cnvs are compatible, update CNVR, remove the old ones and update DT
          cnvrs <- rbind(cnvrs,
                         data.table("r_ID" = paste0(arm, "-", n), "chr" = reg_arm[1],
                                    "start" = min(my_reg[[1]][2], my_reg_m[[1]][2]),
                                    "end" = max(my_reg[[1]][3], my_reg_m[[1]][3])))
          cnvrs <- cnvrs[!r_ID %in% c(my_reg[[2]], my_reg_m[[2]])]
          DT[cnvr %in% c(my_reg[[2]], my_reg_m[[2]]), cnvr := paste0(arm, "-", n)]
          n <- n + 1
          cat("CNVR updated \n")
          # at this point all the for loops are interrupted
          b <- TRUE
          break
        }
      } # fi
      if (b) break
    }
    if (b) break
  }

  return(list(cnvrs, DT, n, b))
}

move_cnvs <- function() {

}

remove_cnvs <- function(DT, prop) {

  for (i in DT[, ix]) {
    my_reg <- get_region(DT[ix == i, ], prop)
    # all cnvs in the cnvr
    reg_list <- get_regions_list(DT[cnvr == DT$cnvr[i]], prop)

    overlaps <- pmin(my_reg[3], reg_list[[3]]) - pmax(my_reg[2], reg_list[[2]]) + 1

    if (!(all(overlaps >= my_reg[4]) &
          all(sapply(overlaps, function(x) x >= reg_list[[4]]))))
      DT$cnvr[i] <- NA_character_
  }

  cat(length(is.na(DT$cnvr)[is.na(DT$cnvr) == T]),
      "CNVs removed from the assigned CNVR\n")

  return(DT)
}

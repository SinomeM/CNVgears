#' Compute Copy Number Variable Regions (CNVRs)
#'
#' @param cnv, a data.table like the one produced by\code{\link{read_results}},
#'  typically the results of \code{\link{inter_res_merge}}.
#' @param g_arms, a data.table containing the genomic location of the genomic
#'  arms. For the assemblies hg38 and hg19 it is provided by the package.
#' @param prop, reciprocal overlap proportion, default 0.3 (30\%).
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
    reg_arm <- get_region(g_arms[arm_ID == arm,])
    # subset compatible cnvs
    DT <- cnvs_cp[chr == reg_arm[1] &
                    between(start, reg_arm[2], reg_arm[3]) &
                    between(end, reg_arm[2], reg_arm[3]), ]
    DT[, ix := 1:nrow(DT)]
    # if no CNVs are present in the arm, just skip
    if (nrow(DT) == 0) next

    # keep track of CNVR and loop number
    n_l <- 1
    n <- 1
    # for each arm keep running until no more cnv is excluded
    while (any(is.na(DT$cnvr))) {
      # create CNVRs/fill CNVRs
      cat("Creating/filling CNVRs, loop #", n_l, "...\n")
      indexes <- DT[is.na(cnvr), ix]
      cat(length(indexes), " CNVs to be assigned\n")

      for (i in indexes) {
        my_reg <- get_region(DT[ix == i, ], prop)
        # possible cnvrs
        if (nrow(cnvrs_tmp[!is.na(r_ID), ]) != 0) {
          cnvr_m <- cnvrs_tmp[!is.na(r_ID),][
            between(start, my_reg[2], my_reg[3], incbounds = TRUE) |
            between(end, my_reg[2], my_reg[3], incbounds = TRUE), ]
        } # this if else is needed for the initial CNVR (of each chr), no better solution ATM
        else cnvr_m <- cnvrs_tmp[!is.na(r_ID), ]
        # no match, initialize new cnvr
        if (nrow(cnvr_m) == 0) {
          cnvrs_tmp <- rbind(cnvrs_tmp, data.table("r_ID" = paste0(arm, "-", n),
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
              cnvrs_tmp[r_ID == r, `:=` (start = min(start, my_reg[2]), end = max(end, my_reg[3]))]
              break
            }
          }
          # cnvr not set mean no match with candidates cnvrs
          if (is.na(DT$cnvr[i])) {
            cnvrs_tmp <- rbind(cnvrs_tmp, data.table("r_ID" = paste0(arm, "-", n),
                                             "chr" = reg_arm[1], "start" = my_reg[2],
                                             "end" = my_reg[3]))
            DT$cnvr[i] <- paste0(arm, "-", n)
            n <- n + 1
          }
        }
      }

      cnvrs_tmp <- cnvrs_tmp[!is.na(r_ID) & r_ID %in% DT$cnvr, ]

      # Check if some CNVRs can be merged
      ## PROBABLY SOME PROBLEMS HERE
      cat("Re-checking CNVRs ...\n")
      if (nrow(cnvrs_tmp > 1)) {
        # this loop should run until no more CNVRs can be merged
        wh <- 1
        while (TRUE == TRUE) {
          if (wh != 1) break
          wh <- 0

          for (i in 1:nrow(cnvrs_tmp)) {
            b <- FALSE
            my_reg <- get_region_with_rID(cnvrs_tmp[i])
            # search compatible cnvrs
            cnvrs_m <-
              cnvrs_tmp[between(start, my_reg[[1]][2], my_reg[[1]][3],
                                incbounds = FALSE) |
                          between(end, my_reg[[1]][2], my_reg[[1]][3],
                                  incbounds = FALSE), ]

            # if there are no hits skip
            if (nrow(cnvrs_m) == 0 ) next

            # if there is a match, check the relative cnvs
            for (mm in 1:nrow(cnvrs_m)) {
              my_reg_m <- get_region_with_rID(cnvrs_m[mm], prop)
              overl <- min(my_reg[[1]][3], my_reg_m[[1]][3]) - max(my_reg[[1]][2],
                                                                   my_reg_m[[1]][2])
              pass_check <- TRUE
              if (overl > my_reg[[1]][4] & overl > my_reg_m[[1]][4]) {
                reg_list <-
                  get_regions_list(DT[cnvr %in% c(my_reg[[2]], my_reg_m[[2]]), ], prop)
                cnvs_tmp <- DT[cnvr %in% c(my_reg[[2]], my_reg_m[[2]]), ]
                # test the CNVs one against all the other
                for (ii in 1:nrow(cnvs_tmp)) {
                  overlaps <- pmin(cnvs_tmp$end[ii], cnvs_tmp$end) -
                              pmax(cnvs_tmp$start[ii], cnvs_tmp$start) + 1
                  # reciprocal overlaps
                  if (!(all(overlaps >= (cnvs_tmp$end[ii]-cnvs_tmp$start[ii]+1) * prop) &
                      all(sapply(overlaps, function(x) x >= reg_list[[4]])))) {
                    pass_check <- FALSE
                    break
                  }
                }
                if (pass_check) {
                  # if all cnvs are compatible, update CNVR, remove the old ones
                  # and update DT
                  cnvrs_tmp <-
                    rbind(cnvrs_tmp,
                          data.table("r_ID" = paste0(arm, "-", n),
                          "chr" = reg_arm[1],
                          "start" = min(my_reg[[1]][2], my_reg_m[[1]][2]),
                          "end" = max(my_reg[[1]][3], my_reg_m[[1]][3])))
                  cnvrs_tmp <- cnvrs_tmp[!r_ID %in% c(my_reg[[2]], my_reg_m[[2]])]
                  DT[cnvr %in% c(my_reg[[2]], my_reg_m[[2]]), cnvr := paste0(arm, "-", n)]
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
        my_reg <- get_region(DT[ix == i, ], prop)
        # all cnvs in the cnvr
        reg_list <- get_regions_list(DT[cnvr == DT$cnvr[i]], prop)
        # this is copy-pasted from above, check twice
        overlaps <- pmin(my_reg[3], reg_list[[3]]) - pmax(my_reg[2], reg_list[[2]]) + 1
        # note that here there is a "!"
        if (!(all(overlaps >= my_reg[4]) &
              all(sapply(overlaps, function(x) x >= reg_list[[4]]))))
          DT$cnvr[i] <- NA_character_
      }
      cat(length(is.na(DT$cnvr)[is.na(DT$cnvr) == T]),
          "CNVs removed from the assigned CNVR\n")
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

get_region_with_rID <- function(my_line, prop = 1) {
  # same as get_region but returns also r_ID if present, in that case reg is a
  # list of vectors, in this way reg[[1]] remain a numeric vector
  chr <- as.integer(my_line$chr)
  st <- as.integer(my_line$start)
  en <- as.integer(my_line$end)
  len <- (en - st +1) * prop
  reg <- c(chr, st, en, len)

  # attach also r_ID if present
  if ("r_ID" %in% colnames(my_line))
    reg <- list(reg, my_line$r_ID)

  return(reg)
}

get_regions_list <- function(my_lines, prop = 1) {
  # same as get_region but with multiple lines, returns a list of vectors
  chr <- as.integer(my_lines$chr)
  st <- as.integer(my_lines$start)
  en <- as.integer(my_lines$end)
  len <- (en - st +1) * prop

  if ("r_ID" %in% colnames(my_lines))
    reg <- list(chr, st, en, len, my_lines$r_ID)
  else if ("cnvr" %in% colnames(my_lines))
      reg <- list(chr, st, en, len, my_lines$cnvr)
  else
    reg <- list(chr, st, en, len)

  return(reg)
}

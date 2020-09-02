
#' Merge adjacent CNV with equal Copy Number
#'
#' \code{merge_calls} screen the segments of a sample and merge adjacent calls
#' with equal GT if close enough
#'
#' This function takes a \code{CVNResults} object and try to merge adjacent calls
#' (having equal GT and being close
#' enough). It is designed to process one sample at the time and is integrated
#' in the function \code{\link{read_results}}. It is not suggested to use it
#' stand-alone. The merging is done starting from the larger calls, on adjacent
#' call with equal GT in a range of +/- \code{length*prop}, trough several
#' iterations until no more calls can be merged.
#'
#' @param DT_in a \code{data.table} in the format of
#'   \code{\link{read_results}} output.
#' @param prop numeric value that multiplied by the length of a call gives
#'   the range where adjacent calls are searched.
#'
#' @export
#'
#' @import data.table

# try to change the functions in that it does require SOME columns,
# not ONLY SOME columns, maybe using some vectors instead of rbind.

merge_calls <- function(DT_in, prop = 0.3) {
  # check if multiple samples are present
  if (length(unique(DT_in$sample_ID)) != 1)
    stop(paste0("Multiple samples detected! Intra results merging require ",
                "an object containig a single sample and it is recommended to perform it",
                "whithin the 'read_results' function.\n"))
  # check colnames and order
  if (!all(colnames(DT_in) == c("chr", "start", "end", "CN", "sample_ID", "GT",
                                "first_P", "last_P", "NP", "len")))
    stop("Wrong columns names or order!\n")

  # similarly to the old version, an infinite while run until necessary
  k <- 1
  while(1 == 1) {
    # cat("\nRound:", k, "\n")
    setorder(DT_in, chr, start)
    DT_out <- data.table()

    # preceed per chr
    for (cc in unique(DT_in$chr)) {
      DT_tmp <- DT_in[chr == cc, ]
      # cat("\nchr:", cc, "\n")
      skippi <- "heck"
      if (nrow(DT_tmp) == 0) next

      for (i in 1:nrow(DT_tmp)) {
        # cat("\n", i, skippi)
        if (i == nrow(DT_tmp)) {
          DT_out <- rbind(DT_out, DT_tmp[i])
          break
        }
        if (as.character(i) == as.character(skippi)) {
          # cat("already merged!\n")
          next
        }
        if (DT_tmp$GT[i] != DT_tmp$GT[i+1]) {
          DT_out <- rbind(DT_out, DT_tmp[i])
          next
        }

        # cat("\nMergable adjacent calls hit...")
        # A-----B....C-------D where A-B and C-D are two segments with equal GT
        A <- DT_tmp$start[i]
        B <- DT_tmp$end[i]
        C <- DT_tmp$start[i+1]
        D <- DT_tmp$end[i+1]
        if ((C-B+1) <= (D-A+1) * prop) {
          # do merge
          # cat("MERGE!")
          start_m <- A
          end_m <- D
          samp_m <- DT_tmp$sample_ID[i]
          # weigthed mean on NP
          cn_m <- round(mean(DT_tmp$CN[i:i+1]*
                               DT_tmp$NP[i:i+1])/sum(DT_tmp$NP[i:i+1]))
          samp_m <- DT_tmp$sample_ID[i]
          gt_m <- DT_tmp$GT[i]
          firstP_m <- min(DT_tmp$first_P[i:i+1])
          lastP_m <- max(DT_tmp$last_P[i:i+1])
          # here NP can be different from sum(tmp$NP)
          np_m <- lastP_m - firstP_m + 1
          # in the same way len can be different from sum(tmp$len)
          len_m <- end_m - start_m + 1
          if(length(cc) != 1 | length(start_m) != 1 | length(end_m) != 1 |
             length(samp_m) != 1 | length(cn_m) != 1 | length(gt_m) != 1 |
             length(firstP_m) != 1 | length(lastP_m) != 1 | length(np_m) != 1 |
             length(len_m) != 1) stop("eu eu")
          DT_out <- rbind(DT_out, data.table("chr" = cc, "start" = start_m,
                                             "end" = end_m, "sample_ID" = samp_m,
                                             "CN" = cn_m, "GT" = gt_m,
                                             "first_P" = firstP_m,
                                             "last_P" = lastP_m,
                                             "NP" = np_m, "len" = len_m))
          skippi <- i + 1
        }
        else {
          DT_out <- rbind(DT_out, DT_tmp[i])
          next
        }
      }
    }
    k <- k+1
    # cat("\n", nrow(DT_in), nrow(DT_out), "\n")
    if (nrow(DT_out) > nrow(DT_in)) stop("OUCH, somethig's really wrong here!")
    if (nrow(DT_out) == nrow(DT_in)) break
    else DT_in <- DT_out
  }
  return(DT_out)
}


## DELETE THIS

merge_calls_old <- function(DT_in, prop = 0.5) {
  # check if multiple samples are present
  if (length(unique(DT_in$sample_ID)) != 1)
    stop(paste0("Multiple samples detected! Intra results merging require ",
                "an object containig a single sample and it is recommended to perform it",
                "whithin the 'read_results' function.\n"))
  # check colnames and order
  if (!all(colnames(DT_in) == c("chr", "start", "end", "CN", "sample_ID", "GT",
                                "first_P", "last_P", "NP", "len")))
    stop("Wrong columns names or order!\n")

  # similarly to the old version, an infinite while run until necessary
  k <- 1
  while(1 == 1) {
    # re-arrange, (re)initialize indexes and stuff
    setorder(DT_in, chr, start)
    DT_in$ix <- 1:nrow(DT_in)
    unmergeble_ix <- vector()
    mergeble_segs <- list()
    n <- 1
    # re-arrange so that the longer ones come first
    setorder(DT_in, -len)
    indexes <- DT_in[, ix]
    for (i in indexes) {
      # indexes in these objects have already been used in this cycle of the 'while'
      used_ix <- c(unlist(mergeble_segs), unmergeble_ix)
      if (i %in% used_ix) next()
      # for each cycle only the adjacent ones are confronted, check if already used
      adjacent_ix <- c(i+1, i-1)[!c(i+1, i-1) %in% used_ix]
      # serch possible oversegmented cnv, note that I use GT instead of CN.
      # For merging a callB with a larger callA, callB$end must be between
      # callA$start - callA$len*prop and callA$start OR callB$start between
      # callA$end and callA$end + callA$len*prop
      tmp <- DT_in[ix %in% adjacent_ix & chr == DT_in$chr[i] & GT == DT_in$GT[i] &
                     ((start > DT_in$end[i] & start < DT_in$end[i] + DT_in$len[i]*prop) |
                        (end > DT_in$start[i] - DT_in$len[i]*prop & end < DT_in$start[i]))]
      # if no "mergeble" segments are present for DT_in[i] keep track of it
      if (nrow(tmp) == 0) unmergeble_ix <- c(unmergeble_ix, i)
      else {
        # if there segments to be merged, add a vector with the sorted ixs to the list
        mergeble_segs[[n]] <- sort(c(tmp$ix, i))
        n <- n + 1
      }
    }
    # apply changes
    if (length(mergeble_segs) == 0) {
      cat("INFO: No more segments to merge after", k, "rounds!",
          "\n#final segments: ", nrow(DT_in), "\n\n")
      break()
    }
    else{
      cat("INFO: intra results merging round", k, "\n#initial segments: ", nrow(DT_in),
          "\n#unmodified segments:", length(unique(unmergeble_ix)), "; #'mergeble' segments:",
          length(unique(unlist(mergeble_segs))), "\n#segments after merging:",
          length(unique(unmergeble_ix)) + length(mergeble_segs), "\n")
      merged_segs <- data.table()
      # for each group of mergeble segments
      for (i in 1:length(mergeble_segs)) {
        # combine the segments values
        tmp <- DT_in[ix %in% mergeble_segs[[i]], ]
        chr_m <- tmp$chr[1]
        start_m <- min(tmp$start)
        end_m <- max(tmp$end)
        samp_m <- tmp$sample_ID[1]
        # weigthed mean on NP
        cn_m <- round(mean(tmp$CN*tmp$NP)/sum(tmp$NP))
        samp_m <- tmp$sample_ID[1]
        gt_m <- tmp$GT[1]
        firstP_m <- min(tmp$first_P)
        lastP_m <- max(tmp$last_P)
        # here NP can be different from sum(tmp$NP)
        np_m <- lastP_m - firstP_m + 1
        # in the same way len can be different from sum(tmp$len)
        len_m <- end_m - start_m + 1
        tmp <- data.table("chr" = chr_m, "start" = start_m, "end" = end_m, "sample_ID" = samp_m,
                          "CN" = cn_m, "GT" = gt_m, "first_P" = firstP_m, "last_P" = lastP_m,
                          "NP" = np_m, "len" = len_m, "ix" = NA)
        merged_segs <- rbind(merged_segs, tmp)
      }
      # substitute merged segments
      DT_in <- rbind(DT_in[ix %in% unmergeble_ix, ], merged_segs)
    }
    k <- k + 1
  }
  DT_in$ix <- NULL
  return(DT_in)
}

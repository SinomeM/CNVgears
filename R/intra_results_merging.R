
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
#' @return a \code{CNVresults}, \code{DT_in} processing result.
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
                "an object containig a single sample and it is recommended to ",
                "perform it whithin the 'read_results' function.\n"))
  # check colnames and order
  if (!all(colnames(DT_in) == c("chr", "start", "end", "CN", "sample_ID", "GT",
                                "first_P", "last_P", "NP", "len")))
    stop("Wrong columns names or order!\n")

  k <- 1
  while(1 == 1) {

    setorder(DT_in, chr, start)
    DT_out <- data.table()

    # proceed per chr
    for (cc in unique(DT_in$chr)) {

      DT_tmp <- DT_in[chr == cc, ]
      merged <- 0

      if (nrow(DT_tmp) == 0) next

      for (i in 1:nrow(DT_tmp)) {

        # last call does not have any consecutive one
        if (i == nrow(DT_tmp)) {
          DT_out <- rbind(DT_out, DT_tmp[i])
          break
        }
        # if merged with the previous skip this round
        if (i == merged) next
        # check if GT is compatible
        if (DT_tmp$GT[i] != DT_tmp$GT[i+1]) {
          DT_out <- rbind(DT_out, DT_tmp[i])
          next
        }

        # A-----B....C-------D where A-B and C-D are two segments with equal GT
        A <- DT_tmp$start[i]
        B <- DT_tmp$end[i]
        C <- DT_tmp$start[i+1]
        D <- DT_tmp$end[i+1]

        if ((C-B+1) <= (D-A+1) * prop) {
          # do merge
          start_m <- A
          end_m <- D
          samp_m <- DT_tmp$sample_ID[i]
          # weighed mean on NP
          cn_m <- round(mean(DT_tmp$CN[i:i+1]*
                               DT_tmp$NP[i:i+1])/sum(DT_tmp$NP[i:i+1]))
          samp_m <- DT_tmp$sample_ID[i]
          gt_m <- DT_tmp$GT[i]
          firstP_m <- min(DT_tmp$first_P[i:i+1])
          lastP_m <- max(DT_tmp$last_P[i:i+1])
          np_m <- lastP_m - firstP_m + 1
          len_m <- end_m - start_m + 1
          # if(length(cc) != 1 | length(start_m) != 1 | length(end_m) != 1 |
          #    length(samp_m) != 1 | length(cn_m) != 1 | length(gt_m) != 1 |
          #    length(firstP_m) != 1 | length(lastP_m) != 1 | length(np_m) != 1 |
          #    length(len_m) != 1) stop("eu eu")
          DT_out <- rbind(DT_out, data.table("chr" = cc, "start" = start_m,
                                             "end" = end_m, "sample_ID" = samp_m,
                                             "CN" = cn_m, "GT" = gt_m,
                                             "first_P" = firstP_m,
                                             "last_P" = lastP_m,
                                             "NP" = np_m, "len" = len_m))
          merged <- i + 1
        }

        else {
          DT_out <- rbind(DT_out, DT_tmp[i])
          next
        }
      }
    }
    k <- k+1
    # if (nrow(DT_out) > nrow(DT_in)) stop("OUCH, somethig's really wrong here!")
    # if something has changed, do another round
    if (nrow(DT_out) == nrow(DT_in)) break
    else DT_in <- DT_out
  }
  return(DT_out)
}

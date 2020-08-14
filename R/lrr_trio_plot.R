
#' Plot markers raw data in a CNV region for a trio
#'
#' @param cnv, a \code{data.table} containig one single line of CNVs calling
#'   results.
#' @param raw_path, \code{character}, the path where raw data processed by
#'   \code{\link{read_NGS_raw}} or \code{\link{read_finalreport_raw}}.
#' @param sample_list, minimal cohort metadata, a \code{data.table} produced by the
#'   function \code{\link{read_metadt}}.
#' @param results, a \code{data.table}, the output of
#'   \code{\link{read_results}}.
#'
#' This function plot the raw data in the regions of interest in order to visually
#' confirm the presence of a good de novo CNV candidate.
#'
#' @return a \code{list} containing the plot.
#' @export
#'
#' @import data.table
#' @import ggplot2

lrr_trio_plot <- function(cnv, raw_path, sample_list, results) {
  if (!is.data.table(cnv) | !is.data.table(sample_list) |
      !is.data.table(results))
    stop("Inputs must be a 'data.table'!\n")
  if (!dir.exists(raw_path))
    stop("Directory with processed raw data don't exist, typo(s)?\n")
  if (nrow(cnv) != 1)
    stop("Please provide a single CNV (one line data.table)\n")

  # CNV infos
  off <- cnv$sample_ID
  st <- cnv$start
  en <- cnv$end
  cc <- cnv$chr
  fam <- sample_list[sample_ID == off, fam_ID]
  moth <- sample_list[fam_ID == fam & role == "mother", sample_ID]
  fath <- sample_list[fam_ID == fam & role == "father", sample_ID]

  # initial larger region to compute the rolling window mean
  len <- en - st + 1
  r_st <- st - len
  r_en <- en + len

  # load data (points and segments)
  p_off <-  load_RDS(raw_path, off, cc, r_st, r_en)
  p_moth <- load_RDS(raw_path, moth, cc, r_st, r_en)
  p_fath <- load_RDS(raw_path, fath, cc, r_st, r_en)

  s_off <- trim_res(results, off, cc, r_st, r_en)
  s_moth <- trim_res(results, moth, cc, r_st, r_en)
  s_fath <- trim_res(results, fath, cc, r_st, r_en)

  setorder(p_off, start)
  setorder(p_moth, start)
  setorder(p_fath, start)
  setorder(s_off, start)
  setorder(s_moth, start)
  setorder(s_fath, start)

  # smaller region to plot
  r_st2 <- st - round(len * 0.2)
  r_en2 <- en + round(len * 0.2)

  # lrr plots
  off_pl <- pl(p_off, cnv, "#FF9999", r_st2, r_en2) # here I pass the CNV
  moth_pl <- pl(p_moth, s_moth, "#66CC33", r_st2, r_en2) # here all segments
  fath_pl <- pl(p_fath, s_fath, "#663333", r_st2, r_en2)

  # CP values only vectors, only whithin the CNV borders
  otmp <- load_RDS(raw_path, off, cc, st, en)$copyratio
  mtmp <- load_RDS(raw_path, moth, cc, st, en)$copyratio
  ftmp <- load_RDS(raw_path, fath, cc, st, en)$copyratio

  # CP ditribution plot
  distr <- pl_distr(otmp, mtmp, ftmp)

  trio_pl <- cowplot::plot_grid(off_pl, moth_pl, fath_pl, distr, nrow = 2,
                                labels = c("A", "B", "C", "D"))

  return(trio_pl)
}

load_RDS <- function(path, samp, cc, st, en) {
  res <- readRDS(file.path(path, paste0(samp, "_chr", cc, ".rds")))[
                    start >= st & end <= en, ]
  return(res)
}

trim_res <- function(DT, samp, cc, st, en) {
  res <- DT[sample_ID == samp & chr == cc, ][
              between(start, st, en, incbounds = T) |
              between(end, st, en, incbounds = T), ]
  return(res)
}

pl <- function(pp, segs, ccol, st, en) {

    # sliding window mean, at the moment these values are not changable
    sstep <- 5
    wind <- 10
    mmean <- evobiR::SlidingWindow("mean", pp$copyratio, wind, sstep)
    # recreate coordinates for the new points
    pp[, center := (end + start)/2 ]
    pp[, cr := c(rep(mmean, each = sstep), rep(NA, nrow(pp)-(length(mmean)*sstep)) )]

    # start of first seg and end of last must be in range
    if (nrow(segs) != 0) {
      if (segs[nrow(segs), end] > en) segs[nrow(segs), end := en]
      if (segs[1, start] < st) segs[1, start := st]
    }

    # reduce the region
    pp <- pp[start >= st & end <= en, ]

    pl_out <- ggplot() +
      geom_point(data = pp, aes(x = center/1000000, y = copyratio),
                 size = 1, colour = ccol) +
      geom_segment(data = segs,
                   aes(x = start/1000000, y = CN/2, xend = end/1000000, yend = CN/2),
                   colour = "red", size = 0.8) +
      geom_line(data = pp, aes(x = center/1000000, y = cr)) +
      xlab("Position (Mbp)") +
      ylab("CopyRatio") +
      theme_bw() +
      scale_x_continuous(labels = scales::comma) +
      ylim(-0.1, 2.1)

    return(pl_out)
}

pl_distr <- function(otmp, mtmp, ftmp) {
  res <- ggplot() +
           geom_density(aes(otmp), colour = "#FF9999") +
           geom_vline(aes(xintercept = mean(otmp)), colour = "#FF9999") +
           geom_density(aes(ftmp), colour = "#663333") +
           geom_vline(aes(xintercept = mean(ftmp)), colour = "#663333") +
           geom_density(aes(mtmp), colour = "#66CC33") +
           geom_vline(aes(xintercept = mean(mtmp)), colour = "#66CC33") +
           geom_vline(aes(xintercept = 1), colour = "orange", linetype="dotted") +
           geom_vline(aes(xintercept = 1.25), colour = "orange", linetype="dotted") +
           geom_vline(aes(xintercept = 0.75), colour = "orange", linetype="dotted") +
           xlab("CopyRatio") +
           ylab("Density") +
           theme_bw() +
           xlim(-0.1, 2.1)
  return(res)
}

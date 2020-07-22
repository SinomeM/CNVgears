

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

  off <- cnv$sample_ID
  st <- cnv$start
  en <- cnv$end
  cc <- cnv$chr
  fam <- sample_list[sample_ID == off, fam_ID]
  moth <- sample_list[fam_ID == fam & role == "mother", sample_ID]
  fath <- sample_list[fam_ID == fam & role == "father", sample_ID]
  len <- en - st + 1
  bb <- round(len * 0.2)
  r_st <- st - bb
  r_en <- en + bb

  p_off <- readRDS(file.path(raw_path, paste0(off, "_chr", cc, ".rds")))[
                                          start >= r_st & end <= r_en, ]
  p_moth <- readRDS(file.path(raw_path, paste0(moth, "_chr", cc, ".rds")))[
                                            start >= r_st & end <= r_en, ]
  p_fath <- readRDS(file.path(raw_path, paste0(fath, "_chr", cc, ".rds")))[
                                            start >= r_st & end <= r_en, ]

  s_off <- results[sample_ID == off & chr == cc, ][between(start, r_st, r_en, incbounds = T) |
                                                     between(end, r_st, r_en, incbounds = T), ]
  s_moth <- results[sample_ID == moth & chr == cc, ][between(start, r_st, r_en, incbounds = T) |
                                                       between(end, r_st, r_en, incbounds = T), ]
  s_fath <- results[sample_ID == fath & chr == cc, ][between(start, r_st, r_en, incbounds = T) |
                                                       between(end, r_st, r_en, incbounds = T), ]

  setorder(p_off, start)
  setorder(p_moth, start)
  setorder(p_fath, start)
  setorder(s_off, start)
  setorder(s_moth, start)
  setorder(s_fath, start)

  pl <- function(pp, segs, ccol) {
    # sliding window mean,
    sstep <- 5
    wind <- 10
    mmean <- evobiR::SlidingWindow("mean", pp$copyratio, wind, sstep)
    # recreate coordinates for the new points
    pp[, center := (end + start)/2 ]
    # sst <- pp[, start][seq(1, by = sstep, to = nrow(pp))]
    # een <- pp[, start][seq(wind, by = sstep, to = nrow(pp))]
    #
    # coord <- pp[, center][seq(1, by = sstep, to = nrow(pp))]
    # DT <- data.table("st" = sst[1:(length(sst)-2)],
    #                  "en" = een[1:(length(een))],
    #                  "cr" = mmean,
    #                  "pos" = coord[2:(length(coord)-1)])
    pp[, cr := c(rep(mmean, each = sstep), rep(NA, nrow(pp)-(length(mmean)*sstep)) )]


    # start of first seg and end of last must be in range
    if (nrow(segs) != 0) {
      if (segs[nrow(segs), end] > r_en) segs[nrow(segs), end := r_en]
      if (segs[1, start] < r_st) segs[1, start := r_st]
    }

    pl_out <- ggplot() +
      geom_point(data = pp, aes(x = center/1000000, y = copyratio),
                 size = 1, colour = ccol) +
      geom_segment(data = segs,
                   aes(x = start/1000000, y = CN/2, xend = end/1000000, yend = CN/2),
                   colour = "red", size = 0.8) +
      geom_line(data = pp, aes(x = center/1000000, y = cr)) +
      # geom_line(data = DT, aes(x = pos/1000000, y = cr),
      #           colour = "#333333", size = 0.4) +
      # geom_segment(data = DT, aes(x = st/1000000, y = cr, xend = en/1000000, yend = cr),
      #           colour = "#333333", size = 0.4) +
      xlab("Position (Mbp)") +
      ylab("CopyRatio") +
      theme_bw() +
      # theme(panel.border=element_blank())+
      scale_x_continuous(labels = scales::comma) +
      ylim(-0.1, 2.1)

    return(pl_out)
  }

  otmp <- p_off$copyratio
  mtmp <- p_moth$copyratio
  ftmp <- p_fath$copyratio

  distr <- ggplot() +
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

  off_pl <- pl(p_off, s_off, ccol = "#FF9999")
  moth_pl <- pl(p_moth, s_moth, ccol = "#66CC33")
  fath_pl <- pl(p_fath, s_fath, ccol = "#663333")

  trio_pl <- cowplot::plot_grid(off_pl, moth_pl, fath_pl, distr, nrow = 2,
                                labels = c("A", "B", "C", "D"))

  return(trio_pl)
}


#' Discover good de novo CNVs candidates
#'
#' \code{denovo_cnvs} screen CNVs calling results to find de novo candidates,
#' given a family based study design.
#'
#' This function
#'
#' @param sample_list, a \code{data.table}, the output of
#'   \code{\link{read_metadt}}.
#' @param markers, a \code{data.table}, the output of
#'   \code{\link{read_NGS_intervals}} or \code{\link{read_finalreport_snps}},
#'   depending on the initial data type.
#' @param results, a \code{data.table}, the output of
#'   \code{\link{read_results}}.
#' @param raw_path, \code{character}, the path where raw data processed by
#'   \code{\link{read_NGS_raw}} or \code{\link{read_finalreport_raw}}.
#'
#' @export
#'
#' @import data.table


cnvs_inheritance <- function(sample_list, markers, results, raw_path,
                             mmethod = 1, alfa = 0.05, min_NP = 10,
                             pre_fine_screen = TRUE, adjust_pval = T,
                             reciprocal_overlap = 0.3) {
  # check input
  if (!is.data.table(results) | !is.data.table(sample_list) |
      !is.data.table(markers))
    stop("Inputs must be a 'data.table'!\n")
  if (!dir.exists(raw_path))
    stop("Directory with processed raw data don't exist, typo(s)?\n")
  if (!"sibling" %in% unique(sample_list[, role]) &
      !"proband" %in% unique(sample_list[, role]))
    stop("No offspring (proband | sibling) in sample table!\n")
  # check also colnames

  offs <- sample_list[role %in% c("proband", "sibling"), sample_ID]
  cat("\nINFO: computing CNVs inheritance in", length(offs),"sample from",
      length(unique(sample_list[, fam_ID])), "families.\n")

  # data.table "set" and ":=" functions act by reference, I create a copy to
  # avoid modifying the original object (perhaps there is a better option?)
  DT <- copy(results)
  rm(results)

  th <- reciprocal_overlap

  # iterate for each offspring (proband | sibling)
  for (samp in offs) {
    off_cnvs <- DT[sample_ID == samp & GT != 0, ]
    if (nrow(off_cnvs) == 0) break()
    fam <- sample_list[sample_ID == samp, fam_ID]
    moth <- sample_list[fam_ID == fam & role == "mother", sample_ID]
    fath <- sample_list[fam_ID == fam & role == "father", sample_ID]
    # complete family can be missing, if a parent is not present in the
    # sample_list issue a warning and skip the trio
    if (length(moth) != 1 | length(fath) != 1) {
      DT[sample_ID == samp & GT != 0, inheritance := "incomplete_trio"]
      warning(paste0("Complete family is missing for sample, ", samp))
      next
    }
    moth_segs <- DT[sample_ID == moth, ]
    fath_segs <- DT[sample_ID == fath, ]

    # screen per chr to reduce the impact of subsetting (?)
    for (cc in unique(off_cnvs$chr)) {

      off_tmp <- off_cnvs[chr == cc, ]

      if (nrow(off_tmp) == 0) next

      moth_tmp_c <- moth_segs[chr == cc, ]
      fath_tmp_c <- fath_segs[chr == cc, ]

      for (i in 1:nrow(off_tmp)) {
        # cnv_points <- off_tmp$first_P[i]:off_tmp$last_P[i]
        GT_tmp <- off_tmp$GT[i]
        moth_tmp <- moth_tmp_c[GT == GT_tmp, ]
        fath_tmp <- fath_tmp_c[GT == GT_tmp, ]
        m <- 0
        p <- 0
        st_tmp <- off_tmp$start[i]
        en_tmp <- off_tmp$end[i]
        len_tmp <- en_tmp - st_tmp + 1

        if (nrow(moth_tmp) != 0 ) {
          for (n in 1:nrow(moth_tmp)) {
            m_st <- moth_tmp$start[n]
            m_en <- moth_tmp$end[n]
            overl <- min(en_tmp, m_en) - max(st_tmp, m_st) + 1
            if (overl >= th*len_tmp & overl >= th*(m_en-m_st+1)) {
              m <- 1
              break
            }
          }
        }
        if (nrow(fath_tmp) != 0 ) {
          for (n in 1:nrow(fath_tmp)) {
            f_st <- fath_tmp$start[n]
            f_en <- fath_tmp$end[n]
            overl <- min(en_tmp, f_en) - max(st_tmp, f_st) + 1
            if (overl >= th*len_tmp & overl >= th*(f_en-f_st+1)) {
              p <- 1
              break
            }
          }
        }

        # return result
        if (m == 1 & p == 0) inh <- "maternal"
        if (m == 0 & p == 1) inh <- "paternal"
        if (m == 1 & p == 1) inh <- "CNP/ancestral/artifact"
        if (m == 0 & p == 0) inh <- "putative_denovo"
        DT[sample_ID == samp & seg_ID == off_tmp$seg_ID[i],
                inheritance := inh]
      }
    }
  }

  cat("\nINFO: fine-screening putative de novo CNVs using intervals/SNPs raw data\n")
  # iterate for each offspring (proband | sibling)
  for (samp in offs) {
    cat("Sample #", samp, "\n")
    if (mmethod == 0) break

    # only the putative de novo CNVs needs to be processed in this way
    off_cnvs <- DT[sample_ID == samp & inheritance == "putative_denovo", ]
    if (nrow(off_cnvs) == 0) next
    fam <- sample_list[sample_ID == samp, fam_ID]
    moth <- sample_list[fam_ID == fam & role == "mother", sample_ID]
    fath <- sample_list[fam_ID == fam & role == "father", sample_ID]

    # raw data is stored per chr thus proceed per chr
    for (cc in unique(off_cnvs$chr)) {
      off_tmp <- off_cnvs[chr == cc, ]
      # check
      if (nrow(off_tmp) == 0) next
      cat("# chr", cc, "\n")

      ## CHANGE THIS ##
      if (cc == 24 | cc == 23) {
        cat("chr Y and X ignored for now ...\n")
        next
      }

      # REMEMBER the mother do not have the RDS for chrY ("24")
      off_ints <- readRDS(file.path(raw_path, paste0(samp, "_chr", cc, ".rds")))
      fath_ints <- readRDS(file.path(raw_path, paste0(fath, "_chr", cc, ".rds")))
      moth_ints <- readRDS(file.path(raw_path, paste0(moth, "_chr", cc, ".rds")))

      # # REMOVE THIS LATER
      # if (!"copyratio" %in% colnames(off_ints)) {
      #   off_ints[, copyratio := 2^log2R]
      #   fath_ints[, copyratio := 2^log2R]
      #   moth_ints[, copyratio := 2^log2R]
      # } ## I forgot to compute this in read_finalreport_raw and I'm testing on SPARK

      for (i in 1:nrow(off_tmp)) {
        # st_tmp <- off_tmp$start[i]
        # st_tmp <- off_tmp$end[i]
        f_P <- off_tmp$first_P[i]
        l_P <- off_tmp$last_P[i]
        gt <- off_tmp$GT[i]
        sid <- off_tmp$seg_ID[i]
        # markers ID (P_ID) is unique for all samples, can be used in all the
        # trio copyratio vectors of the trio for the region of interest
        off_ints_tmp <- off_ints[P_ID >= f_P & P_ID <= l_P, copyratio]
        fath_ints_tmp <- fath_ints[P_ID >= f_P & P_ID <= l_P, copyratio]
        moth_ints_tmp <- moth_ints[P_ID >= f_P & P_ID <= l_P, copyratio]

        # Initial screening: count the number of points with a GT (from the
        # copyratios) compatible with the offspring CNV, the presence of a
        # substantial proportion of such points could indicate a missed or
        # splitted call in the parent. This is enanched by the segments-level
        # screening method, i.e. 30% RECIPROCAL overlap, quite stringent in a
        # certain way, again, in particular for events inherited but splitted in
        # the offspring or in the parent or only partially inherited
        if (pre_fine_screen == TRUE) {
          if (gt == 1) {
            t_cp <- 0.75
            mlen <- length(moth_ints_tmp[moth_ints_tmp < t_cp])
            flen <- length(fath_ints_tmp[fath_ints_tmp < t_cp])
            # more than 75% (at the moment) points of a parent are compatible with
            # the offspring call GT, it could be a missed or splitted call in the
            # parent. 75% is a pretty high threshold, but given that this is only
            # a preliminary screening (for the marker-based part) I think it's OK
            # to not drop to many putative de novo events
            if (mlen >= 0.75 * length(off_ints_tmp))
              DT[sample_ID == samp & seg_ID == sid,
                      inheritance := "p.maternal"]
            if (flen >= 0.75 * length(off_ints_tmp))
              DT[sample_ID == samp & seg_ID == sid,
                      inheritance := "p.paternal"]
            if (mlen >= 0.75 * length(off_ints_tmp) &
                flen >= 0.75 * length(off_ints_tmp))
              DT[sample_ID == samp & seg_ID == sid,
                      inheritance := "p.CNP/ancestral/artifact"]
          }
          if (gt == 2){
            t_cp <- 1.25
            mlen <- length(moth_ints_tmp[moth_ints_tmp > t_cp])
            flen <- length(fath_ints_tmp[fath_ints_tmp > t_cp])
            # see above
            if (mlen >= 0.75 * length(off_ints_tmp))
              DT[sample_ID == samp & seg_ID == sid,
                      inheritance := "p.maternal"]
            if (flen >= 0.75 * length(off_ints_tmp))
              DT[sample_ID == samp & seg_ID == sid,
                      inheritance := "p.paternal"]
            if (mlen >= 0.75 * length(off_ints_tmp) &
                flen >= 0.75 * length(off_ints_tmp))
              DT[sample_ID == samp & seg_ID == sid,
                      inheritance := "p.CNP/ancestral/artifact"]
          }
        }

        # if inheritance has been update by the first screening skip the other
        # passages
        if (DT[sample_ID == samp & seg_ID == sid,
                    inheritance] != "putative_denovo")
          next

        ## Fine screening

        # Possibilities: 1. compare the two means, 2. count the points in the
        # off that are outside the region mean+/-2*SD, 3. ...

        # Approach 1, compare the two means
        if (mmethod == 1) {
          # Here a minimum number of points is required, at the moment default
          # is 10
          if (off_tmp$NP[i] < min_NP) {
            DT[sample_ID == samp & seg_ID == sid,
                    `:=` (inheritance = "too_small")]
            next
          }

          if (gt == 1) {
            # test if the mean is lower in the offspring
            mpval <- wilcox.test(off_ints_tmp, moth_ints_tmp, exact = F)$p.value
            fpval <- wilcox.test(off_ints_tmp, fath_ints_tmp, exact = F)$p.value
          }
          if (gt == 2) {
            # test if the mean is greater in the offspring
            mpval <- wilcox.test(off_ints_tmp, moth_ints_tmp, exact = F)$p.value
            fpval <- wilcox.test(off_ints_tmp, fath_ints_tmp, exact = F)$p.value
          }

          # it also possible to obtain the confidence interval and, if desired,
          # discard results also based on that (a finer equivalent of the lower
          # limit on the number of points)

          if (mpval >= alfa)
            DT[sample_ID == samp & seg_ID == sid,
                    `:=` (inheritance = "p.maternal",
                          m_pval = mpval, p_pval = fpval)]
          if (fpval >= alfa)
            DT[sample_ID == samp & seg_ID == sid,
                    `:=` (inheritance = "p.paternal",
                          m_pval = mpval, p_pval = fpval)]
          if (mpval >= alfa & fpval >= alfa)
            DT[sample_ID == samp & seg_ID == sid,
                    `:=` (inheritance = "p.CNP/ancestral/artifact",
                          m_pval = mpval, p_pval = fpval)]
          # at the moment the p-value is for the probability of being inherited
          if (mpval < alfa & fpval < alfa)
            DT[sample_ID == samp & seg_ID == sid,
                    `:=` (inheritance = "denovo",
                          m_pval = mpval, p_pval = fpval)]
        }

        # Approach 2
        if (mmethod == 2) {

          # Here a minimum number of points is required, at the moment 10
          if (off_tmp$NP[i] < 10) {
            DT[sample_ID == samp & seg_ID == sid,
                    `:=` (inheritance = "too_small")]
            next
          }

          # It is possible that this method it too stringent, two entire
          # standard deviations are a lot (in particular for events with a low
          # number of points)

          omean <- mean(off_ints_tmp)
          osd <- sd(off_ints_tmp)
          mmean <- mean(moth_ints_tmp)
          msd <- sd(moth_ints_tmp)
          fmean <- mean(fath_ints_tmp)
          fsd <- sd(fath_ints_tmp)

          # here i recycle the name mlen and flen but they are different objects
          mlen <- length(off_ints_tmp[off_ints_tmp >= mmean - 2*msd |
                                        off_ints_tmp <= mmean + 2*msd])
          flen <- length(off_ints_tmp[off_ints_tmp >= fmean - 2*fsd |
                                        off_ints_tmp <= fmean + 2*fsd])
          # less than X% (50% at the moment) points from offspring are in the
          # tails, it could be inherited
          if (mlen >= 0.5 * length(off_ints_tmp))
            DT[sample_ID == samp & seg_ID == sid,
                    inheritance := "p.maternal"]
          if (flen >= 0.5 * length(off_ints_tmp))
            DT[sample_ID == samp & seg_ID == sid,
                    inheritance := "p.paternal"]
          if (mlen >= 0.5 * length(off_ints_tmp) &
              flen >= 0.5 * length(off_ints_tmp))
            DT[sample_ID == samp & seg_ID == sid,
                    inheritance := "p.CNP/ancestral/artifact"]
          # at the moment no p-value is given with this approach
          if (mlen < 0.5 * length(off_ints_tmp) &
              flen < 0.5 * length(off_ints_tmp))
            DT[sample_ID == samp & seg_ID == sid,
                    inheritance := "denovo"]

        }

        # Approach 2 BIS (only on SD)
        if (mmethod == 2.1) {
          if (off_tmp$NP[i] < 10) {
            DT[sample_ID == samp & seg_ID == sid,
                    `:=` (inheritance = "too_small")]
            next
          }
          omean <- mean(off_ints_tmp)
          osd <- sd(off_ints_tmp)
          mmean <- mean(moth_ints_tmp)
          msd <- sd(moth_ints_tmp)
          fmean <- mean(fath_ints_tmp)
          fsd <- sd(fath_ints_tmp)

          mlen <- length(off_ints_tmp[off_ints_tmp >= mmean - msd |
                                        off_ints_tmp <= mmean + msd])
          flen <- length(off_ints_tmp[off_ints_tmp >= fmean - fsd |
                                        off_ints_tmp <= fmean + fsd])
          # less than X% (50% at the moment) points from offspring are in the
          # tails, it could be inherited
          if (mlen >= 0.5 * length(off_ints_tmp))
            DT[sample_ID == samp & seg_ID == sid,
                    inheritance := "p.maternal"]
          if (flen >= 0.5 * length(off_ints_tmp))
            DT[sample_ID == samp & seg_ID == sid,
                    inheritance := "p.paternal"]
          if (mlen >= 0.5 * length(off_ints_tmp) &
              flen >= 0.5 * length(off_ints_tmp))
            DT[sample_ID == samp & seg_ID == sid,
                    inheritance := "p.CNP/ancestral/artifact"]
          if (mlen < 0.5 * length(off_ints_tmp) &
              flen < 0.5 * length(off_ints_tmp))
            DT[sample_ID == samp & seg_ID == sid,
                    inheritance := "denovo"]

        }

      }
    }
  }
  if (adjust_pval == TRUE)
    DT[!is.na(m_pval) & !is.na(p_pval),
              `:=` (adj_m_pval = p.adjust(m_pval[!is.na(m_pval)]),
                    adj_p_pval = p.adjust(p_pval[!is.na(p_pval)])) ]

  return(DT)
}

#' Compute CNVs inheritance
#'
#' \code{cnvs_inheritance} compute CNVs inheritance pattern and search good de novo
#' CNVs candidates
#'
#' Given a trio (mother, father, offspring) this function computes inheritance patterns
#' of the offspring's CNVs. This is done both by comparing the actual calls in the
#' trio (CNVs-level) and the raw data (markers-level). For this reason it is suggested
#' to use it separately on each method and then combine the results. As an example
#' if the focus are de novo CNVs, select only those and then combine the outputs
#' using \code{\link{inter_res_merge}}. In this way it also possible to further
#' increase the confidence in the results, e.g. by filtering out all de novo candidates
#' called by a single method.
#'
#' Internally the function is structured in several steps. [...]
#'
#' @param sample_list a \code{data.table}, the output of
#'   \code{\link{read_metadt}}.
#' @param markers a \code{data.table}, the output of
#'   \code{\link{read_NGS_intervals}} or \code{\link{read_finalreport_snps}},
#'   depending on the initial data type.
#' @param results a \code{CNVresults} object, the output of
#'   \code{\link{read_results}} or \code{\link{inter_res_merge}}.
#' @param raw_path \code{character}, the path where raw data processed by
#'   \code{\link{read_NGS_raw}} or \code{\link{read_finalreport_raw}}.
#' @param mmethod fine-screening method.
#' @param alfa p value alfa value, default in 0.05.
#' @param min_NP minimum number of points to try fine screening a de novo
#'   candidate
#' @param adjust_pval logical compute adjusted p-value?
#' @param reciprocal_overlap minimum reciprocal overlap for inherited CNVs
#'   detection
#'
#' @return a \code{CNVresults} object with some additional columns: inheritance
#'   and related p-values.
#'
#' @export
#'
#' @import data.table


cnvs_inheritance <- function(sample_list, markers, results, raw_path,
                             mmethod = 1, alfa = 0.05, min_NP = 10,
                             adjust_pval = TRUE, reciprocal_overlap = 0.3) {
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
  # avoid modifying the original object
  DT <- copy(results)
  rm(results)

  th <- reciprocal_overlap

  # iterate for each offspring (proband | sibling)
  for (samp in offs) {

    # select the trio CNVs
    trio <- select_cnvs(DT, sample_list, samp)
    if (length(trio) == 1)
      if (!trio[[1]]) next
    off_cnvs <- trio[[1]]

    # screen per chromosome
    for (cc in unique(off_cnvs[, chr])) {

      off_tmp <- off_cnvs[chr == cc, ]

      if (nrow(off_tmp) == 0) next

      moth_tmp_chr <- trio[[2]][chr == cc, ]
      fath_tmp_chr <- trio[[3]][chr == cc, ]

      for (i in 1:nrow(off_tmp)) {

        GT_tmp <- off_tmp$GT[i]
        moth_tmp <- moth_tmp_chr[GT == GT_tmp, ]
        fath_tmp <- fath_tmp_chr[GT == GT_tmp, ]

        reg_tmp <- get_region(off_tmp[i], th)

        if (nrow(moth_tmp) != 0 ) m <- check_overlap(moth_tmp, reg_tmp, th)
        else m <- 0
        if (nrow(fath_tmp) != 0 ) p <- check_overlap(fath_tmp, reg_tmp, th)
        else p <- 0

        ### QUESTE IF DIVENTANO UNA FUNZIONE
        # interpret results
        if (m == 1 & p == 0) inh <- "maternal"
        if (m == 0 & p == 1) inh <- "paternal"
        if (m == 1 & p == 1) inh <- "CNP/ancestral/artifact"
        if (m == 0 & p == 0) inh <- "putative_denovo"
        ###
        ### seg_ID forse e' meglio che diventi un "ix" creato internamente alla funzione
        DT[sample_ID == samp & seg_ID == off_tmp$seg_ID[i],
                inheritance := inh]
      }
    }
  }

  # skip markers-level comparison
  if (mmethod == 0) return(DT)

  cat("\nINFO: fine-screening putative de novo CNVs using intervals/SNPs raw data\n")
  # iterate for each offspring (proband | sibling)
  for (samp in offs) {

    cat("Sample #", samp, "\n")

    # only the putative de novo CNVs needs to be processed in this way
    off_cnvs <- DT[sample_ID == samp & inheritance == "putative_denovo", ]
    if (nrow(off_cnvs) == 0) next
    fam <- sample_list[sample_ID == samp, fam_ID]
    moth <- sample_list[fam_ID == fam & role == "mother", sample_ID]
    fath <- sample_list[fam_ID == fam & role == "father", sample_ID]

    # raw data is stored per chr thus proceed per chr
    for (cc in unique(off_cnvs$chr)) {

      off_tmp <- off_cnvs[chr == cc, ]
      if (nrow(off_tmp) == 0) next
      if (cc == 24 | cc == 23) {
        cat("chr Y and X ignored for now ...\n")
        next
      }
      else cat("# chr", cc, "\n")

      # REMEMBER the mother do not have the RDS for chrY ("24")
      off_ints <- readRDS(file.path(raw_path, paste0(samp, "_chr", cc, ".rds")))
      fath_ints <- readRDS(file.path(raw_path, paste0(fath, "_chr", cc, ".rds")))
      moth_ints <- readRDS(file.path(raw_path, paste0(moth, "_chr", cc, ".rds")))

      for (i in 1:nrow(off_tmp)) {

        my_reg <- get_region(off_tmp[i])
        gt <- off_tmp$GT[i]
        sid <- off_tmp$seg_ID[i]

        off_ints_tmp <- off_ints[start >= my_reg[2] & end <= my_reg[3], copyratio]
        fath_ints_tmp <- fath_ints[start >= my_reg[2] & end <= my_reg[3], copyratio]
        moth_ints_tmp <- moth_ints[start >= my_reg[2] & end <= my_reg[3], copyratio]

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


          ### { MAKE THIS A FUNCTION
          # check fo missing values, otherwise the test returns ERROR if there
          # are to many of them
          if (length(moth_ints_tmp[is.na(moth_ints_tmp)]) < 0.5*length(moth_ints_tmp) &
              length(fath_ints_tmp[is.na(fath_ints_tmp)]) < 0.5*length(fath_ints_tmp) &
              length(off_ints_tmp[is.na(off_ints_tmp)]) < 0.5*length(off_ints_tmp)) {
            alt <- ifelse(gt==1, "less", "greater")
            # test if the mean is less/greater in the offspring, depending on GT
            mpval <- wilcox.test(off_ints_tmp, moth_ints_tmp, exact = FALSE,
                                 alternative = alt)$p.value
            fpval <- wilcox.test(off_ints_tmp, fath_ints_tmp, exact = FALSE,
                                 alternative = alt)$p.value

            if (mpval >= alfa) inh <- "p.maternal"
            if (fpval >= alfa) inh <- "p.paternal"
            if (mpval >= alfa & fpval >= alfa) inh <- "p.CNP/ancestral/artifact"
            if (mpval < alfa & fpval < alfa) inh <- "denovo"
          }
          ### }
        }

        # Approach 2, count the points
        if (mmethod == 2) {

          # Here a minimum number of points is required, at the moment 10
          if (off_tmp$NP[i] < 10) {
            DT[sample_ID == samp & seg_ID == sid,
                    `:=` (inheritance = "too_small")]
            next
          }

          ### {
          mmean <- mean(moth_ints_tmp)
          msd <- sd(moth_ints_tmp)
          fmean <- mean(fath_ints_tmp)
          fsd <- sd(fath_ints_tmp)

          # count the points of the offspring that are within the region
          # mean+/-x*SD of both parents, where x is either 1 or 2 ATM
          if (mmmethod == 2) {
            mlen <- length(off_ints_tmp[off_ints_tmp >= mmean - 2*msd |
                                          off_ints_tmp <= mmean + 2*msd])
            flen <- length(off_ints_tmp[off_ints_tmp >= fmean - 2*fsd |
                                          off_ints_tmp <= fmean + 2*fsd])
          }
          else {
            mlen <- length(off_ints_tmp[off_ints_tmp >= mmean - msd |
                                          off_ints_tmp <= mmean + msd])
            flen <- length(off_ints_tmp[off_ints_tmp >= fmean - fsd |
                                          off_ints_tmp <= fmean + fsd])
          }
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
          ### }
        }

        ### LA data.table FINALE SI PUO' AGGIORNARE QUI!!!!
        DT[sample_ID == samp & seg_ID == sid, `:=`
           (inheritance = inh, m_pval = mpval, p_pval = fpval)]
      }
    }
  }

  if (adjust_pval == TRUE)
    DT[!is.na(m_pval) & !is.na(p_pval),
       `:=` (adj_m_pval = p.adjust(m_pval[!is.na(m_pval)]),
             adj_p_pval = p.adjust(p_pval[!is.na(p_pval)])) ]

  return(DT)
}


## SUBFUNCTIONS

select_cnvs <- function(DT, cohort, samp) {

  off_cnvs <- DT[sample_ID == samp & GT != 0, ]

  if (nrow(off_cnvs) == 0) return(list(FALSE))
  else {
    fam <- cohort[sample_ID == samp, fam_ID]
    moth <- cohort[fam_ID == fam & role == "mother", sample_ID]
    fath <- cohort[fam_ID == fam & role == "father", sample_ID]
  }
    # complete family can be missing
    if (length(moth) != 1 | length(fath) != 1) {
      DT[sample_ID == samp & GT != 0, inheritance := "incomplete_trio"]
      warning("Complete family is missing for sample, ", samp)
      return(list(FALSE))
    }
    else
      return(list(off_cnvs, DT[sample_ID == moth, ], DT[sample_ID == fath, ]))
}

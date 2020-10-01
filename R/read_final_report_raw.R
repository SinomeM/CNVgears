#' Read Illumina array raw data
#'
#' \code{read_finalreport} handles inputs of data in FinalReport-like format
#'
#' This function handles the input, pre-processing and temporary storage (as RDS
#' files) of the the markers-level raw data for each sample starting from
#' FinalReport-like files (any plain text file with columns header can be read).
#' Input must be one file per sample.
#'
#' @param DT_path, character, the path to the directory with the raw data (i.e.
#'   the splitted final reports)
#' @param rds_path, path to the directory where RDS should be stored.
#' @param pref,suff, eventual prefix an suffix (e.g. ".txt"). If not necessary
#'   must be set to \code{NA}.
#' @param sample_list, a \code{data.table}, the output of
#'   \code{\link{read_metadt}}.
#' @param markers, a \code{data.table}, the output of
#'   \code{\link{read_NGS_intervals}} or \code{\link{read_finalreport_snps}},
#'   depending on the initial data type.
#' @param chr_col, name of the column containing the chromosome information in
#'   the input file. Default is \code{"Chr"}.
#' @param pos_col, name of the column containing the SNPs position information
#'   in the input file. Default is \code{"Position"}.
#' @param LRR_col, name of the column containing the LRR information in the
#'   input file. Default is \code{"Log R Ratio"}.
#' @param BAF_col, name of the column containing the BAF information in the
#'   input file. Default is \code{"B Allele Freq"}.
#'
#' @export

read_finalreport_raw <- function(DT_path, rds_path, pref, suff,
                                 sample_list, markers, chr_col = "Chr",
                                 pos_col = "Position", LRR_col = "Log R Ratio",
                                 BAF_col = "B Allele Freq") {
  # check inputs
  if (!dir.exists(DT_path)) stop("File do not exist, typo(s)?\n")

  # create the directory for the RDS if don't exist
  dir.create(rds_path, showWarnings = FALSE)

  for (i in 1:nrow(sample_list)) {
    cat("Reading data for sample #:", i ,"\n")
    DT_path_samp <- file.path(DT_path, paste0(pref[!is.na(pref)], sample_list$sample_ID[i],
                           suff[!is.na(suff)]))
    # read file, excluding eventual header using fread()
    columns <- c(chr_col, pos_col, LRR_col, BAF_col)
    DT <- fread(DT_path_samp, skip = chr_col, select = columns)
    colnames(DT) <- c("chr", "start", "log2R", "BAF")

    # this is (MOSTLY) copy-pasted from read_NGS_raw.R, could be a function (like cnr_uniform)?
    sex <- sample_list$sex[i]
    # standardize chr
    DT <- chr_uniform(DT)
    if (sex == 2) DT <- DT[chr %in% as.character(1:23)]
    # standardize
    DT[, `:=` (start = as.integer(start), end = as.integer(start), log2R = as.numeric(log2R))]
    # col order
    DT <- DT[, .(chr, start, end, log2R, BAF)]

    # sort
    setorder(DT, chr, start)
    # compute and add these columns: P_ID (match from markers), P_CN, seg_ID, copyratio  ..?
    for (cc in unique(DT$chr)) {
      DT[chr == cc, P_ID := markers[chr == cc, ][match(DT[chr == cc, start],
                                                   markers[chr == cc, start]), P_ID]]
    }
    DT[, `:=` (P_CN = round((2^log2R)*2), copyratio = 2^log2R)]

    # # add seg_ID GT etc
    # sids <- results[sample_ID == sample_list$sample_ID[i], seg_ID]
    # for (sid in sids) {
    #   res_line <- results[sample_ID == sample_list$sample_ID[i] & seg_ID == sid]
    #   last_P <- res_line[, last_P]
    #   first_P <- res_line[, first_P]
    #   # GT <- res_line[, GT]
    #   # CN <- res_line[, CN]
    #   # DT[P_ID >= first_P & P_ID <= last_P, `:=` (seg_ID = sid, seg_GT = GT, seg_CN = CN)]
    #   DT[P_ID >= first_P & P_ID <= last_P, seg_ID := sid]
    # }

    # save in RDS one chromosome at time in the specified path
    if (sex == 1) chrs <- as.character(1:24)
    else chrs <- as.character(1:23)
    cat("Writing data for sample #:", i ,"\n")
    for (cc in chrs) {
      rds_file <- file.path(rds_path, paste0(sample_list$sample_ID[i], "_chr", cc, ".rds"))
      # cat("INFO: saving", rds_file, "\n")
      saveRDS(DT[chr == cc, ], file = rds_file)
    }
  }
}

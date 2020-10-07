#' Read raw copyratio/LRR data for NGS intervals
#'
#' \code{read_NGS_raw}
#'
#' This function handles the input, pre-processing and temporary storage (as RDS
#' files) of the the markers-level raw data for each sample. Those are mostly
#' required for CNVs inheritance detection. Any marker-level log2R-like measure
#' is supported, given that, for somatic chromosomes:
#' $$LRR = log2R = log2(CopyRatio) = log2(numeric_CN/2*)$$
#' Input must be one file per sample.
#'
#' @param DT_path, path to the input file.
#' @param chr_col, name of the column containing the chromosome information in
#'   the input file.
#' @param pref,suff, eventual prefix an suffix (e.g. ".txt"). If not necessary
#'   must be set to \code{NA}.
#' @param start_col, name of the column containing the start information in the
#'   input file.
#' @param end_col, name of the column containing the end information in the
#'   input file.
#' @param raw_type, character, it describes the format in which the LRR-like
#'   data is give. Can be either "log2R", "numeric_CN" or "copyratio", see
#'   description.
#' @param sample_list, a \code{data.table}, the output of
#'   \code{\link{read_metadt}}.
#' @param markers, a \code{data.table}, the output of
#'   \code{\link{read_NGS_intervals}} or \code{\link{read_finalreport_snps}},
#'   depending on the initial data type.
#' @param rds_path, path to the directory where RDS should be stored.
#' @param raw_col, name of the column containing the marker-level raw data in
#'   the input file.
#'
#' @return nothing, this function saves the results on disk.
#'
#' @export
#'
#' @import data.table

# LRR | log2CopyRatio = log2(numeric_CN/2*)
# with Illumina array LRR is given, with NGS-based pipeline I have seen either
# LRR, numeric_CN/2 and numeric_CN.
# Honestly I think the easiest thing is to convert everything in LRR to have a
# standardized format. Another option is to use the Copy Ratio that can be visualized
# easily.
# *: a semplification, it should be "observed value/expected value"

# num_CN  num_CN/2  log2(num_CN/2)
# 0       0         -Inf
# 1       0.5       -1
# 2       1         0
# 3       1.5       ~0.585
# 4       2         1
# 5       2.5       ~1.322
# 6       3         ~1.585


# read and save RDS for all the samples in sample_list

read_NGS_raw <- function(DT_path, rds_path, pref, suff,
                         chr_col, start_col, end_col,
                         raw_col, raw_type,
                         sample_list, markers) {
  # check input
  if (missing(DT_path) | missing(chr_col) | missing(start_col) | missing(end_col) |
      missing(raw_col) | missing(raw_type) | missing(sample_list) | missing(rds_path) |
      missing(pref) | missing(suff) | missing(markers))
    stop("Missing parameters!\n")
  if (!raw_type %in% c("log2R", "copyratio", "numeric_CN"))
    stop("Wrong 'raw_type' format!\n")
  # create the directory for the RDS if don't exist
  dir.create(rds_path, showWarnings = FALSE)

  for (i in 1:nrow(sample_list)) {
    cat("Reading data for sample #:", i ,"\n")
    DT_path_samp <- file.path(DT_path, paste0(pref[!is.na(pref)], sample_list$sample_ID[i],
                           suff[!is.na(suff)]))
    # read file, excluding eventual header using grep whithin fread()
    columns <- c(chr_col, start_col, end_col, raw_col)
    DT <- fread(DT_path_samp, skip = chr_col, select = columns)
    # (CHANGE HERE IF WE WANT TO USE COPY RATIOs)
    # the "raw" column is called log2R in any case here, I correct the content dowstream
    colnames(DT) <- c("chr", "start", "end", "log2R")

    sex <- sample_list$sex[i]
    # standardize chr
    DT <- chr_uniform(DT)
    if (sex == 2)
      DT <- DT[chr %in% as.character(1:23)]

    # standardize
    DT[, `:=` (start = as.integer(start), end = as.integer(end), log2R = as.numeric(log2R))]
    # sort
    setorder(DT, chr, start)

    # standardize to log2 if necessary (CHANGE HERE IF WE WANT TO USE COPY RATIOs)
    # somatic chrs
    if (raw_type == "copyratio")
      DT[chr %in% as.character(1:24), `:=` (copyratio = log2R, log2R = log2(log2R))]
    if (raw_type == "log2R")
      DT[chr %in% as.character(1:24), copyratio := 2^log2R]
    if (raw_type == "numeric_CN") {
      # numeric CN can't be < 0 it's a normalization artifact, eventually normalize to 0
      DT[log2R < 0, log2R := 0]
      DT[chr %in% as.character(1:22), `:=` (copyratio = log2R/2, log2R = log2(log2R/2))]
      # sex chrs
      if (sex == 1)
        DT[chr %in% as.character(23:24), `:=` (copyratio = log2R/1, log2R = log2(log2R/1))]
                                            # useless to divide by 1, I keep it for clarity
      if (sex == 2)
        DT[chr == "23", `:=` (copyratio = log2R/2, log2R = log2(log2R/2))]
    }

    # compute and add these columns: P_ID (match from markers), P_CN, seg_CN, (seg_ID),  ..?
    for (k in unique(DT$chr)) {
      DT[chr == k, P_ID := markers[chr == k, ][match(DT[chr == k, start],
                                                   markers[chr == k, start]), P_ID]]
    # this step is uneccesary if we are sure that the number of markers (in "raw data") is exactly
    # the same in every sample and exactly the same of "markers" (that is used to annotate the calls)
    }
    DT[, P_CN := round((2^log2R)*2)]

    # save in RDS one chromosome at time in the specified path
    if (sex == 1) chrs <- as.character(1:24)
    else chrs <- as.character(1:23)
    cat("Writing data for sample #:", i ,"\n")
    for (n in chrs) {
      rds_file <- file.path(rds_path, paste0(sample_list$sample_ID[i], "_chr", n, ".rds"))
      # cat("INFO: saving", rds_file, "\n")
      saveRDS(DT[chr == n, ], file = rds_file)
    }
  }
}

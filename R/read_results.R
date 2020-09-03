#' Read CNVs calling or semgentation results
#'
#' \code{read_results} takes the results of a CNVs calling pipeline and return
#' them in a standardized object.
#'
#' This function aims to convert a variety of possible types of CNVs
#' calling/segmentation pipelines and/or algorithms results into a standardized
#' format in order to easily integrate with the other functions in this
#' package. Currently two main files type and two main file-organization structures
#' are considered, for a total of four generic situations:
#' \itemize{
#' \item VCF files, one per sample (e.g. the results of GATK gCNV pipeline);
#' \item VCF file, all sample of a cohort in the same file (not yet fully implemented);
#' \item TSV/CSV file, one file per sample (e.g. the results of GATK ModSeg pipeline,
#'   or the results of running "manually" PennCNV);
#' \item TSV/CSV file, all samples of a cohort in the same file (e.g. the results
#'   of EnsembleCNV).
#' }
#' If multiple files containing results for multiple samples are present (e.g. the
#' results of PennCNV joint calling on trios) at the moment it is recommended that
#' the user concatenated those individual file in a single one prior loading them
#' with \code{read_results}.
#' Note that any line occurring before the columns header are automatically skipped
#' by \code{fread}.
#'
#' @param DT_path, path to the directory containing the individual files, if
#'   \code{res_type} is set to "directory" or to the single file if \code{res_type}
#'   is set to "file".
#' @param pref,suff, eventual prefix an suffix (e.g. ".txt") to the files to be used
#'   when \code{res_type} is set to "directory". If not necessary must be set to
#'   \code{NA}.
#' @param sample_list, minimal cohort metadata, a \code{data.table} produced by the
#'   function \code{\link{read_metadt}}.
#' @param res_type, can be either "directory" or "file", indicates whether the
#'   function must expect a single file for all samples or one file per sample.
#' @param DT_type, can be either "VCF" or "TSV/CSV", indicate the file type.
#' @param chr_col, name of the column containing the chromosome information
#'   in the input data.
#' @param start_col, name of the column containing the start information in
#'   the input data.
#' @param end_col, name of the column containing the end information in the
#'   input data.
#' @param CN_col, name of the column containing the Copy Number information
#'   in the input data.
#' @param samp_ID_col, name of the column containing the sample ID information in
#'   the input file, required if \code{res_type} is set to "file".
#' @param markers, a \code{data.table} containing the marker list, the output
#'   \code{\link{read_finalreport_snps}} with \code{DT_type} set to "markers"
#'   or \code{\link{read_NGS_intervals}}.
#' @param end_vcf, name of the field containing the segment end information in the
#'   VCF file(s), passed to the function \code{\link{read_vcf}}.
#' @param CN_vcf, name of the field containing the segment copy number
#'   information in the VCF file(s), passed to the function \code{\link{read_vcf}}.
#' @param do_merge, logical, indicates whether the function \code{\link{merge_calls}}
#'   should be automatically called for each sample (strongly suggested).
#' @param method_ID, character identifying the method (algorithms/pipeline), one letter
#'    code is strongly encouraged (e.g. "P" for PennCNV and "M" for GATK ModSeg).
#'    Numeric are converted to character.
#'
#' @export
#'
#' @import data.table

# When it everything will be OK it could be nice to parallelize
# the data input

# consider describing better the path etc

read_results <- function(DT_path, res_type, DT_type, pref = NA, suff = NA,
                         sample_list, markers,
                         chr_col, start_col, end_col, CN_col, samp_ID_col,
                         end_vcf = "END", CN_vcf = "CN",
                         do_merge = TRUE, merge_prop = 0.5, method_ID) {
  # check parameters
  if (missing(DT_path) | missing(res_type) | missing(DT_type))
    stop("Missing inputs!\n")
  if (missing(markers)) stop("Missing markers object!\n")
  if (missing(sample_list)) stop("Missing sample_list object!\n")
  if (!res_type %in% c("file", "directory")) stop("Wrong res_type format!\n")
  if (!DT_type %in% c("VCF", "TSV/CSV")) stop("Wrong DT_type format!\n")
  if (res_type == "directory" & (missing(pref) | missing(suff)))
    stop("res_type set to 'directory' but no 'pref' | 'suff' specified!\n")
  if (DT_type == "TSV/CSV" & (missing(chr_col) | missing(start_col) |
      missing(end_col))) stop("Missing parameters for DT_type == 'TSV/CSV'!\n")
  if (!is.logical(do_merge))
    stop("Wrong 'do_merge' format!\n")
  if (missing(method_ID))
    stop("Please specify an Id for this method. One letter code is encouraged!\n")

  # One single file for the resutls of all the samples ------------------------
  # if res_type is "file" one single file is expected, then the additional
  # columns and the data standardization must be computed per sample
  if (res_type == "file") {
    if (DT_type == "TSV/CSV") {
      if (missing(samp_ID_col))
        stop("res_type 'directory' selected but no samp_ID_col provided!\n")

      columns <- c(chr_col, start_col, end_col, CN_col, samp_ID_col)
      DT <- fread(DT_path, skip = chr_col, select = columns)
      colnames(DT) <- c("chr", "start", "end", "CN", "sample_ID")
    }
    if (res_type == "VCF") {
      ## TO BE TESTED ##
      DT <- read_vcf(DT_path = DT_path, samples = sample_list$sample_ID)
    }

    DT_out <- data.table()
    for (i in 1:nrow(sample_list)) {
      sid <- sample_list$sample_ID[i]
      sex <- sample_list$sex[i]

      if (nrow(DT[sample_ID == sid, ]) > 0) {
        # HERE if DT_in = DT[sample_ID == sid, ] has no row there is an error
        tmp <- DT_uniform_internal(DT_in = DT[sample_ID == sid, ],
                                   markers = markers, sex = sex)
        # merge segments
        if (do_merge)
          tmp <- merge_calls(tmp, prop = merge_prop)

        # add seg_ID for this sample
        tmp[, seg_ID := 1:nrow(tmp)]
        DT_out <- rbind(DT_out, tmp)
      }
    }
  }

  # One file per sample -------------------------------------------------------
  # if res_type is "directory"  everything can be computed per file and
  # the results concatenated in a single data.table
  if (res_type == "directory") {
      DT_out <- data.table()
      for (i in 1:nrow(sample_list)) {
        path <- file.path(DT_path, paste0(pref[!is.na(pref)],
                       sample_list$sample_ID[i], suff[!is.na(suff)]))
        if (DT_type == "TSV/CSV") {
          # to be tested
          columns <- c(chr_col, start_col, end_col, CN_col)
          # if there is an header like in GATK results, skip it
          DT <- fread(path, skip = chr_col, select = columns)[
                  , sample_ID := sample_list$sample_ID[i]]
        }
        if (DT_type == "VCF") {
          DT <- read_vcf(DT_path =  path, end_vcf = end_vcf, CN_vcf = CN_vcf)[
                  , sample_ID := sample_list$sample_ID[i]]
        }
        colnames(DT) <- c("chr", "start", "end", "CN", "sample_ID")

        # uniform data
        DT <- DT_uniform_internal(DT_in = DT, markers = markers,
                                  sex = sample_list$sex[i])
        # merge segments
        if (do_merge)
          DT <- merge_calls(DT, prop = merge_prop)

        # add seg_ID for this sample
        DT[, seg_ID := 1:nrow(DT)]

        DT_out <- rbind(DT_out, DT)
      }
  }
  # finally add the method_ID column
  DT_out[, meth_ID := as.character(method_ID)]
  # return
  class(DT_out) <- c("CNVresults", class(DT_out))
  return(DT_out)
}


# DT uniform subfunction
DT_uniform_internal <- function(DT_in, markers, sex) {
  if (missing(markers) | missing(sex))
    stop("Missing parameter to DT_uniform!\n")

  # Input: "chr" "start" "end" "CN" ;
  # Compute "GT" and "NP", sort and check columns data type
  DT_in[, `:=` (chr = as.character(chr),
                start = as.integer(start),
                end = as.integer(end))]

  DT_in <- chr_uniform(DT_in)

  # Some segmentation algorithms might use "+"/"-" notation for CN calls
  # for duplications/deletions without a precise CN, at the moment these
  # are converted in 3/1 but this might change
  # Here I assume also that everything not "+" or "-" is the normal CN indicator
  DT_in[, CN := gsub(" ", "", CN)]
  if ("+" %in% unique(DT_in$CN) | "-" %in% unique(DT_in$CN)) {
    # Autosomes
    DT_in[chr %in% as.character(1:22) & CN != "+" & CN != "-", CN := "2"][
      chr %in% as.character(1:22) & CN == "+", CN := "3"][
        chr %in% as.character(1:22) & CN == "-", CN := "1"]
    if (sex == 1){
      DT_in[chr %in% as.character(23:24) & CN != "+" & CN != "-", CN := "1"][
        chr %in% as.character(23:24) & CN == "+", CN := "2"][
          chr %in% as.character(23:24) & CN == "-", CN := "0"]
    }
    if (sex == 2) {
      # also drop calls in Y
      DT_in <- DT_in[chr %in% as.character(1:23)][
        chr == "23" & CN != "+" & CN != "-", CN := "2"][
          chr == "23" & CN == "+", CN := "3"][
            chr == "23" & CN == "-", CN := "1"]
    }
  }
  # now it is possible to check data format
  DT_in[, CN := as.integer(CN)]
  # GT = 0 (normal CN), GT = 1 (deletion), GT = 2 (duplication)
  # Autosomes
  DT_in[chr %in% as.character(1:22) & CN == 2, GT := 0][
    chr %in% as.character(1:22) & CN < 2, GT := 1][
      chr %in% as.character(1:22) & CN > 2, GT := 2]
  # Sex chromosomes
  if (sex == 1) {
    # male
    DT_in[chr %in% as.character(23:24) & CN == 1, GT := 0][
      chr %in% as.character(23:24) & CN == 0, GT := 1][
        chr %in% as.character(23:24) & CN > 1, GT := 2]
  }
  if (sex == 2) {
    # female
    # drop calls in Y
    DT_in <- DT_in[chr %in% as.character(1:23)][
      chr == "23" & CN == 2, GT := 0][
        chr == "23" & CN < 2, GT := 1][
          chr == "23" & CN > 2, GT := 2]
  }
  # Find mark_ID of the fisrt and last markers and compute "NP"
  for (cc in unique(DT_in$chr)) {
    # intricate but works correctly
    DT_tmp <- DT_in[chr == cc, ]
    mm_tmp <- markers[chr == cc, ]
    DT_in[chr == cc, `:=` (first_P = mm_tmp[match(DT_tmp[, start],
                                                  mm_tmp[, start]), P_ID],
                           last_P = mm_tmp[match(DT_tmp[, end],
                                                 mm_tmp[, end]), P_ID])]
  }
  setorder(DT_in, chr, start)
  DT_in[,`:=` (NP = last_P - first_P + 1, len = end - start + 1)]

  return(DT_in)
}

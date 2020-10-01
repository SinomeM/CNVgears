

## Old function no longer needed. This file will be deleted soon.


#' Uniform data i/o
#'
#' \code{DT_uniform} returns the input (CNVs calling/segmentation results,
#' markers list or raw data) in a standardized format.
#'
#' This function handles all the initial data manipulation in order to obtain a
#' standardized format regardless of the source (array data or NGS), both in
#' terms of actual data columns presence, and in terms of columns names and
#' content format and type, it represent the first step of all the possible
#' analysis. It uses \code{data.table} to speed up the computations on large
#' dataset, such as the markers on a SNPs array (usually between 500k and
#' 1,500k) or the merged results of a CNV calling/segmentation pipeline on
#' thousands of samples. Depending on the parameter \code{DT_type} the final
#' \code{data.table} is consituted of the following columns:
#' \itemize{
#' \item \code{chr} (as character), \code{start} (as interger), \code{end} (as
#' integer), \code{GT} (genotype, as integer, 0 for normal CN, 1 for deletions,
#' and 2 for duplications), \code{CN} (Copy Number, as integer), \code{segID}
#' (an univoque integer used as ID of the segment/ CNV), \code{NP} (Number of
#' Points, SNPs or intervals, that constitute the segment/CNV),
#' \code{first_P_ID} (integer ID of the fisrt point of the segment/CNV),
#' \code{last_P_ID}, \code{len} (length of the segment/CNV in bp), \code{sample}
#' (the sample ID), when \code{DT_type} is "results_array" or "results_NGS";
#' \item \code{chr} (as character), \code{start} (as interger), \code{end} (as
#' integer), \code{ID}, \code{name}, when \code{DT_type} is "markers_array" or
#' "markers_NGS";
#' \item \code{chr} (as character), \code{start} (as interger),
#' \code{end} (as integer), \code{ID}, \code{LRR}, \code{BAF}, \code{mark_ID},
#' when \code{DT_type} is "raw_array";
#' \item \code{chr} (as character), \code{start} (as interger),
#' \code{end} (as integer), \code{ID}, raw value depending on
#' \code{raw_col_type}, when \code{DT_type} is "raw_NGS";
#' }
#'
#' This function handles the standardization of most the main data type in this
#' package: the CNV calling /segmentation results, the markers table, the
#' finalreports-like tables for array data and the raw LRR-like values for NGS
#' data. It is intended to be used whithin other functions that handles more
#' directly the data input/output. Please note that each possible \code{DT_type}
#' has its own specific requirements both in terms of actual columns in the data
#' being passed in, and in terms of other objects required in order to being
#' processed.
#'
#' @param DT_in, input data as \code{data.table}.
#' @param DT_type, type of input to process, can be either:
#'   "results", "markers_array", "markers_NGS", "raw_array", "raw_NGS". Define
#'   the function's behaviour and the subsequent required parameters.
#' @param chr_col, name of the column containig the chromosome information
#'   in the input data, required by the following \code{DT_type}:
#'   "results", "markers_array", "markers_NGS", "raw_array", "raw_NGS".
#' @param start_col, name of the column containig the start information in
#'   the input data, required by the following \code{DT_type}: "results_array",
#'   "results_NGS", "markers_array", "markers_NGS", "raw_array", "raw_NGS".
#' @param end_col, name of the column containig the end information in the
#'   input data, required by the following \code{DT_type}: "results_array",
#'   "results_NGS", "markers_NGS", "raw_NGS".
#' @param CN_col, name of the column containig the Copy Number information
#'   in the input data, required by the following \code{DT_type}:
#'   "results".
#' @param LRR_col, name of the column containig the LRR (Log R Ratio)
#'   information in the input data, required by the following \code{DT_type}:
#'   "raw_array".
#' @param BAF_col, name of the column containig the BAF (B Allele
#'   Frquency) information in the input data, required by the following
#'   \code{DT_type}: "raw_array".
#' @param raw_col, name of the column containig the Copy Ratio or LRR-like
#'   information in the input data, required by the following \code{DT_type}:
#'   "raw_NGS".
#' @param raw_col_type, specify the format of \code{raw_col}, can be
#'   either: "CR", "LRR_like" or "log2CR". Required by the following
#'   \code{DT_type}: "raw_NGS".
#' @param markers, a \code{data.table} containing the informations on the
#'   markers used in the analysis (SNPs for arrays, intervals for NGS). It must
#'   have a specific format and it is produced by this function via the
#'   \code{DT_type} "markers_array" or "markers_NGS". It is required by the
#'   following \code{DT_type}: "results_array", "results_", "raw_NGS"
#' @param sex, Sex of the sample(s) being processed. The format is 0 for
#'   males and 1 for females. Can be a simgle \code{integer} if one single
#'   sample is present or a \code{data.frame} with the columns "ID" and "sex".
#'   It is required by the following \code{DT_type}: "results".
#'
#' @return a \code{data.table}, see the description for the details.
#'
#'
#' @import data.table

# At the moment every column but chr, start, end, CN, ID (sample or marker) is discarded from
# the original file, I have planned a keep_all param but at the moment is not a priority.

# At the moment for "markers_array" the column "mark_ID_col" is always required, this could
# (would?) be changed, since it can not be always present in the data.
# For "markers_NGS" the situation is the opposite.

# the anotation of the raw data (segID segGT etc), mostly needed for inheritance patterns
# determination is moved to the function "markers_annotator"

# "chr" columns is now a character, %in% as.character(1:24), I feel like it wasn't necessary
# to convert it to integer, at least not ALWAYS.

# Results are precessed in the same way for array and NGS: "GT" coputed from "CN" and "sex",
# "NP" computed from "start", "end" and the markers object.
# The difference will be between CNVs calling and segmentation results, in other functions.

# At the moment CN +/- is converted to 3/1


# IMPORTANT !!!
# Once the major "read" functions are tested this function can be scomponed and moved to the
# relative "read" function, while the description can be moved to the main man (at least in part)

dt_uniform <- function(DT_in, DT_type, markers, chr_col, start_col, end_col, CN_col,
                       LRR_col, BAF_col, raw_col, raw_col_type, mark_ID_col, sample_ID_col,
                       sample_ID, sex, keep_all = FALSE) {
  # check argument consistency
  if (!is.data.table(DT_in)) stop("DT_in is not a data.table object!\n")
  if (!DT_type %in% c("results", "markers_array", "markers_NGS", "raw_array",
                      "raw_NGS")) stop("Invalid DT_type format!\n")

  ## SNP array markers & NGS intervals ----------------------------------------
  # output cols: ("SNP_ID"), "chr"*, "start"*, "end"*, "P_ID"*
  # *:integer
  if (DT_type == "markers_array" | DT_type == "markers_NGS") {
    # if NGS, "mark_ID_col" must be explicitly NA, if array "end_col" must be explicitly NA
    if (missing(mark_ID_col) | missing(chr_col) | missing(start_col) |
        missing(end_col)) stop("Missing parameter to DT_uniform!\n")
    # select the three necessary columns from the input, set standard columns names
    # and data format, sort and assign a numeric ID, only three step are different
    # between "markers_array" and "markers_NGS"
    if (DT_type == "markers_array") columns <- c(mark_ID_col, chr_col, start_col)
    if (DT_type == "markers_NGS") columns <- c(chr_col, start_col, end_col)
    DT_in <- DT_in[ , ..columns]
    if (DT_type == "markers_array") colnames(DT_in) <- c("SNP_ID", "chr", "start")
    if (DT_type == "markers_NGS") colnames(DT_in) <- c("chr", "start", "end")

    DT_in <- chr_uniform(DT_in)

    # for SNPs array "start" and "end" are the same number, it is reduntant but it makes
    # thinfs easier downstream
    if (DT_type == "markers_array") DT_in[, c("start", "end") := list(as.integer(start),
                                                                      as.integer(start))]
    if (DT_type == "markers_NGS") DT_in[, c("start", "end") := list(as.integer(start),
                                                                    as.integer(end))]
    setorder(DT_in, chr, start)
    DT_in$P_ID <- 1:nrow(DT_in)
  }

  ## CNVs calling/segmentation results ----------------------------------------
  if (DT_type == "results") {
    if (missing(markers) | missing(sex) | missing(sample_ID_col))
      stop("Missing parameter to DT_uniform!\n")

    # Input: "chr" "start" "end" "CN" ;
    # Compute "GT" and "NP", sort and check columns data type
    DT_in[, chr := as.character(chr)][, start := as.integer(start)][
        , end := as.integer(end)]

    DT_in <- chr_uniform(DT_in)

    # Some segmentation algorithms migth use "+"/"-" notation for CN calls
    # for duplications/deletions without a precise CN, at the moment these
    # are coverted in 3/1 but this might change
    # Here I assume also that everithing not "+" or "-" is the normal CN indicator
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
                       chr == "23" & CN == 2, GT := 0][
                         chr == "23" & CN == 2, GT := 0]
    }
    # Find mark_ID of the fisrt and last markers and compute "NP"
    for (i in as.character(1:24)) {
      # intricated but works correctly
      DT_in[chr == i, `:=` (first_P = markers[chr == i][match(DT_in[chr == i, start],
                                                              markers[chr == i, start]), P_ID],
                            last_P = markers[chr == i][match(DT_in[chr == i, end],
                                                             markers[chr == i, end]), P_ID])]
    }
    DT_in[,`:=` (NP = last_P - first_P + 1, len = end - start + 1)]
  }

  # Return results ------------------------------------------------------------
  if (DT_type == "markers_array" | DT_type == "markers_NGS" | DT_type == "results")
    return(DT_in)
}

#' Read sample file with minimal metadata
#'
#' \code{read_medaDT} handles the input of the sample table (sampleID, sex, role,
#' famID) in a standardized format.
#'
#' This function is needed in the first step of virtually every analysis. The input data
#' must have at least the following columns:
#' \itemize{
#' \item sample ID, self describing;
#' \item sex, ideally in 1/2 format, for males and females, however also "male"/"female"
#'   or "Male"/"Female" are accepted;
#' \item role, role of the sample in the family, either "father", "mother", "proband" or
#' "sibling";
#' \item family ID, self describing.
#' }
#' Actual name and order of the columns in the file is not relevant since they are
#' passed to the function via parameter.
#' Since the function in this package are optimized for family based studies, family ID and
#' role information for each sample are required, however if the user is interested only,
#' as an example, in CNVRs computation, genic content annotation or identification of calls
#' in IG regions and does not have such information "role" and famID can be "NA". Note that
#' doing so some functions won't be usable.
#'
#' @param DT_path path to the input file.
#' @param sample_ID_col name of the columns containing the sample ID in the original file;
#' @param sex_col name of the columns containing the sex information in the original file;
#' @param role_col name of the columns containing the role information ID in the original file;
#' @param fam_ID_col name of the columns containing the family ID in the original file;
#'
#' @return cohort metatadata object, a \code{data.table}. Will be of the
#'   \code{SampleList} class in future versions.
#'
#' @export

# CHANGE INPUT TO PED

read_metadt <- function(DT_path, sample_ID_col = "sample_ID", sex_col = "sex",
                        role_col = "role", fam_ID_col = "fam_ID") {
  # raed DT and drop eventual uneeded columns
  columns <- c(sample_ID_col, sex_col, role_col, fam_ID_col)
  DT <- fread(DT_path, select = columns)
  colnames(DT) <- c("sample_ID", "sex", "role", "fam_ID")

  # check sex and role
  if (!all(unique(DT$sex) %in%
           c("1", "2", "male", "female", "Male", "Female")))
    stop("Bad \"sex\" format!\n")
  if (any(unique(DT$role) %in% c("Father", "Mother", "Proband", "Sibling"))) {
    DT[role == "Father", role := "father"][role == "Mother", role := "mother"][
      role == "Proband", role := "proband"][role == "Sibling", role := "sibling"]
  }

  # if role is NA print warning, check it and skip the rest
  if (any(is.na(unique(DT$role)) | (unique(DT$role) == "NA"))) {
    cat("WARNING: \"role\" not specified!\n")
    DT[, role := NA]
  }
  else {
    if (!all(unique(DT$role) %in%
             c("father", "mother", "proband", "sibling", "Father", "Mother",
               "Proband", "Sibling"))) stop("Bad \"role\" format!\n")
    # convert "male"/"female" or "Male"/"Female" to "1"/"2" and
    # "Father", "Mother" etc to "father", "mother" if needed
    # there may be an easiest way
    if (any(unique(DT$sex) %in% c("male", "female", "Male", "Female"))) {
      DT[sex == "male" | sex == "Male", sex := 1][
        sex == "female" | sex == "Female", sex := 2]
    }
  }

  return(DT)
}

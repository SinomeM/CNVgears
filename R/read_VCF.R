#' Convert a VCF file of genomics segments into a \code{data.table}
#'
#' \code{read_vcf} read a VCF file into a \code{data.table}
#'
#' This function use \code{readVcf} from \code{VariantAnnotation} to read VCF
#' files, then it select only the necessary columns (for the purpose of CNVs
#' calling results analysis) and convert it to a \code{data.table}. Can also be
#' used to check the names of the necessary fields (end and copy number) if not
#' already known, using the parameter \code{explore}. By default it expect a
#' file containing data for a single sample (e.g. the results of gCNV from
#' GATK), but it can process files containing multiple samples if a character
#' vector containing the IDs is given to the parameter \code{samples}.
#'
#' @param DT_path path to the file.
#' @param end_vcf name of the field containing the segment end information.
#' @param CN_vcf name of the field containing the segment copy number
#'   information.
#' @param samples NA by default, if a character vector is provided is used to
#'   identify and select samples in a VCF containing multiple ones.
#' @param explore logic, \code{FALSE} by default. If \code{TRUE} the file in
#'   \code{DT_path} is not loaded, instead, several infos about the VCF fields
#'   are printed.
#'
#' @export

# Single sample,    Tested         OK!
# Multiple samples, to be finished


read_vcf <- function(DT_path, end_vcf = "END", CN_vcf = "CN", samples = NA,
                     explore = FALSE) {
  # just explore the header (e.g. if not sure on the right on col names)
  if (explore == TRUE) {
    cat("\nVCF header content\n")
    print(VariantAnnotation::scanVcfHeader(DT_path))
    cat("\nVCF 'geno' content \n")
    print(VariantAnnotation::geno(VariantAnnotation::scanVcfHeader(DT_path)))
    cat("\nVCF 'geno' content description\n")
    print(VariantAnnotation::info(VariantAnnotation::scanVcfHeader(DT_path)))
    cat("\nVCF 'fixed' content\n")
    print(VariantAnnotation::fixed(VariantAnnotation::scanVcfHeader(DT_path)))
  }
  # actually read VCF
  if (explore == FALSE) {
    # check defaults
    if (end_vcf != "END" | CN_vcf != "CN")
      cat("WARNING: default VCF field idenfier modified!\n")
    # single sample
    if (is.na(samples)) {
      cat("\nINFO: Reading VCF of a single sample\n")
      vcf <- VariantAnnotation::readVcf(DT_path,
                                        param = VariantAnnotation::ScanVcfParam(info = end_vcf,
                                                                                geno = CN_vcf))
      DT <- data.table(as.data.table(DelayedArray::rowRanges(vcf))[,1:2],
                       as.data.table(VariantAnnotation::info(vcf))[, ..end_vcf],
                       as.data.table(VariantAnnotation::geno(vcf))$value)
      colnames(DT) <- c("chr", "start", "end", "CN")
      rm(vcf)
    }
    # multiple samples in the same VCF
    if (!is.na(samples)) {
      cat("\nINFO: Reading VCF of multiple samples\n")
      # ... #
    }
    return(DT)
  }
}



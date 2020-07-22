
#' Generate blacklist for telomeric and centromeric regions
#'
#' @param DT_in, a \code{data.table} with start, and centromere position for
#'   each chromosome. For the assemblies hg18, hg19 and hg38 this is provided in
#'   the package.
#' @param telom, logical, should be telomeric blacklist be produced?
#' @param centrom, logical, should be centromeric blacklist be produced?
#' @param region, integer. How many basepairs large should the regions be? The
#'   centromeric region will be twice this large.
#'
#' @return, a \code{data.table} that can be passed to
#' \code{\link{cleaning_filter}} as blacklist.
#' @export
#'

telom_centrom <- function(DT_in, telom = TRUE, centrom = TRUE,
                          region = 50000) {
  if (!is.data.table(DT_in))
    stop("input must be a data.table")
  if (telom == TRUE) {
    if (!any(c("chr", "start", "end") %in% colnames(DT_in)))
      stop("input is not in the required format")
  }
  if (centrom == TRUE) {
    if (!any(c("chr", "centromere") %in% colnames(DT_in)))
      stop("input is not in the required format")
  }
  if (telom == FALSE & centrom == FALSE)
    warning("The function won't do anithing...")

  DT_out <- data.table()

  if (telom == TRUE) {
    for (i in 1:nrow(DT_in)) {
      tmp <- data.table("chr" = DT_in$chr[i],
                        "start" = c(DT_in$start[i], DT_in$end[i]-region),
                        "end" = c(DT_in$start[i]+region, DT_in$end[i]))
      DT_out <- rbind(DT_out, tmp)
    }
  }

  if (centrom == TRUE) {
    for (i in 1:nrow(DT_in)) {
      tmp <- data.table("chr" = DT_in$chr[i],
                        "start" = DT_in$centromere[i]-region,
                        "end" = DT_in$centromere[i]+region)
      DT_out <- rbind(DT_out, tmp)
    }
  }

  return(DT_out)
}


write_BED <- function(DT_in, name = NA, file_name) {
  # chrom chromStart chromEnd Name

  if (!is.na(name)) {
    tmp <- DT_in[, c("chr", "start", "end", name), with = FALSE][
        name := paste0()] # ... #
  }
  else {
    tmp <- DT_in[, .(chr, start, end)]

  }

}

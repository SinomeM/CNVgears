
# Collection of internal functions used by more than one functions in the package

get_region <- function(my_line, prop = 1) {
  # given a line consisting of a single CNV, returns a vector constaing chr,
  # start, end , length * prop as integers
  chr <- as.integer(my_line$chr)
  st <- as.integer(my_line$start)
  en <- as.integer(my_line$end)
  len <- (en - st +1) * prop
  reg <- c(chr, st, en, len)

  return(reg)
}

get_region_with_rID <- function(my_line, prop = 1) {
  # same as get_region but returns also r_ID if present, in that case reg is a
  # list of vectors, in this way reg[[1]] remain a numeric vector
  chr <- as.integer(my_line$chr)
  st <- as.integer(my_line$start)
  en <- as.integer(my_line$end)
  len <- (en - st +1) * prop
  reg <- c(chr, st, en, len)

  # attach also r_ID if present
  if ("r_ID" %in% colnames(my_line))
    reg <- list(reg, my_line$r_ID)

  return(reg)
}

get_regions_list <- function(my_lines, prop = 1) {
  # same as get_region but with multiple lines, returns a list of vectors
  chr <- as.integer(my_lines$chr)
  st <- as.integer(my_lines$start)
  en <- as.integer(my_lines$end)
  len <- (en - st +1) * prop

  if ("r_ID" %in% colnames(my_lines))
    reg <- list(chr, st, en, len, my_lines$r_ID)
  else if ("cnvr" %in% colnames(my_lines))
      reg <- list(chr, st, en, len, my_lines$cnvr)
  else
    reg <- list(chr, st, en, len)

  return(reg)
}

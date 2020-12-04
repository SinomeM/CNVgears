
#' Annotate genic load
#'
#' @param DT_in, a \code{data.table} consisting of CNVs calling results.
#' @param biotypes, character vector of Genecode biotypes, default value is
#'   "protein_coding".
#' @param mart, user specified \code{biomaRt::useMart()} object. Used if DT_in
#'   is assembly (e.g. hg19) require older Ensembl releases. If NA the latest
#'   release available by \code{biomaRt::useMart()} is used.
#'
#'   This function takes a \code{CNVresults} object as input and add two additional
#'   columns representing the genic content of each call, i.e. "genes" and "n_genes".
#'   The genes considered can be changed using the \code{biotypes} parameter depending
#'   on which types of genes is the user interested in. At the moment, genes a
#'   stored as a "-" separated list of Ensembl IDs.
#'
#' @return the \code{CNVresults} object \code{DT_in} with additional columns:
#'   genes and n_genes.
#'
#' @export
#'
#' @import data.table
#'
#' @examples
#' \dontrun{
#' DT <- genic_load(penn_22)
#' }

genic_load <- function(DT_in, biotypes = "protein_coding", mart = NULL) {
  if (!is.data.table(DT_in))
    stop("DT_in must be a data.table!")

  # data.table "set" and ":=" functions act by reference, I create a copy to
  # avoid modifying the original object (perhaps there is a better option?)
  DT <- copy(DT_in)
  rm(DT_in)

  if (length(mart) == 0)
  mart <- biomaRt::useMart("ensembl",
                           dataset="hsapiens_gene_ensembl")

  genes <- biomaRt::getBM(attributes=c("chromosome_name", "start_position",
                                       "end_position", "ensembl_gene_id",
                                       "gene_biotype"),
                          filters = c("biotype"), values = biotypes,
                          mart = mart)

  setDT(genes)
  genes <- genes[, gene_biotype := NULL]
  setnames(genes, c("chromosome_name", "start_position", "end_position"),
                  c("chr", "start", "end"))
  genes <- genes[chr %in% c(1:22, "X", "Y"), ]  # TEMPORARY SOLUTION !!! #
  genes <- chr_uniform(genes)
  setorder(genes, chr, start)
  setorder(DT, chr, start)

  # screen per chromosome
  DT[, ix := 1:.N]
  for (cc in unique(DT$chr)) {
    DT_tmp <- DT[chr == cc, ]
    genes_tmp <- genes[chr == cc, ]
    for (i in 1:nrow(DT_tmp)) {
      st <- DT_tmp$start[i]
      en <- DT_tmp$end[i]
      index <- DT_tmp$ix[i]
      hits <- genes_tmp[between(start, st, en) | between(end, st, en) |
                          (start < st & end > en)]
      n_hits <- nrow(hits)
      if (n_hits > 0){
        DT[ix == index, `:=` (n_genes = n_hits,
                              genes = paste0(hits$ensembl_gene_id,
                                             collapse = "-"))]
      }
    }
  }

  DT[, ix := NULL]
  return(DT)
}


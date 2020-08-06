
#' Retrieve genomic regions of consecutive immunoglobulin genes
#'
#' @param biotype, vector of character, is used to filter the genes based on
#'   \code{gene_biotype} from to \code{biomaRt::getBM}. By default all the
#'   immunoglobulin biotypes al listed in Genecode
#'   <https://www.gencodegenes.org/pages/biotypes.html>. If the user specify one
#'   or more values only those are used.
#' @param n_genes, integer, number of minim consecutive genes required to
#'   construct an "immunoglobulin region", default is 5.
#' @param mart, user provided mart (from \code{biomaRt::useMart} function). For
#'   older assemblies the user must manually retrieve the correct mart via
#'   \code{biomaRt}.
#'
#'   \code{immuno_regions} is used to construct a \code{data.table} of genomic
#'   regions where all the consecutive gene have an immunoglobulin (or user
#'   specified). The output can be used as a blacklist in
#'   \code{\link{cleaning_filter}}.
#'
#' @return a list with two element, the first is a \code{data.table} containing
#'   the actual regions (to be passed to \code{\link{cleaning_filter}} for
#'   example), the second is the full genes list (for the required biotypes),
#'   again as \code{data.table}.
#'
#' @export
#'
#' @import data.table

# returns a list, the first element can be passed to filter_calls.R and both can be
# written as BED for IGV visualization

immuno_regions <- function(biotype = NULL,
                           n_genes = 5, mart = NULL) {

  if (length(biotype) == 0)
    biotype <- c("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene",
                 "IG_V_gene", "TR_C_gene", "TR_J_gene", "TR_V_gene",
                 "TR_D_gene", "IG_pseudogene", "IG_C_pseudogene",
                 "IG_J_pseudogene", "IG_V_pseudogene",
                 "TR_V_pseudogene", "TR_J_pseudogene")

  # retrive IG genes and standardize chr
  if (length(mart) == 0) ensembl <- biomaRt::useMart("ensembl",
                                               dataset="hsapiens_gene_ensembl")
  else ensembl <- mart

  genes <- biomaRt::getBM(attributes=c("ensembl_gene_id", "chromosome_name",
							          "start_position", "end_position", "gene_biotype"),
      				     filters = "chromosome_name", values = c(1:22, "X", "Y"),
                         mart = ensembl)
  setDT(genes)
  colnames(genes) <- c("ensembl_gene_id", "chr", "start", "end", "gene_biotype")
  genes <- chr_uniform(genes)
  setorder(genes, chr, start)
  genes$id <- 1:nrow(genes)

  # compute IG region
  genes_ig <- genes[gene_biotype %in% biotype, ]
  ids <- genes_ig$id
  # search regions of consecutive IG genes
  breaks <- c(0, which(diff(ids) != 1), length(ids))
  tmp <- sapply(seq(length(breaks) - 1),
                function(i) ids[(breaks[i] + 1):breaks[i+1]])
  IG_regions <- data.table()
  for (i in 1:length(tmp)) {
    if (length(tmp[[i]]) > n_genes)
      IG_regions <- rbind(IG_regions,
                          data.table(first = tmp[[i]][1],
                                     last = tmp[[i]][length(tmp[[i]])],
                                     n_genes = length(tmp[[i]])))
  }
  IG_regions[, `:=` (chr = mapply(FUN = function(x)
                                    genes_ig$chr[genes_ig$id == x], first),
                     start = mapply(FUN = function(x)
                                      genes_ig$start[genes_ig$id == x], first),
                     end = mapply(FUN = function(x)
                                    genes_ig$end[genes_ig$id == x], last))][
                       , `:=` (first = NULL, last = NULL)]

  setorder(IG_regions, chr, start)
  setorder(genes_ig, chr, start)

  return(list(IG_regions, genes_ig))
}


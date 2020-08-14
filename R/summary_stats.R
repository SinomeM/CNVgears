#' Explore CNV calling results prior filtering
#'
#' @param object, a \code{data.table}, the output of
#'   \code{\link{read_results}}.
#' @param sample_list, a \code{data.table}, the output of
#'   \code{\link{read_metadt}}.
#' @param markers, a \code{data.table}, the output of
#'   \code{\link{read_NGS_intervals}} or \code{\link{read_finalreport_snps}},
#'   depending on the initial data type.
#' @param plots_path, path where the plots should be saved, if \code{NA} no plot is
#'   produced.
#'
#'   This function produce several summary statistics on the CNVs results in
#'   input. Ideally, it should be used interactively, together with
#'   \code{\link{cleaning_filter}}. Some information is printend on the console
#'   (mostly \code{summary()} on several characteristics of the results and the
#'   cohort), several plots can be produced and saved in the user specified
#'   location. Finally the function also return a \code{data.table} of sample-level
#'   summary.
#'
#' @return a \code{data.table} with summary statistics for the samples in the
#'   cohort.
#' @export
#'
#' @import data.table
#' @import ggplot2

# the function also returns the object sstat that contains number of CNVs etc
# per sample cab be used to decide if a sample is oversegmented and need to be
# excluded (some lines already printed int the console)

# add the option to save as RDS and PNG!

summary.CNVresults <- function(object, sample_list, markers, plots_path = NA) {
  # check input
  if (!is.data.table(object) | !is.data.table(sample_list) |
      !is.data.table(markers))
    stop("Inputs must be a 'data.table'!\n")
  if (length(unique(object$sample_ID)) != length(unique(sample_list$sample_ID)))
    cat("WARNING: number of samples in sample list and results differ!\n")
  if (!all(unique(sample_list$sample_ID) %in% unique(object$sample_ID)))
    cat("WARNING: not all the samples in sample_list are present in results.\n")

  # Very basic summary stats
  cat("\nBasic summary statistics...\n",
      "\n# samples:", nrow(sample_list),
      "\n# families:", length(unique(sample_list$fam_ID)),
      "\n# probands:", length(sample_list[role == "proband", sample_ID]),
      "\n# segments (normal genotype, deletions, duplications):", nrow(object),
      "\n# CNVs (deletions, duplications):", nrow(object[GT != 0, ]),
      "\n# deletions:", nrow(object[GT == 1, ]),
      "\n# duplications:", nrow(object[GT == 2, ]),
      "\n# CNVs per sample:",
      round(nrow(object[GT != 0, ]) / nrow(sample_list), 1), "\n")

  # Summary len & NP
  cat("\nSummary of CNVs length, len, and number of points (SNPs or intervals),",
      "NP in the results.\n", "\nCNVs length\n")
  print(summary(object[GT != 0, len]))
  cat("\nDeletions length\n")
  print(summary(object[GT == 1, len]))
  cat("\nDuplications length\n")
  print(summary(object[GT == 2, len]))
  if ("NP" %in% colnames(object)) {
    cat("\nCNVs NP\n")
    print(summary(object[GT != 0, NP]))
    cat("\nDeletions NP\n")
    print(summary(object[GT == 1, NP]))
    cat("\nDuplications NP\n")
    print(summary(object[GT == 2, NP]))
  }

  if (!is.na(plots_path)) {
    # create directory if not existent
    suppressWarnings(dir.create(plots_path))
    cat("\nSaving plots in the user specified directory...\n")
    # Basic plots #
    # CNVs per Mbp per chr (barplot)*, CNVs per NP per chr (barplot)*,
    # len distribution (kernel) of most most segmented chrs (CNVs per kbp) (?)

    # the chromosomes will not be sorted correctly FIX IT !!!!

    # Segments per chr, these two plots will be different only if the results are
    # from a segmentation algorithm
    pdf(file = file.path(plots_path, "nSeg_per_chr.pdf"), paper = "a4r")
    print(
    cowplot::plot_grid(
      ggplot( data = object, aes(x = reorder(as.character(chr), as.integer(chr)),
                                  fill = as.character(CN)) ) +
        geom_bar(colour = "black") +
        theme_bw() +
        labs(x = "chr", y = "# segments", title = "Number of segments per chromosomes") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(size=10), axis.text = element_text(size = 6),
              axis.title = element_text(size = 8), legend.text = element_text(size = 6),
              legend.title = element_text(size = 8), legend.key.size = unit(.5, "cm")),

      ggplot( data = object[GT != 0, ], aes(x = reorder(as.character(chr), as.integer(chr)),
                                             fill = as.character(CN))) +
        geom_bar(colour = "black") +
        theme_bw() +
        labs(x = "chr", y = "# segments",
             title = "Number of CNVs (segments with GT != 0) per chromosomes") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(size=10), axis.text = element_text(size = 6),
              axis.title = element_text(size = 8), legend.text = element_text(size = 6),
              legend.title = element_text(size = 8), legend.key.size = unit(.5, "cm")),
      nrow = 2)
    )
    dev.off()
    if ("NP" %in% colnames(object)) {
      # Points per chr (from intervals/markers)
      pdf(file = file.path(plots_path, "nMarkers_per_chr.pdf"), paper = "a4r")
      print(
        ggplot( data = markers, aes(x = reorder(as.character(chr),as.integer(chr)),
                                    fill = reorder(as.character(chr),as.integer(chr)))) +
          geom_bar(colour = "black") +
          theme_bw() +
          labs(x = "chr", y = "# points", title = "Number of markers points per chromosome") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                plot.title = element_text(size=10), axis.text = element_text(size = 6),
                axis.title = element_text(size = 8), legend.text = element_text(size = 6),
                legend.title = element_text(size = 8), legend.key.size = unit(.5, "cm")) +
          guides(fill=FALSE)
      )
      dev.off()
    }
    # Length and NP distributions
    pdf(file = file.path(plots_path, "lengths_and_NP_distributions.pdf"), paper = "a4r")
    print(
    cowplot::plot_grid(
      ggplot( object, aes(log10(len)) ) +
        geom_histogram(bins = 50, colour="black", aes(fill=as.character(CN))) +
        labs(x = "log10(length)", y = "# segments",
             title = "Length distribution of all segments") +
        theme_bw() +
        theme(plot.title = element_text(size=10), axis.text = element_text(size = 6),
              axis.title = element_text(size = 8), legend.text = element_text(size = 6),
              legend.title = element_text(size = 8), legend.key.size = unit(.5, "cm")),
      ggplot( object, aes(log10(NP)) ) +
        geom_histogram(bins = 50, colour="black", aes(fill=as.character(CN))) +
        labs(x = "log10(NP)", y = "# segments",
             title = "NP (number of points) distribution in all segments") +
        theme_bw() +
        theme(plot.title = element_text(size=10), axis.text = element_text(size = 6),
              axis.title = element_text(size = 8), legend.text = element_text(size = 6),
              legend.title = element_text(size = 8), legend.key.size = unit(.5, "cm")),
      nrow = 2)
    )
    dev.off()
  }

  # compute per sample:
  # CNVs, mean NP, mean len,
  cat("\nNumber of CNVs, mean length and NP computed per sample, could help",
      "identifig over-segmented samples.",
      "\nFirst and last 10 lines (most segmented samples first)...\n")

  if ("NP" %in% colnames(object))
    sstat <- object[GT != 0, .(n_cnvs = .N, mean_NP = mean(NP, na.rm = TRUE),
                                mean_len = mean(len, na.rm = TRUE)),
                     by = sample_ID]
  else
    sstat <- object[GT != 0, .(n_cnvs = .N, mean_len = mean(len, na.rm = TRUE)),
                     by = sample_ID]

  setorder(sstat, -n_cnvs)
  print(head(sstat, 10))
  print(tail(sstat, 10))
  cat("\nSummary on number of CNVs in the entire cohort..\n")
  print(summary(sstat$n_cnvs))

  if (!is.na(plots_path)) {
    cat("\nSaving plots in the user specified directory...\n")

    pdf(file = file.path(plots_path, "nCNV_per_sample_distribution.pdf"), paper = "a4r")
    print(
    ggplot(sstat, aes(x = n_cnvs)) +
      geom_histogram(bins = 50, colour = "black") +
      labs(title = "Number of CNVs per sample distibution.") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size=10), axis.text = element_text(size = 6),
            axis.title = element_text(size = 8), legend.text = element_text(size = 6),
            legend.title = element_text(size = 8), legend.key.size = unit(.5, "cm"))
    )
    dev.off()

    pdf(file = file.path(plots_path, "meanLen_on_nCNVs_per_sample.pdf"), paper = "a4r")
    print(
    ggplot(sstat, aes(x = mean_len, y = n_cnvs)) +
      geom_point() +
      theme_bw()
    )
    dev.off()
  }

  return(sstat)
}

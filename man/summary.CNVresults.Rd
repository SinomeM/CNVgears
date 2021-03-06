% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_stats.R
\name{summary.CNVresults}
\alias{summary.CNVresults}
\title{Explore CNV calling results prior filtering}
\usage{
\method{summary}{CNVresults}(object, sample_list, markers, plots_path = NA, ...)
}
\arguments{
\item{object, }{a \code{data.table}, the output of
\code{\link{read_results}}.}

\item{sample_list, }{a \code{data.table}, the output of
\code{\link{read_metadt}}.}

\item{markers, }{a \code{data.table}, the output of
\code{\link{read_NGS_intervals}} or \code{\link{read_finalreport_snps}},
depending on the initial data type.}

\item{plots_path, }{path where the plots should be saved, if \code{NA} no plot is
produced.}

\item{...}{compatibility

  This function produce several summary statistics on the CNVs results in
  input. Ideally, it should be used interactively, together with
  \code{\link{cleaning_filter}}. Some information is printed on the console
  (mostly \code{summary()} on several characteristics of the results and the
  cohort), several plots can be produced and saved in the user specified
  location. Finally the function also return a \code{data.table} of sample-level
  summary.}
}
\value{
a \code{data.table} with summary statistics for the samples in the
  cohort.
}
\description{
Explore CNV calling results prior filtering
}
\examples{
DT <- summary(penn_22, cohort_examples, markers_examples)
}

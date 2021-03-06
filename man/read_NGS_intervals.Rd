% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_NGS_intervals.R
\name{read_NGS_intervals}
\alias{read_NGS_intervals}
\title{Read genomic intervals}
\usage{
read_NGS_intervals(
  DT_path,
  chr_col = "CONTIG",
  start_col = "START",
  end_col = "END"
)
}
\arguments{
\item{DT_path, }{path to the input file.}

\item{chr_col, }{name of the column containing the chromosome information in
the input file.}

\item{start_col, }{name of the column containing the start information in the
input file.}

\item{end_col, }{name of the column containing the end information in the input
file.}
}
\value{
a \code{data.table}, will be of \code{Markers} class in future
  versions.
}
\description{
\code{read_intervals} handles inputs of data used as the markers in CNVs
calling/segmentation using NGS data (WES or WGS)
}
\details{
This function is used to load data in interval list or BED like formats into
a \code{data.table} that integrates with the other functions of the package.
This is usually done at the beginning of a project involving CNVs
calling/segmentation on NGS data (WES or WGS) pipelines' results analysis.
The function should automatically skip any eventual header.
The parameters default values are for file in GATK interval list like format.
}
\examples{
read_NGS_intervals(DT_path = system.file("extdata", "markers_WES_example.tsv",
package = "CNVgears"), chr_col = "chr", start_col = "start", end_col = "end")
}

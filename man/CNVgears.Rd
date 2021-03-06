% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CNVgears-package.R
\docType{package}
\name{CNVgears}
\alias{CNVgears}
\title{CNVgears: A package to analyze CNVs calling/segmentation results}
\description{
\code{CNvgears} provides several functions to analyze the results of CNVs
calling and/or segmentation on SNPs arrays or NGS data.
}
\details{
The CNVgears package provides several functions useful in order to perform a
series of analysis the result of CNVs calling or segmentation pipelines or
algorithms, on both Ilummina SNP array (e.g. PennCNV, iPattern or
EnsembleCNV) and NGS data (e.g. ModSeg and gCNV pipelines from GATK), in an
integrated framework. To do so all the data is imported in a standardized
manner, allowing the user to perform analysis and data manipulation regardless
of the initial raw data type, from (among the others) CNVRs creation and
exclusion of immunoglobulin regions, to de novo CNVs discovery and genic
content annotation.

It has been originally developed for the CNVs characterization of the Italian
Autism Network (ITAN) collection (DOI: 10.1186/s12888-018-1937-y).
}
\section{Analysis pipelines examples}{

Here are briefly illustrated some workflow examples that can be done either
interactively on sequentially. See the vignettes for further details.

Staring from the results of gCNV and ModSeg pipelines on WES data in a cohort
of families:
\enumerate{
  \item load the intervals list (using \code{\link{read_NGS_intervals}});
  \item load samples table with minimal metadata (sample ID, sex, role, family ID);
  \item load the segmentation results of all the samples in the cohort, for each
  pipeline separately;
  \item merge adjacent segments (with equal CN);
  \item filter out CNVs in immunoglobulin (IG) regions;
  \item find eventual oversegmented samples (can be marked or excluded from the
  analysis);
  \item find replicated segments in the pipelines and merge the results into a single
  \code{data.table};
  \item create the CNVRs;
  \item exclude common CNVs based on the CNVRs frequency;
  \item annotate genic contents of the CNVs
  \item find the inheritance pattern of a selected subset of events (or the whole
  dataset) in the offspring, based on the segments of the parents;
  \item fine-screen putative de novo calls using the per-interval raw data (copy
  ratio or LRR like) of the trio;
  \item visualize the good de novo candidate per point raw data in the family to
  visually confirm the inheritance pattern.
}
}

\section{CNVgears functions}{

The CNVgears functions are organized in groups:
\itemize{
  \item input/output
  \item filtering
  \item CNVRs
  \item inter results comparison and merging
  \item de novo discovery/inheritance pattern detection
  \item plotting
}
}


# README #

CNVgears is an R package that implement a framework to integrate and analyze 
CNVs calling/segmentation results from multiple pipeline/algorithms and data
type (SNPs array and NGS). It is particularly suited for studies of family 
based cohorts. 

We provide functions to standardize several step of CNVs calling results 
analysis, including:

* Different methods/programs/pipelines results integration/merging;
* CNVs filtering;
* CNVs inheritance pattern detection;
* Copy Number Variable Regions (CNVRs) creation;
* Putative de novo calls visual confirmation;
* Basic CNVs annotation;
* Precessed data export in standard formats.

To install the latest version from R: 
`devtools::install_github("sinomem/cnvgears")`

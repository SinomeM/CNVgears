# README #

cnvgeaRs is an R package that implement a framework to analyze CNV calling/segmentation results from multiple 
pipeline/algorithms and data type. It is particularly suited for studies of family based cohorts. 

### Illustrative pipeline ###

The following points briefly illustrate how to use the main function of the package in the suggested order. 
Most of the step are completely optional, however some depends on previous steps. For instance CNVRs must be 
created if a user want to extract the rare call from his cohort's results.

* Load samples metadata from a PED file [read_metadt()];
* Load individual markers (SNPs/intervals depending on the raw data type) list for each methods [read_finalreport_snps() or read_NGS_intervals(), depending on data type];
* Load individual methods results for all samples [read_results()];
* Load marker-level raw data for all samples [read_finalreport_raw() or read_NGS_raw()];
* Perform intra-methods calls merging for each sample (i.e. merge eventual over-segmented calls) [merge_call()];
* Explore the results summary statistics in order to decide the best filtering parameters [summary_stat()];
* Filter the results based on several possible criteria, between the other, number of points (markers) per
  call, length, number of calling methods, telomeric and centromeric regions, immunoglobulin regions [cleaning_filter()];
* Filter eventual over-segmented samples [cleaning_filter()]; 
* Perform inter-methods calls merging for each sample (i.e. detect detect the calls made by more than one 
  method and combine them) in order to obtain a single results dataset [inter_res_merge()];
* Compute Copy Number Variable Regions (CNVRs) [cnvrs_create()]; 
* Annotate genomic locus of the calls and search for CNVs in known disease-linked loci (e.g. in studies on 
  ASD or schizophrenia) [genomic_locus()];
* Extract confident rare calls using CNVRs (and DGV) [not fully developed yet];
* Annotate genic content of selected CNVs [genic_load()];
* Detect CNVs inheritance pattern, in family based studies (this should be done per data type) [cnvs_inheritance()]; 
* Visually confirm good de novo CNVs candidate with LRR/BAF plot of the region of interest in the family [lrr_trio_plot()];
* Export the selected results as BED to be imported (for example) in IGV (integrative genome viewer) [not fully developed yet];

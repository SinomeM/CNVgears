
# Retrieve genomic information for the major assemblies
#
# Generates chromosome start, end centrosome location, as well as chromosomal
# arm start end datatsets in three different assemblies.
#
# The important related objects are already bundle in the package.
#
# internal function, does not need to be exported

chr_st_en_etc <- function() {

  # hg38
  bands <- fread("https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz")
  bands <- bands[V1 %in% paste0("chr", c(1:22, "X", "Y")),]
  merge(bands[,(.SD[1]),by = V1][, .(V1,V2)],
        bands[,(.SD[.N]),by = V1][, .(V1,V3)]) -> hg38_start_end
  colnames(hg38_start_end) <- c("chr", "start", "end")
  DT <- data.table()
  for (i in 1:nrow(bands)) {
    if (substr(bands$V4[i], 1,1) == "p" &
        substr(bands$V4[i+1], 1,1) == "q")
      DT <- rbind(DT, data.table("chr"= bands$V1[i], "centromere" = bands$V3[i]))
  }
  hg38_start_end_centromeres <- merge(hg38_start_end, DT)
  hg38_start_end_centromeres <- chr_uniform(hg38_start_end_centromeres)

  hg38_chr_arms <- data.table()
  for (i in 1:nrow(hg38_start_end_centromeres)) {
    hg38_chr_arms <- rbind(hg38_chr_arms, data.table("chr" = hg38_start_end_centromeres$chr[i],
                                                     "arm_ID" = paste0(hg38_start_end_centromeres$chr[i],
                                                                       c("p", "q")),
                                                     "start" = c(0, hg38_start_end_centromeres$centromere[i] + 1),
                                                     "end" = c(hg38_start_end_centromeres$centromere[i],
                                                               hg38_start_end_centromeres$end[i])))
  }
  rm(i)

  usethis::use_data(hg38_chr_arms)
  usethis::use_data(hg38_start_end_centromeres)
  rm(hg38_chr_arms)
  rm(hg38_start_end_centromeres)

  # hg19
  bands <- fread("https://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz")
  bands <- bands[V1 %in% paste0("chr", c(1:22, "X", "Y")),]
  merge(bands[,(.SD[1]),by = V1][, .(V1,V2)],
        bands[,(.SD[.N]),by = V1][, .(V1,V3)]) -> hg19_start_end
  colnames(hg19_start_end) <- c("chr", "start", "end")
  DT <- data.table()
  for (i in 1:nrow(bands)) {
    if (substr(bands$V4[i], 1,1) == "p" &
        substr(bands$V4[i+1], 1,1) == "q")
      DT <- rbind(DT, data.table("chr"= bands$V1[i], "centromere" = bands$V3[i]))
  }
  hg19_start_end_centromeres <- merge(hg19_start_end, DT)
  hg19_start_end_centromeres <- chr_uniform(hg19_start_end_centromeres)

  hg19_chr_arms <- data.table()
  for (i in 1:nrow(hg19_start_end_centromeres)) {
    hg19_chr_arms <- rbind(hg19_chr_arms, data.table("chr" = hg19_start_end_centromeres$chr[i],
                                                     "arm_ID" = paste0(hg19_start_end_centromeres$chr[i],
                                                                       c("p", "q")),
                                                     "start" = c(0, hg19_start_end_centromeres$centromere[i] + 1),
                                                     "end" = c(hg19_start_end_centromeres$centromere[i],
                                                               hg19_start_end_centromeres$end[i])))
  }
  rm(i)

  usethis::use_data(hg19_chr_arms)
  usethis::use_data(hg19_start_end_centromeres)
  rm(hg19_chr_arms)
  rm(hg19_start_end_centromeres)

  # hg18
  bands <- fread("https://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/cytoBand.txt.gz")
  bands <- bands[V1 %in% paste0("chr", c(1:22, "X", "Y")),]
  merge(bands[,(.SD[1]),by = V1][, .(V1,V2)],
        bands[,(.SD[.N]),by = V1][, .(V1,V3)]) -> hg18_start_end
  colnames(hg18_start_end) <- c("chr", "start", "end")
  DT <- data.table()
  for (i in 1:nrow(bands)) {
    if (substr(bands$V4[i], 1,1) == "p" &
        substr(bands$V4[i+1], 1,1) == "q")
      DT <- rbind(DT, data.table("chr"= bands$V1[i], "centromere" = bands$V3[i]))
  }
  hg18_start_end_centromeres <- merge(hg18_start_end, DT)
  hg18_start_end_centromeres <- chr_uniform(hg18_start_end_centromeres)

  hg18_chr_arms <- data.table()
  for (i in 1:nrow(hg18_start_end_centromeres)) {
    hg18_chr_arms <- rbind(hg18_chr_arms, data.table("chr" = hg18_start_end_centromeres$chr[i],
                                                     "arm_ID" = paste0(hg18_start_end_centromeres$chr[i],
                                                                       c("p", "q")),
                                                     "start" = c(0, hg18_start_end_centromeres$centromere[i] + 1),
                                                     "end" = c(hg18_start_end_centromeres$centromere[i],
                                                               hg18_start_end_centromeres$end[i])))
  }
  rm(i)

  usethis::use_data(hg18_chr_arms)
  usethis::use_data(hg18_start_end_centromeres)
  rm(hg18_chr_arms)
  rm(hg18_start_end_centromeres)
}

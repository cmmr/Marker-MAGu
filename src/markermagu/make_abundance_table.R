#!/usr/bin/env Rscript


suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))


args = commandArgs(trailingOnly=TRUE)

if (length(args)!=5) {
  stop("four arguments must be supplied \
       (*marker-magu.unique_alignment.coverm.tsv, read number, \
       out directory/name then sampleID).", call.=FALSE)
} else if (length(args)==5) {
  sprintf("arguments found. Running.")
}

## load coverm table
test_coverm_dt2 <- fread(args[1], sep = "\t", header = T, 
                         col.names = c("contig", "length", 
                                       "covered_bases", "reads_aligned"))

TOTAL_READS <- args[2]

TOTAL_READS <- as.numeric(TOTAL_READS)


## split contig field
test_coverm_dt_split2 <- test_coverm_dt2[, c("gene_name", "lineage", "genome") := 
                                           tstrsplit(contig, ";", fixed=TRUE)] 


if (args[5] == "relaxed"){

  summary_dt2 <- test_coverm_dt_split2 %>%
    group_by(lineage) %>%
    summarize(total_genes = n(),
              ## detected genes = genes with 1 or more aligned read
              detected_genes = sum(reads_aligned >= 1), 
              total_length = sum(length), 
              total_aligned_reads = sum(reads_aligned)) %>%
    ## sets threshold of detection at 33.3% genes detected
    filter(detected_genes/total_genes >= 0.333, 
           ## 3 or more detected genes (stringent with SGBs with few markers)
           detected_genes >= 3,
           ## 4 or more total marker genes for SGB
           total_genes >= 4, 
           ## 10 or more total aligned reads
           total_aligned_reads >= 10) %>%
    ## RPKM and relative abundance metric
    mutate(RPKM = 
             (
               (total_aligned_reads/(total_length/1000)) / (TOTAL_READS/1000000)
             ),
           rel_abundance = round(RPKM / sum(RPKM), 5)
    ) %>%
    ## only reporting SGBs with rel_abundance of at least 0.0001
    filter(rel_abundance >= 0.0001)
} else {
  summary_dt2 <- test_coverm_dt_split2 %>%
    group_by(lineage) %>%
    summarize(total_genes = n(),
              ## detected genes = genes with 1 or more aligned read
              detected_genes = sum(reads_aligned >= 1), 
              total_length = sum(length), 
              total_aligned_reads = sum(reads_aligned)) %>%
    ## sets threshold of detection at 75% genes detected
    filter(detected_genes/total_genes >= 0.75, 
           ## 4 or more total marker genes for SGB
           total_genes >= 4, 
           ## 10 or more total aligned reads
           total_aligned_reads >= 10) %>%
    ## RPKM and relative abundance metric
    mutate(RPKM = 
             (
               (total_aligned_reads/(total_length/1000)) / (TOTAL_READS/1000000)
             ),
           rel_abundance = round(RPKM / sum(RPKM), 5)
    ) %>%
    ## only reporting SGBs with rel_abundance of at least 0.0001
    filter(rel_abundance >= 0.0001)
}
summary_dt2$sampleID <- as.character(args[4])

write.table(summary_dt2, sprintf("%s.detected_species.tsv", args[3]), 
            quote = F, row.names = F, sep = "\t")





#!/usr/bin/env Rscript


suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))


args = commandArgs(trailingOnly=TRUE)

if (length(args)!=4) {
  stop("Two arguments must be supplied (*marker-magu.unique_alignment.coverm.tsv, read number, out directory/name then sampleID).", call.=FALSE)
} else if (length(args)==4) {
  # default output file
  sprintf("arguments found. Running.")
}

test_coverm_dt2 <- fread(args[1], sep = "\t", header = T, col.names = c("contig", "length", "covered_bases", "reads_aligned"))

TOTAL_READS <- args[2]

TOTAL_READS <- as.numeric(TOTAL_READS)


## split contig field
test_coverm_dt_split2 <- test_coverm_dt2[, c("gene_name", "lineage", "genome") := tstrsplit(contig, ";", fixed=TRUE)] 



summary_dt2 <- test_coverm_dt_split2 %>%
  group_by(lineage) %>%
  summarize(total_genes = n(), 
            detected_genes = sum(reads_aligned >= 1), 
            total_length = sum(length), 
            total_aligned_reads = sum(reads_aligned)) %>%
  filter(detected_genes/total_genes >= 0.75, total_genes >= 4, total_aligned_reads >= 10) %>%
  mutate(RPKM = 
           (
             (total_aligned_reads/(total_length/1000)) / (TOTAL_READS/1000000)
           ),
         rel_abundance = round(RPKM / sum(RPKM), 5)
  ) %>%
  filter(rel_abundance >= 0.0001)

summary_dt2$sampleID <- as.character(args[4])

write.table(summary_dt2, sprintf("%s.detected_species.tsv", args[3]), quote = F, row.names = F, sep = "\t")





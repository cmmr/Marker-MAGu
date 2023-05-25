#!/usr/bin/env Rscript

library(data.table)
library(dplyr)

## only argument is the directory where the *.detected_species.tsv files are located
args = commandArgs(trailingOnly=TRUE)

## list of relevant files in the directory
detection_files <- list.files(path = sprintf("%s", args[1]), 
                        pattern = "*.detected_species.tsv", full.names = TRUE)

## load and combine all as one table
combined_table <- rbindlist(lapply(detection_files, fread))

## write table to current working directory
write.table(combined_table, 
            file = sprintf("%s.combined_profile.tsv", args[1]),
            sep = "\t", row.names = F, quote = F)

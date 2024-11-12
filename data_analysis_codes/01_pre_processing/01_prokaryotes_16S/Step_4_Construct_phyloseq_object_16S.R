# Load Libraries----
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(tidyverse); packageVersion("tidyverse")

## Setup working dictionary first----
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/Pre_processing_16S")

## load theme for ggplot
theme_set(theme_bw())

#new script
asv_tab <- readRDS("ASV_tab_all_batches_bayesian_20240101.rds")
asv_tax <- readRDS("ASV_tax_tab_all_batches_bayesian_20240101.rds")

# Load metadata
sample.info <- read.csv("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/Metadata/Fish_pond_metadata_3.csv")
row.names(sample.info) <- sample.info[,"Sample"] #phyloseq must have sample names as the row names of the table
#should subset to remove extra sample  names column

# combine asv sequence table, asv tax table and metadata into a phyloseq object
ps1 <- phyloseq(otu_table(asv_tab, taxa_are_rows=TRUE), 
               sample_data(sample.info), 
               tax_table(asv_tax))
ps1 # 38250 taxa and 1066 samples

#add reads number into the metadata
sample_data(ps1)$Read_depth <- sample_sums(ps1)

ps <- ps1

#save RDS file of phyloseq object
saveRDS(ps, "phyloseq_all_samples_20240105.rds")


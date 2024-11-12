# 01/04/2024
# Sanjit Debnath
# This script is to construct the phyloseq object

# Load library
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(tidyverse); packageVersion("tidyverse")

## Setup working dictionary first----
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/pre_processing_18S/primer_removed")

# set seeds
set.seed(1234)

# read asv table and tax tables for 90-150bp size amplicon
asv_tab <- readRDS("asv_tab_18S_90-150bp_silva_20240331.rds")# can use this for silva and pr2
asv_tax.silva <- readRDS("asv_tax_tab_18S_90-150bp_silva_20240331.rds")
asv_tax.pr2 <- readRDS("asv_tax_table_18S_90_150bp_pr2_v5.0.0_20240331.rds")

# Read metadata
sample.info <- read.csv("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/Metadata/Fish_pond_metadata_6_18S_20240403.csv")
row.names(sample.info) <- sample.info[,"Sample"] #phyloseq must have sample names as the row names of the table
#should subset to remove extra sample  names column

# construct phyloseq object
ps.silva <- phyloseq(otu_table(asv_tab, taxa_are_rows=TRUE), 
                     sample_data(sample.info), 
                     tax_table(asv_tax.silva))
ps.silva #9345 taxa and 955 samples
# due to some reason couldn't add "pcr_neg2_2". after trying several times, gave up

# save the phyloseq object for downstream analusis
saveRDS(ps.silva, "phyloseq_18S_90-150bp_silva_20240403.rds")

# construct phyloseq object
ps.pr2 <- phyloseq(otu_table(asv_tab, taxa_are_rows=TRUE), 
                   sample_data(sample.info), 
                   tax_table(asv_tax.pr2))
ps.pr2 #9345 taxa and 955 samples

# save the phyloseq object for downstream analusis
saveRDS(ps.pr2, "phyloseq_18S_90_150bp_pr2_v5.0.0_20240403.rds")


# from these, I will consider phyloseq object ps.pr2, because if only classified with pr2 is also enough





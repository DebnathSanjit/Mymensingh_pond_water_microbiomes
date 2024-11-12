# This script is to check the sequencing profile after all types of filtering

# Load libraries
library(phyloseq); packageVersion("phyloseq")
library(tidyverse); packageVersion("tidyverse")
library(microViz); packageVersion("microViz")

# Prokaryotic----
# load phyloseq object
ps <- readRDS("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/statistical_analysis/phyloseq_metadata_6_v3_20240419.rds")
ps # 10523 taxa and 891 samples

## Overall 16s dataset----
ps.16s.pt <- ps %>% 
  ps_filter(
    Crop_species != "Shing",
    Crop_species != "Gulsha.Carp", 
    Crop_species != "Gulsha.Pabda"
  ) %>% tax_fix()
ps.16s.pt  # 10378 taxa and 812 samples 

# Check the total reads in the phyloseq object
total_reads.16s2 <- sum(otu_table(ps.16s.pt))
total_reads.16s2 # 10358751
read_per_sam.16s2 <- total_reads.16s2/812
read_per_sam.16s2 #12757.08

## check summary
summary(sample_sums(ps.16s.pt))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2001    6057   10135   12757   15077  156489 

## Count the number of phyla and genera----
### Convert phyloseq object to phylum level----
ps.phyla <- tax_glom(ps.16s.pt, taxrank = "Phylum")
ps.phyla # 54 taxa and 812 samples

### Convert phyloseq object to family level----
ps.family <- tax_glom(ps.16s.pt, taxrank = "Family")
ps.family # 368 taxa and 812 samples
# 609 taxa and 812 samples

### Convert phyloseq object to genus level----
ps.genus <- tax_glom(ps.16s.pt, taxrank = "Genus")
ps.genus # 611 taxa and 812 samples
# 1092 taxa and 812 samples 

## Check sequencing profile species wise for supplementary table----
## subset only pangasius pond samples####
pseq_pang <- ps.16s.pt %>% 
  ps_filter(
    Crop == "Pangasius"
  ) %>% tax_fix()
pseq_pang # 8091 taxa and 421 samples

### check the total reads in the phyloseq object----
total_reads.pseq_pang16s2 <- sum(otu_table(pseq_pang))
total_reads.pseq_pang16s2 # 4671008

## check summary
summary(sample_sums(pseq_pang))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2076    5368    8932   11095   13722  156489 

## Count the number of phyla and genera----
### Convert phyloseq object to phylum level----
pseq_pang.phyla <- tax_glom(pseq_pang, taxrank = "Phylum")
pseq_pang.phyla # 52 taxa and 421 samples

### Convert phyloseq object to family level----
pseq_pang.family <- tax_glom(pseq_pang, taxrank = "Family")
pseq_pang.family # 353 taxa and 421 samples
# 577 taxa and 421 samples

### Convert phyloseq object to genus level----
pseq_pang.genus <- tax_glom(pseq_pang, taxrank = "Genus")
pseq_pang.genus # 547 taxa and 421 samples
# 996 taxa and 421 samples

## subset only tilapia pond samples####
pseq_tila <- ps.16s.pt %>% 
  ps_filter(
    Crop == "Tilapia"
  ) %>% tax_fix()
pseq_tila # 8380 taxa and 391 samples

### check the total reads in the phyloseq object----
total_reads.pseq_tila16s2 <- sum(otu_table(pseq_tila))
total_reads.pseq_tila16s2 # 5687743

## check summary
summary(sample_sums(pseq_tila))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2001    6872   11005   14547   17814   83792 

### Convert phyloseq object to phylum level----
pseq_tila.phyla <- tax_glom(pseq_tila, taxrank = "Phylum")
pseq_tila.phyla # 52 taxa and 391 samples

### Convert phyloseq object to family level----
pseq_tila.family <- tax_glom(pseq_tila, taxrank = "Family")
pseq_tila.family # 341 taxa and 391 samples
#  566 taxa and 391 samples

### Convert phyloseq object to genus level----
pseq_tila.genus <- tax_glom(pseq_tila, taxrank = "Genus")
pseq_tila.genus # 557 taxa and 391 samples
# 1004 taxa and 391 samples 

# Microeukaryotic----
# load phyloseq object
ps.18s <- readRDS("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/statistical_analysis/phyloseq_18S_filtered_with_tree_pr2_90-150bp_20240416.rds")
ps.18s # 5390 taxa and 872 samples

## Overall 16s dataset----
ps.18s.pt <- ps.18s %>% 
  ps_filter(
    Crop_species != "Shing",
    Crop_species != "Gulsha.Carp", 
    Crop_species != "Gulsha.Pabda"
  ) %>% tax_fix()
ps.18s.pt  # 5210 taxa and 797 samples 

# Check the total reads in the phyloseq object
total_reads.18s2 <- sum(otu_table(ps.18s.pt))
total_reads.18s2 # 7676770
read_per_sam.18s2 <- total_reads.18s2/797
read_per_sam.18s2 # 9632.083

## check summary
summary(sample_sums(ps.18s.pt))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2002    5569    8520    9632   12261   76009 

## Count the number of phyla and genera----
### Convert phyloseq object to phylum level----
ps.18s.division <- tax_glom(ps.18s.pt, taxrank = "Division")
ps.18s.division # 19 taxa and 797 samples
#  26 taxa and 797 samples

### Convert phyloseq object to family level----
ps.18s.family <- tax_glom(ps.18s.pt, taxrank = "Family")
ps.18s.family # 226 taxa and 797 samples
# 314 taxa and 797 samples

### Convert phyloseq object to genus level----
ps.18s.genus <- tax_glom(ps.18s.pt, taxrank = "Genus")
ps.18s.genus # 360 taxa and 797 samples
#  534 taxa and 797 samples

## Check sequencing profile species wise for supplementary table----
## subset only pangasius pond samples####
pseq_pang.18s <- ps.18s.pt %>% 
  ps_filter(
    Crop == "Pangasius"
  ) %>% tax_fix()
pseq_pang.18s # 3947 taxa and 421 samples

# to check the total reads in the phyloseq object
total_reads.pseq_pang18s2 <- sum(otu_table(pseq_pang.18s))
total_reads.pseq_pang18s2 # 3806513

## check summary
summary(sample_sums(pseq_pang.18s))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2002    4870    8068    9042   11801   60811 

## Count the number of phyla and genera----
### Convert phyloseq object to phylum level----
pseq_pang.18s.division <- tax_glom(pseq_pang.18s, taxrank = "Division")
pseq_pang.18s.division #  18 taxa and 421 samples
# 24 taxa and 421 samples

### Convert phyloseq object to family level----
pseq_pang.18s.family <- tax_glom(pseq_pang.18s, taxrank = "Family")
pseq_pang.18s.family #  214 taxa and 421 samples
#  292 taxa and 421 samples 

### Convert phyloseq object to genus level----
pseq_pang.18s.genus <- tax_glom(pseq_pang.18s, taxrank = "Genus")
pseq_pang.18s.genus # 314 taxa and 421 samples
# 473 taxa and 421 samples

## subset only tilapia pond samples####
pseq_tila.18s <- ps.18s.pt %>% 
  ps_filter(
    Crop == "Tilapia"
  ) %>% tax_fix()
pseq_tila.18s # 4463 taxa and 376 samples

# to check the total reads in the phyloseq object
total_reads.pseq_tila18s2 <- sum(otu_table(pseq_tila.18s))
total_reads.pseq_tila18s2 # 3870257

## check summary
summary(sample_sums(pseq_tila.18s))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2067    6136    9114   10293   12819   76009 

### Convert phyloseq object to phylum level----
pseq_tila.18s.division <- tax_glom(pseq_tila.18s, taxrank = "Division")
pseq_tila.18s.division #  19 taxa and 376 samples
# 26 taxa and 376 samples

### Convert phyloseq object to family level----
pseq_tila.18s.family <- tax_glom(pseq_tila.18s, taxrank = "Family")
pseq_tila.18s.family # 219 taxa and 376 samples
#  305 taxa and 376 samples

### Convert phyloseq object to genus level----
pseq_tila.18s.genus <- tax_glom(pseq_tila.18s, taxrank = "Genus")
pseq_tila.18s.genus # 338 taxa and 376 samples
# 498 taxa and 376 samples 

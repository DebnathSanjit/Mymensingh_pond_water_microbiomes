# Date: 20241006
# Author: Sanjit Debnath
# PERMANOVA Model for effect of multiple variable on overall microbial community (prokaryotes), Factor with biggest R2 value has the biggest effect on overall microbial community

# Load Libraries
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(microViz); packageVersion("microViz")
library(vegan); packageVersion("vegan") # needed for PERMANOVA test
library(tidyverse); packageVersion("tidyverse")
library(pairwiseAdonis); packageVersion("pairwiseAdonis")
library(writexl); packageVersion("writexl")
# Libraries for tests
#library(lme4)
#library(lmerTest)
#library(variancePartition);packageVersion("variancePartition")

## Setup working dictionary first
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/statistical_analysis")
getwd ()

# Set theme for ggplot
theme_set(theme_bw())
set.seed(1234)

# Prokaryotes
# Load the phyloseq object----
ps <- readRDS("phyloseq_metadata_6_v3_20240419.rds")
ps #10523 taxa and 891 samples

# Add season variables as i require
ps <- ps %>% 
  ps_mutate(Season2 = case_when(
    Season %in% c("01.Monsoon", "04.Monsoon") ~ "Monsoon",
    Season %in% c("02.Winter", "05.Winter") ~ "Winter",
    Season %in% c("03.Pre.Monsoon", "06.Pre.monsoon") ~ "Pre-monsoon",
    TRUE ~ as.character(Season)),
    Season3 =  case_when(
      Season %in% "01.Monsoon" ~ "Monsoon1",
      Season %in% "04.Monsoon" ~ "Monsoon2",
      Season %in% "02.Winter" ~ "Winter1",
      Season %in% "05.Winter" ~ "Winter2",
      Season %in% "03.Pre.Monsoon" ~ "Pre_monsoon1",
      Season %in% "06.Pre.monsoon" ~ "Pre_monsoon2",
      TRUE ~ as.character(Season))) 
ps # 10523 taxa and 891 samples

# Subset in such a way that this object contains all the samples covering all the variables
pseq2 <- ps %>% 
  ps_filter(
    Crop_species != "Shing" & Crop_species != "Gulsha.Carp" & Crop_species != "Gulsha.Pabda",
    Culture_system != "P.no.fish",
    Sampling_point > 9) %>% 
  ps_mutate(Culture_system2 = ifelse(grepl("\\.", Crop_species), "Polyculture", "Monoculture")) %>%  
  ps_mutate(Culture_system3 = ifelse(grepl("\\.", Crop_species), "01. Polyculture", Crop_species)) 
pseq2 # 8922 taxa and 418 samples

# Prepare metadata for permanova
metadata_all <- sample_data(pseq2) %>%
  data.frame() #%>%
#  tibble()

# Check for any missing values
metadata_no_na <- na.omit(metadata_all) # Temperature_C has one missing value
missing_rows <- which(is.na(metadata_all$Temperature_C))
print(missing_rows)

# Specify the value you want to add
specific_value <- 30.70  # The value you want to insert
row_index <- 15  # The row number where the value is missing

# Method 2: Using base R
metadata_all[row_index, "Temperature_C"] <- specific_value

# Calculate distance matrix
wuni_dist_all = phyloseq::distance(pseq2, method = "wunifrac")

# PERMANOVA----
PERM_all_wuni <- adonis2(wuni_dist_all ~ Season3 + Upazila + Reported_disease + Culture_system +
                             Temperature_C + Salinity_ppm + 
                            pH + DO_mg_l, data = metadata_all)
PERM_all_wuni

# Convert PERMANOVA result to a data frame
permanova_df <- as.data.frame(PERM_all_wuni)
permanova_df <- cbind(Factor = rownames(permanova_df), permanova_df)
rownames(permanova_df) <- NULL  # Remove row names
# Save to Excel
#write_xlsx(permanova_df, "PERMANOVA_results_all_factor_20241006.xlsx")

##### Pairwise adonis----
pair.season <- pairwise.adonis(wuni_dist_all, metadata_all$Season3,  p.adjust.m = "BH")
pair.up <- pairwise.adonis(wuni_dist_all, metadata_all$Upazila,  p.adjust.m = "BH")
pair.dis <- pairwise.adonis(wuni_dist_all, metadata_all$Reported_disease,  p.adjust.m = "BH")
pair.cs <- pairwise.adonis(wuni_dist_all, metadata_all$Culture_system,  p.adjust.m = "BH")

# Convert each pairwise results object to a data frame and add row names as a column
# Assume pairwise results contain a significance (e.g., p-value) column
season_df <- as.data.frame(pair.season)
upazila_df <- as.data.frame(pair.up)
disease_df <- as.data.frame(pair.dis)
culture_system_df <- as.data.frame(pair.cs)

# If the significance column (e.g., p-value) is missing, make sure to add it here
# Add row names as the "Comparison" column to each data frame for easier reference
season_df <- cbind(Comparison = rownames(season_df), season_df)
upazila_df <- cbind(Comparison = rownames(upazila_df), upazila_df)
disease_df <- cbind(Comparison = rownames(disease_df), disease_df)
culture_system_df <- cbind(Comparison = rownames(culture_system_df), culture_system_df)

# Remove row names to keep the tables clean
rownames(upazila_df) <- NULL
rownames(season_df) <- NULL
rownames(disease_df) <- NULL
rownames(culture_system_df) <- NULL

# Organize the data frames into a list, where each item will become a sheet in Excel
pairwise_results_list <- list(
  "Upazila" = upazila_df,
  "Season" = season_df,
  "Disease" = disease_df,
  "Culture System" = culture_system_df
)

# Export to an Excel workbook with significance values included
#write_xlsx(pairwise_results_list, "PairwiseAdonis_with_significance_for_all_20241007.xlsx")









# 11/05/2024
# Sanjit Debnath
# This script is to plot rarefaction curve for prokaryotes and microeukaryotes dataset

# Load libraries
library(phyloseq); packageVersion("phyloseq")
library(tidyverse); packageVersion("tidyverse")
library(microViz); packageVersion("microViz")
library(ggplot2); packageVersion("ggplot2")
library(microbiome); packageVersion("microbiome")# data analysis and visualisation (here used for rarefaction curve)
library(cowplot); packageVersion("cowplot")

# set seed
set.seed(1234)

# Prokaryotes----
# load phyloseq object
ps <- readRDS("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/statistical_analysis/phyloseq_metadata_6_v3_20240419.rds")
ps # 10523 taxa and 891 samples

# Subser only pangasius and tilapia
ps.16s.pt <- ps %>% 
  ps_filter(
    Crop_species != "Shing",
    Crop_species != "Gulsha.Carp", 
    Crop_species != "Gulsha.Pabda"
  )
ps.16s.pt  # 10378 taxa and 812 samples 

## Add reads number into the metadata
#sample_data(ps)$Read_depth <- sample_sums(ps)

## Make a barplot of reads----
# Extract sample data from Phyloseq object
sample_data.16s <- sample_data(ps.16s.pt)

# Extract Sample and Read_depth columns
sample_read_depth.16s <- sample_data.16s[, c("Sample", "Read_depth")]

# Convert to dataframe
sample_read_depth_df.16s <- as.data.frame(sample_read_depth.16s)

# Calculate the mean of the reads
mean_reads.16s <- mean(sample_read_depth_df.16s$Read_depth)

# Create bar plot
readplot16s <- ggplot(sample_read_depth_df.16s, aes(x = Sample, y = Read_depth)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Sample", 
       y = "Read Depth", 
       title = "Read Depth of prokaryotes samples") +
  theme(
    axis.text.x = element_blank(), # to remove sample names from x-axis
    #axis.text.x = element_text(size = 12),  # Adjust font size for x-axis labels
    axis.text.y = element_text(size = 14),  # Adjust font size for y-axis labels
    axis.title.x = element_text(size = 18),  # Adjust font size for x-axis title
    axis.title.y = element_text(size = 18),  # Adjust font size for y-axis title
    plot.title = element_text(size = 18),   # Adjust font size for plot title
    panel.grid.major = element_blank(),     # Remove major gridlines
    panel.grid.minor = element_blank()) +   # Remove minor gridlines
  ylim(0, 170000) +
  geom_hline(yintercept = mean_reads.16s, linetype = "dotted", color = "red", linewidth = 1)  # Add a dotted line representing the mean reads
readplot16s

## Make a rarefraction curve----
# follow this: https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/alpha-diversities.html
# from the summary above, we can see there is a large difference (over 37*)in the number of reads
# by plotting rarefaction curve we can check how has the richness captured in the sequencing effort

# Rarefraction curve
otu_tab.16s <- t(abundances(ps.16s.pt))  # Transpose the abundance table
rcurve.16s <- vegan::rarecurve(otu_tab.16s, 
                                step = 50, label = FALSE, 
                                sample = min(rowSums(otu_tab.16s), 
                                             col = "blue", cex = 0.6)) # 11:32 - 11:37, takes around 5 minutes
# This one doesn't change the color but gives the vertical line

# Microeukaryotes----
# load phyloseq object
ps.18s <- readRDS("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/statistical_analysis/phyloseq_18S_filtered_with_tree_pr2_90-150bp_20240416.rds")
ps.18s # 5390 taxa and 872 samples

# Subset pangasius and tilapia
ps.18s.pt <- ps.18s %>% 
  ps_filter(
    Crop_species != "Shing",
    Crop_species != "Gulsha.Carp", 
    Crop_species != "Gulsha.Pabda"
  )
ps.18s.pt  # 5210 taxa and 797 samples 

## Add reads number into the metadata
#sample_data(ps.18s)$count_per_sample <- sample_sums(ps.18s) 

## Make a barplot of reads----
# Extract sample data from Phyloseq object
sample_data.18s <- sample_data(ps.18s.pt)

# Extract Sample and count_per_sample columns
sample_read_depth.18s <- sample_data.18s[, c("Sample", "count_per_sample")]

# Convert to dataframe
sample_read_depth_df.18s <- as.data.frame(sample_read_depth.18s)

# Calculate the mean of the reads
mean_reads.18s <- mean(sample_read_depth_df.18s$count_per_sample)

# Create the bar plot
readplot18s <- ggplot(sample_read_depth_df.18s, aes(x = Sample, y = count_per_sample)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Sample", 
       y = "Read Depth", 
       title = "Read Depth of microeukaryotes samples") +
  theme(
    axis.text.x = element_blank(), # to remove sample names from x-axis
    #axis.text.x = element_text(size = 12),  # Adjust font size for x-axis labels
    axis.text.y = element_text(size = 14),  # Adjust font size for y-axis labels
    axis.title.x = element_text(size = 18),  # Adjust font size for x-axis title
    axis.title.y = element_text(size = 18),  # Adjust font size for y-axis title
    plot.title = element_text(size = 18),   # Adjust font size for plot title
    panel.grid.major = element_blank(),     # Remove major gridlines
    panel.grid.minor = element_blank()) +   # Remove minor gridlines
  ylim(0, 170000) +
  geom_hline(yintercept = mean_reads.18s, linetype = "dotted", color = "red", linewidth = 1)  # Add a dotted line representing the mean reads
readplot18s

## Make a rarefraction curve----
otu_tab.18s <- t(abundances(ps.18s.pt)) # Transpose the abundance table
rcurve.18s <- vegan::rarecurve(otu_tab.18s, 
                                step = 50, label = FALSE, 
                                sample = min(rowSums(otu_tab.18s), 
                                             cex = 0.6, col = "blue")) # 11:39 - 11:40, not take much time

# Combine plot----
# Barplot
comb.barplot <- cowplot::plot_grid(readplot16s,
                                   readplot18s,
                                   labels = "AUTO")
comb.barplot # Not using this plot, since I have these information in the supplementary table of sequence profile

# Rarefaction curve
# The rarefaction curve generated with microbiome package is not a ggplot. So, this can not be converted using grid.arrange from gridExtra package, or as ggplot from ggplotify package
# So save both plots and combine using inkscape

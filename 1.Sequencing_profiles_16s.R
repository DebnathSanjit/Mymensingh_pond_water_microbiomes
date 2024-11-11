# Load libraries----
# Install and load the required packages
#if (!requireNamespace("ShortRead", quietly = TRUE)) {
#  install.packages("ShortRead")
#}

library(ShortRead)

# Set the working directory
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique")

# set seed
set.seed(1234)

## Step 1: Check counts before filtering----
# Show the path where the sequencing files are
path <- paste0(getwd(), "/Sequences/16S/B1")

# List all files in the specified directory with the pattern ".fq.gz"
files <- list.files(path, pattern = "r1.fq.gz", full.names = TRUE)

# Initialize a variable to store the total number of reads
total_reads <- 0

# Loop through each file, read the FASTQ, and count the number of reads
for (file in files) {
  reads <- readFastq(file)
  if (length(reads) > 0) {
    total_reads <- total_reads + length(reads)
  } else {
    warning(paste("No reads found in file:", file))
  }
}

# Print the total number of reads
print(total_reads)
B1_reads.16s <- 1786628

# Show the path where the sequencing files are
path <- paste0(getwd(), "/Sequences/16S/B2")

# List all files in the specified directory with the pattern ".fq.gz"
files <- list.files(path, pattern = "r1.fq.gz", full.names = TRUE)

# Initialize a variable to store the total number of reads
total_reads <- 0

# Loop through each file, read the FASTQ, and count the number of reads
for (file in files) {
  reads <- readFastq(file)
  if (length(reads) > 0) {
    total_reads <- total_reads + length(reads)
  } else {
    warning(paste("No reads found in file:", file))
  }
}

# Print the total number of reads
print(total_reads)
B2_reads.16s <- 1378463

# Show the path where the sequencing files are
path <- paste0(getwd(), "/Sequences/16S/B6")

# List all files in the specified directory with the pattern ".fq.gz"
files <- list.files(path, pattern = "r1.fq.gz", full.names = TRUE)

# Initialize a variable to store the total number of reads
total_reads <- 0

# Loop through each file, read the FASTQ, and count the number of reads
for (file in files) {
  reads <- readFastq(file)
  if (length(reads) > 0) {
    total_reads <- total_reads + length(reads)
  } else {
    warning(paste("No reads found in file:", file))
  }
}

# Print the total number of reads
print(total_reads)
B6_reads.16s <- 2509465

# Show the path where the sequencing files are
path <- paste0(getwd(), "/Sequences/16S/B7")

# List all files in the specified directory with the pattern ".fq.gz"
files <- list.files(path, pattern = "r1.fq.gz", full.names = TRUE)

# Initialize a variable to store the total number of reads
total_reads <- 0

# Loop through each file, read the FASTQ, and count the number of reads
for (file in files) {
  reads <- readFastq(file)
  if (length(reads) > 0) {
    total_reads <- total_reads + length(reads)
  } else {
    warning(paste("No reads found in file:", file))
  }
}

# Print the total number of reads
print(total_reads)
B7_reads.16s <- 3257161

# Show the path where the sequencing files are
path <- paste0(getwd(), "/Sequences/16S/B8")

# List all files in the specified directory with the pattern ".fq.gz"
files <- list.files(path, pattern = "r1.fq.gz", full.names = TRUE)

# Initialize a variable to store the total number of reads
total_reads <- 0

# Loop through each file, read the FASTQ, and count the number of reads
for (file in files) {
  reads <- readFastq(file)
  if (length(reads) > 0) {
    total_reads <- total_reads + length(reads)
  } else {
    warning(paste("No reads found in file:", file))
  }
}

# Print the total number of reads
print(total_reads)
B8_reads.16s <- 3716925

# Show the path where the sequencing files are
path <- paste0(getwd(), "/Sequences/16S/B9")

# List all files in the specified directory with the pattern ".fq.gz"
files <- list.files(path, pattern = "r1.fq.gz", full.names = TRUE)

# Initialize a variable to store the total number of reads
total_reads <- 0

# Loop through each file, read the FASTQ, and count the number of reads
for (file in files) {
  reads <- readFastq(file)
  if (length(reads) > 0) {
    total_reads <- total_reads + length(reads)
  } else {
    warning(paste("No reads found in file:", file))
  }
}

# Print the total number of reads
print(total_reads)
B9_reads.16s <- 4327541

# count reads without batch 4 as it contains samples with incubation, which I haven't considered
total_reads.16s <- B1_reads.16s + B2_reads.16s + B6_reads.16s + B7_reads.16s + B8_reads.16s + B9_reads.16s
total_reads.16s #16976183
total_reads.16s <- 16976183

# number of total samples
B1.16s <- (194/2)
B1.16s #97
B2.16s <- (154/2)
B2.16s #77
B6.16s <- (254/2)
B6.16s #127
B7.16s <- (350/2)
B7.16s #175
B8.16s <- (514/2)
B8.16s #257
B9.16s <- (460/2)
B9.16s #230
totalsample.16s <- B1.16s+B2.16s+B6.16s+B7.16s+B8.16s+B9.16s
totalsample.16s # 963
totalsample.16s <- 963

# check reads per sample before filtering
read_per_sam.16s <- total_reads.16s/totalsample.16s
read_per_sam.16s #17628.44

## Step 2: Check counts after filtering----
# after all the filtering, removing contaminants, we need to check the read numbers and remove samples with low reads
# to do that, following codes can be used
# To see the number of reads of each sample
library(phyloseq); packageVersion("phyloseq")
library(tidyverse); packageVersion("tidyverse")
library(microViz); packageVersion("microViz")

# load phyloseq object
ps <- readRDS("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/statistical_analysis/phyloseq_metadata_6_v3_20240419.rds")
ps # 10523 taxa and 891 samples

## easiest way to do it is -
##add reads number into the metadata
#sample_data(ps)$Read_depth <- sample_sums(ps)

## to check the total reads in the phyloseq object
#total_reads <- sum(otu_table(ps))
#total_reads # 11425934
#avg_read_sample <- total_reads/981
#avg_read_sample # 11647

# as I will not include those samples without pangasius or tilapia, remove them and check again
ps.16s.pt <- ps %>% 
  ps_filter(
    Crop_species != "Shing",
    Crop_species != "Gulsha.Carp", 
    Crop_species != "Gulsha.Pabda"
  )
ps.16s.pt  # 10378 taxa and 812 samples 

# to check the total reads in the phyloseq object
total_reads.16s2 <- sum(otu_table(ps.16s.pt))
total_reads.16s2 # 10358751
read_per_sam.16s2 <- total_reads.16s2/812
read_per_sam.16s2 #12757.08

## check summary
summary(sample_sums(ps.16s.pt))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2001    6057   10135   12757   15077  156489 

## Step 3: make a barplot of reads----
# Extract sample data from Phyloseq object
sample_data.16s <- sample_data(ps.16s.pt)

# Extract Sample and Read_depth columns
sample_read_depth.16s <- sample_data.16s[, c("Sample", "Read_depth")]

# Convert to dataframe
sample_read_depth_df.16s <- as.data.frame(sample_read_depth.16s)

# Calculate the mean of the reads
mean_reads.16s <- mean(sample_read_depth_df.16s$Read_depth)

# load library
library(ggplot2); packageVersion("ggplot2")

# Create bar plot
readplot16s <- ggplot(sample_read_depth_df.16s, aes(x = Sample, y = Read_depth)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Sample", 
       y = "Read Depth", 
       title = "a) Read Depth of Samples for prokaryotes (16S)") +
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
readplot16s

# Add text label for the average value
#text(x = 1, y = mean_reads.18s, labels = round(mean_reads.18s, digits = 2), pos = 3, col = "red")  # Adjust x and y coordinates as needed

# save with ggsave
ggsave("bar_plot.png", plot = readplot16s, width = 8, height = 6, dpi = 300)
ggsave("bar_plot.svg", plot = readplot16s, width = 8, height = 6, dpi = 300)

# combine 16s and 18s barplot
combined.readplot <- cowplot:::plot_grid(readplot16s, readplot18s, nrow = 2)
combined.readplot

# save with ggsave
ggsave("combined_bar_plot.png", plot = combined.readplot, width = 8, height = 6, dpi = 300)
ggsave("combined_bar_plot.svg", plot = combined.readplot, width = 8, height = 6, dpi = 300)

## Step 4: make a rarefraction curve----
# follow this: https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/alpha-diversities.html
# from the summary above, we can see there is a large difference (over 37*)in the number of reads
# by plotting rarefaction curve we can check how has the richness captured in the sequencing effort
# load library
library(microbiome); packageVersion("microbiome")# data analysis and visualisation

# rarefraction curve
otu_tab.16s <- t(abundances(ps.16s.pt))
rcurve.16s <- vegan::rarecurve(otu_tab.16s, 
                               step = 50, label = FALSE, 
                               #sample = min(rowSums(otu_tab), # for this, "col ="  doesn't work, but without this vertical line doesn't appear
                               col = "skyblue", cex = 0.6)

# this one doesn't change the color but gives the vertical line
rcurve.16s2 <- vegan::rarecurve(otu_tab.16s, 
                                step = 50, label = FALSE, 
                                sample = min(rowSums(otu_tab.16s), 
                                col = "blue", cex = 0.6))

# this one change the color but doesn't gives the vertical line
rcurve.16s3 <- vegan::rarecurve(otu_tab.16s, 
                                step = 50, label = FALSE, 
                                #sample = min(rowSums(otu_tab.16s), 
                                col = "blue", cex = 0.6)

# this can not be converted using grid.arrange from gridExtra package, or as.ggplot from ggplotify package
# so just save it as it is, then combine using powerpoint or inkscape

## Step 5: count the number of phyla and genera----
### Convert phyloseq object to phylum level----
ps.phyla <- tax_glom(ps.16s.pt, taxrank = "Phylum")
ps.phyla # 54 taxa and 812 samples

### Convert phyloseq object to family level----
ps.family <- tax_glom(ps.16s.pt, taxrank = "Family")
ps.family # 368 taxa and 812 samples

### Convert phyloseq object to genus level----
ps.genus <- tax_glom(ps.16s.pt, taxrank = "Genus")
ps.genus # 611 taxa and 812 samples

## Step 6: check sequencing profile species wise for supplementary table----
# subset only pangasius pond samples
pseq_pang <- ps.16s.pt %>% 
  ps_filter(
    Crop == "Pangasius"
  )
pseq_pang # 8091 taxa and 421 samples

### check the total reads in the phyloseq object----
total_reads.pseq_pang16s2 <- sum(otu_table(pseq_pang))
total_reads.pseq_pang16s2 # 4671008

## check summary
summary(sample_sums(pseq_pang))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2076    5368    8932   11095   13722  156489 

# subset only tilapia pond samples
pseq_tila <- ps.16s.pt %>% 
  ps_filter(
    Crop == "Tilapia"
  )
pseq_tila # 8380 taxa and 391 samples

### check the total reads in the phyloseq object----
total_reads.pseq_tila16s2 <- sum(otu_table(pseq_tila))
total_reads.pseq_tila16s2 # 5687743

## check summary
summary(sample_sums(pseq_tila))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2001    6872   11005   14547   17814   83792 


## Step 7: count the number of phyla and genera----
#### Pangasius
### Convert phyloseq object to phylum level----
pseq_pang.phyla <- tax_glom(pseq_pang, taxrank = "Phylum")
pseq_pang.phyla #  52 taxa and 421 samples

### Convert phyloseq object to family level----
pseq_pang.family <- tax_glom(pseq_pang, taxrank = "Family")
pseq_pang.family # 353 taxa and 421 samples

### Convert phyloseq object to genus level----
pseq_pang.genus <- tax_glom(pseq_pang, taxrank = "Genus")
pseq_pang.genus # 547 taxa and 421 samples

#### Tilapia
### Convert phyloseq object to phylum level----
pseq_tila.phyla <- tax_glom(pseq_tila, taxrank = "Phylum")
pseq_tila.phyla #  52 taxa and 391 samples

### Convert phyloseq object to family level----
pseq_tila.family <- tax_glom(pseq_tila, taxrank = "Family")
pseq_tila.family # 341 taxa and 391 samples

### Convert phyloseq object to genus level----
pseq_tila.genus <- tax_glom(pseq_tila, taxrank = "Genus")
pseq_tila.genus # 557 taxa and 391 samples

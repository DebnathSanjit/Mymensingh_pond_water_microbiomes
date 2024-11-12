# 12/09/2024 (updated from 23/04/2024)

# This script calculated the relative abundance, standard error, mean of standard errof of taxa
# at different taxonomic level. then also pairwise compare the taxa 

# Load Libraries
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(microViz); packageVersion("microViz")
library(tidyverse); packageVersion("tidyverse")
library(writexl); packageVersion("writexl")
# Libraries for tests
library(ggpubr); packageVersion("ggpubr")
library(rstatix); packageVersion("rstatix")
library(dunn.test); packageVersion("dunn.test")

## Setup working dictionary first
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/statistical_analysis")
getwd ()

# Prokaryotes----
# Load the phyloseq object
ps <- readRDS("phyloseq_metadata_6_v3_20240419.rds")
ps #10523 taxa and 891 samples

# Disease----
pseq_pang <- ps %>%
  ps_filter(
    Crop == "Pangasius",
    Crop_species != "Shing") %>% 
  tax_fix() 
pseq_pang # 8091 taxa and 421 samples

## Phylum---- 
# Aggregate at phyla level
dis.phyla <- tax_glom(pseq_pang, taxrank = "Phylum")
dis.phyla # 54 taxa and 421 samples

# calculate relative abundance
dis.phy.rel <- transform_sample_counts(dis.phyla, function(x) x / sum(x))# do not do this if want to see absulate abundance

# Convert to tibble
tb.dis <- psmelt(dis.phy.rel) %>%
  as_tibble()

### Mean and standard deviation----
tb1.dis <- tb.dis %>%
  group_by(Phylum, Reported_disease) %>%
  summarize(Mean = mean(Abundance), SD = sd(Abundance)) %>%
  ungroup()

### Standard error of the mean (SEM)----
tb2.dis <- tb1.dis %>%
  group_by(Phylum) %>%
  mutate(SEM = SD / sqrt(n())) %>%
  ungroup()

### Convert mean, SD, and SEM to percentages----
tb3.dis <- tb2.dis %>%
  mutate(
    Mean_percent = Mean * 100,
    SD_percent = SD * 100,
    SEM_percent = SEM * 100
  ) %>%
  select(Phylum, Reported_disease, Mean_percent, SD_percent, SEM_percent)

### Pairwise comparison (Kruskal-Wallis)of taxa----
# To perform the KW followed by dunn test, follow code from jamie, as the code form chatgpt seems giving very high p value
# Agglomerate at phylum level
V4.phy.dis <- pseq_pang %>% tax_glom("Phylum", NArm = FALSE) 

# Convert to relative abundance, (could use tidyPS(., TRUE) function from above)
V4.phy.ra.dis <- transform_sample_counts(V4.phy.dis, function(otu) otu/sum(otu))

# Keep only those grodiss with a mean relative abundance greater than 0.01%
V4.main.phy.dis <- filter_taxa(V4.phy.ra.dis, function(x) mean(x) > 0.0001, TRUE)

# Melt phyloseq object into a single tibble
V4.phy.melted.dis <- as_tibble(psmelt(V4.main.phy.dis), rownames = "Rows") %>%
  mutate(Abundance_adj = ifelse(Abundance == 0, NA, Abundance))   # Replace abundance zeros with NA

# Perform Kruskal-Wallis test for each phylum using the combined variable
v4.phy.kw.dis <- V4.phy.melted.dis %>%
  group_by(Phylum) %>%
  kruskal_test(Abundance_adj ~ Reported_disease)

# Perform Dunn's post-hoc test for each phylum where Kruskal-Wallis was significant
v4.phy.dunn.dis <- V4.phy.melted.dis %>%
  group_by(Phylum) %>% #also grodis by sample type to
  dunn_test(Abundance_adj ~ Reported_disease, p.adjust.method = "BH") #%>% 
#dplyr::mutate(FDR = p.adjust(.$p, method = "BH", 30)) #choosing number of comparisons as 30 as we have 10 phyla and three sample types from which comparisons are made
# p adjustment is needed based on how many times we are running the test.

### Combine relative abundance and test results together----
# Reshape tb3 to wide format, separating by Crop (tilaasius and Tilapia)
tb3_wide.dis <- tb3.dis %>%
  pivot_wider(
    names_from = Reported_disease,
    values_from = c(Mean_percent, SD_percent, SEM_percent)
  )

# Merge tb3_wide with v4.phy.dunn by Phylum
v4.phy.combined_data.dis <- v4.phy.dunn.dis %>%
  left_join(tb3_wide.dis, by = "Phylum")

# Save the combined data as an Excel file
#writexl::write_xlsx(v4.phy.combined_data.dis, "combined RA and KW results for disease_20241001.xlsx")
# also save tb3 for RA results
#writexl::write_xlsx(tb3.up, "RA of phyla disease_20241001.xlsx")

## Family----
# Aggregate at phyla level
dis.fam <- tax_glom(pseq_pang, taxrank = "Family")
dis.fam # 577 taxa and 421 samples

# calculate relative abundance
dis.fam.rel <- transform_sample_counts(dis.fam, function(x) x / sum(x))# do not do this if want to see absulate abundance

# Convert to tibble
tb.dis.fam <- psmelt(dis.fam.rel) %>%
  as_tibble()

### Mean and standard deviation----
tb1.dis.fam <- tb.dis.fam %>%
  group_by(Family, Reported_disease) %>%
  summarize(Mean = mean(Abundance), SD = sd(Abundance)) %>%
  ungroup()

### Standard error of the mean (SEM)----
tb2.dis.fam <- tb1.dis.fam %>%
  group_by(Family) %>%
  mutate(SEM = SD / sqrt(n())) %>%
  ungroup()

### Convert mean, SD, and SEM to percentages----
tb3.dis.fam <- tb2.dis.fam %>%
  mutate(
    Mean_percent = Mean * 100,
    SD_percent = SD * 100,
    SEM_percent = SEM * 100
  ) %>%
  select(Family, Reported_disease, Mean_percent, SD_percent, SEM_percent)

### Pairwise comparison (Kruskal-Wallis)of taxa----
# To perform the KW followed by dunn test, follow code from jamie, as the code form chatgpt seems giving very high p value
## Agglomerate at Family level
#V4.fam.dis <- pseq_pang %>% tax_glom("Family", NArm = FALSE) 
#
## Convert to relative abundance, (could use tidyPS(., TRUE) function from above)
#V4.fam.ra.dis <- transform_sample_counts(V4.fam.dis, function(otu) otu/sum(otu))
#
# Keep only those grodiss with a mean relative abundance greater than 0.1%
V4.main.fam.dis <- filter_taxa(dis.fam.rel, function(x) mean(x) > 0.001, TRUE)

# Melt phyloseq object into a single tibble
V4.fam.melted.dis <- as_tibble(psmelt(V4.main.fam.dis), rownames = "Rows") %>%
  mutate(Abundance_adj = ifelse(Abundance == 0, NA, Abundance))   # Replace abundance zeros with NA

# Perform Kruskal-Wallis test for each Family using the combined variable
v4.fam.kw.dis <- V4.fam.melted.dis %>%
  group_by(Family) %>%
  kruskal_test(Abundance_adj ~ Reported_disease)

# Perform Dunn's post-hoc test for each Family where Kruskal-Wallis was significant
v4.fam.dunn.dis <- V4.fam.melted.dis %>%
  group_by(Family) %>% #also grodis by sample type to
  dunn_test(Abundance_adj ~ Reported_disease, p.adjust.method = "BH") #%>% 
#dplyr::mutate(FDR = p.adjust(.$p, method = "BH", 30)) #choosing number of comparisons as 30 as we have 10 phyla and three sample types from which comparisons are made
# p adjustment is needed based on how many times we are running the test.

### Combine relative abundance and test results together----
# Reshape tb3 to wide format, separating by Crop (tilaasius and Tilapia)
tb3_wide.dis.fam <- tb3.dis.fam %>%
  pivot_wider(
    names_from = Reported_disease,
    values_from = c(Mean_percent, SD_percent, SEM_percent)
  )

# Merge tb3_wide with v4.fam.dunn by Family
v4.fam.combined_data.dis <- v4.fam.dunn.dis %>%
  left_join(tb3_wide.dis.fam, by = "Family")

# Save the combined data as an Excel file
#writexl::write_xlsx(v4.fam.combined_data.dis, "combined RA and KW results for disease family_20241005.xlsx")
# also save tb3 for RA results
#writexl::write_xlsx(tb3.up.fam, "RA of family disease_20241005.xlsx")


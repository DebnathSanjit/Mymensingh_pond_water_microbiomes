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

# 5. Upazila----
pseq_up <- ps %>%
  ps_filter(
    Crop_species != "Shing",
    Crop_species != "Gulsha.Carp" & Crop_species != "Gulsha.Pabda"
  ) %>% tax_fix()
#ps_mutate(
#  Upazila = recode(Upazila, # recode from dplyr package to replace the existing names
#                   "Muktagacha" = "MG",
#                   "Tarakanda" = "TK",
#                   "Jamalpur_Sadar" = "JS"))
pseq_up # 10378 taxa and 812 samples

## Phylum----
# Aggregate at phyla level
up.phyla <- tax_glom(pseq_up, taxrank = "Phylum")
up.phyla # 54 taxa and 812 samples

# calculate relative abundance
up.phy.rel <- transform_sample_counts(up.phyla, function(x) x / sum(x))# do not do this if want to see absulate abundance

# Convert to tibble
tb.up <- psmelt(up.phy.rel) %>%
  as_tibble()

### Mean and standard deviation----
tb1.up <- tb.up %>%
  group_by(Phylum, Upazila) %>%
  summarize(Mean = mean(Abundance), SD = sd(Abundance)) %>%
  ungroup()

### Standard error of the mean (SEM)----
tb2.up <- tb1.up %>%
  group_by(Phylum) %>%
  mutate(SEM = SD / sqrt(n())) %>%
  ungroup()

### Convert mean, SD, and SEM to percentages----
tb3.up <- tb2.up %>%
  mutate(
    Mean_percent = Mean * 100,
    SD_percent = SD * 100,
    SEM_percent = SEM * 100
  ) %>%
  select(Phylum, Upazila, Mean_percent, SD_percent, SEM_percent)

### Pairwise comparison (Kruskal-Wallis)of taxa----
# To perform the KW followed by dunn test, follow code from jamie, as the code form chatgpt seems giving very high p value
# Agglomerate at phylum level
V4.phy.up <- pseq_up %>% tax_glom("Phylum", NArm = FALSE) 

# Convert to relative abundance, (could use tidyPS(., TRUE) function from above)
V4.phy.ra.up <- transform_sample_counts(V4.phy.up, function(otu) otu/sum(otu))

# Keep only those groups with a mean relative abundance greater than 0.01%
V4.main.phy.up <- filter_taxa(V4.phy.ra.up, function(x) mean(x) > 0.0001, TRUE)

# Melt phyloseq object into a single tibble
V4.phy.melted.up <- as_tibble(psmelt(V4.main.phy.up), rownames = "Rows") %>%
  mutate(Abundance_adj = ifelse(Abundance == 0, NA, Abundance))   # Replace abundance zeros with NA

# Perform Kruskal-Wallis test for each phylum using the combined variable
v4.phy.kw.up <- V4.phy.melted.up %>%
  group_by(Phylum) %>%
  kruskal_test(Abundance_adj ~ Upazila)

# Perform Dunn's post-hoc test for each phylum where Kruskal-Wallis was significant
v4.phy.dunn.up <- V4.phy.melted.up %>%
  group_by(Phylum) %>% #also group by sample type to
  dunn_test(Abundance_adj ~ Upazila, p.adjust.method = "BH") #%>% 
#dplyr::mutate(FDR = p.adjust(.$p, method = "BH", 30)) #choosing number of comparisons as 30 as we have 10 phyla and three sample types from which comparisons are made
# p adjustment is needed based on how many times we are running the test.

### Combine relative abundance and test results together----
# Reshape tb3 to wide format, separating by Crop (tilaasius and Tilapia)
tb3_wide.up <- tb3.up %>%
  pivot_wider(
    names_from = Upazila,
    values_from = c(Mean_percent, SD_percent, SEM_percent)
  )

# Merge tb3_wide with v4.phy.dunn by Phylum
v4.phy.combined_data.up <- v4.phy.dunn.up %>%
  left_join(tb3_wide.up, by = "Phylum")

# Save the combined data as an Excel file
#writexl::write_xlsx(v4.phy.combined_data.up, "combined RA and KW results for upazila_20241001.xlsx")
# also save tb3 for RA results
#writexl::write_xlsx(tb3.up, "RA of phyla upazila_20241001.xlsx")

## Family----
# Aggregate at family level
up.fam <- tax_glom(pseq_up, taxrank = "Family")
up.fam # 609 taxa and 812 samples

# calculate relative abundance
up.fam.rel <- transform_sample_counts(up.fam, function(x) x / sum(x))# do not do this if want to see absulate abundance

# Convert to tibble
tb.up.fam <- psmelt(up.fam.rel) %>%
  as_tibble()

### Mean and standard deviation----
tb1.up.fam <- tb.up.fam %>%
  group_by(Family, Upazila) %>%
  summarize(Mean = mean(Abundance), SD = sd(Abundance)) %>%
  ungroup()

### Standard error of the mean (SEM)----
tb2.up.fam <- tb1.up.fam %>%
  group_by(Family) %>%
  mutate(SEM = SD / sqrt(n())) %>%
  ungroup()

### Convert mean, SD, and SEM to percentages----
tb3.up.fam <- tb2.up.fam %>%
  mutate(
    Mean_percent = Mean * 100,
    SD_percent = SD * 100,
    SEM_percent = SEM * 100
  ) %>%
  select(Family, Upazila, Mean_percent, SD_percent, SEM_percent)

### Pairwise comparison (Kruskal-Wallis)of taxa----
# To perform the KW followed by dunn test, follow code from jamie, as the code form chatgpt seems giving very high p value
## Agglomerate at Family level
#V4.fam.up <- pseq_up %>% tax_glom("Family", NArm = FALSE) 
#
## Convert to relative abundance, (could use tidyPS(., TRUE) function from above)
#V4.fam.ra.up <- transform_sample_counts(V4.fam.up, function(otu) otu/sum(otu))

# Keep only those groups with a mean relative abundance greater than 0.1%
V4.main.fam.up <- filter_taxa(up.fam.rel, function(x) mean(x) > 0.001, TRUE)

# Melt phyloseq object into a single tibble
V4.fam.melted.up <- as_tibble(psmelt(V4.main.fam.up), rownames = "Rows") %>%
  mutate(Abundance_adj = ifelse(Abundance == 0, NA, Abundance))   # Replace abundance zeros with NA

# Perform Kruskal-Wallis test for each Family using the combined variable
v4.fam.kw.up <- V4.fam.melted.up %>%
  group_by(Family) %>%
  kruskal_test(Abundance_adj ~ Upazila)

# Perform Dunn's post-hoc test for each Family where Kruskal-Wallis was significant
v4.fam.dunn.up <- V4.fam.melted.up %>%
  group_by(Family) %>% #also group by sample type to
  dunn_test(Abundance_adj ~ Upazila, p.adjust.method = "BH") #%>% 
#dplyr::mutate(FDR = p.adjust(.$p, method = "BH", 30)) #choosing number of comparisons as 30 as we have 10 family and three sample types from which comparisons are made
# p adjustment is needed based on how many times we are running the test.

### Combine relative abundance and test results together----
# Reshape tb3 to wide format, separating by Crop (tilaasius and Tilapia)
tb3_wide.up.fam <- tb3.up.fam %>%
  pivot_wider(
    names_from = Upazila,
    values_from = c(Mean_percent, SD_percent, SEM_percent)
  )

# Merge tb3_wide with v4.fam.dunn by Family
v4.fam.combined_data.up <- v4.fam.dunn.up %>%
  left_join(tb3_wide.up.fam, by = "Family")

# Save the combined data as an Excel file
#writexl::write_xlsx(v4.fam.combined_data.up, "combined RA and KW results for upazila family_20241005.xlsx")
# also save tb3 for RA results
#writexl::write_xlsx(tb3.up.fam, "RA of family for upazila_20241005.xlsx")










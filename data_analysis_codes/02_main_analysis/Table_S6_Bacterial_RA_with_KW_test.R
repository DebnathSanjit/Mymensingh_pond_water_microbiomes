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

# 1. Pangasius monoculture vs polyculture----
pseq_pang.mp <- ps %>% 
  ps_filter(
    Sampling_point > 8, 
    Crop_species != "Shing" & Crop_species != "P.no.fish",
    Crop == "Pangasius",
  ) %>% 
  ps_mutate(Culture_system2 = ifelse(grepl("\\.", Crop_species), "Polyculture", "Monoculture")) %>%  #add a column in metadata named Culture_system2 where it will look for any patern with "." in the Crop_species column and if it's true, it will mark as Polyculture, if false, Monoculture
  ps_mutate(Culture_system3 = ifelse(grepl("\\.", Crop_species), "01. Polyculture", Crop_species)) %>% #add a column in metadata named Culture.system where it will look for any patern with "." in the Crop_species column and if it's true, it will mark as Poly
  tax_fix()
pseq_pang.mp # 6518 taxa and 227 samples

## Phylum level----
# Aggregate at phyla level
pang.mp.phyla <- tax_glom(pseq_pang.mp, taxrank = "Phylum")
pang.mp.phyla # 50 taxa and 227 samples

# calculate relative abundance
pang.mp.phy.rel <- transform_sample_counts(pang.mp.phyla, function(x) x / sum(x))# do not do this if want to see absulate abundance

# Convert to tibble
tb.pang.mp <- psmelt(pang.mp.phy.rel) %>%
  as_tibble()

### Mean and standard deviation----
tb1.pang.mp <- tb.pang.mp %>%
  group_by(Phylum, Culture_system) %>%
  summarize(Mean = mean(Abundance), SD = sd(Abundance)) %>%
  ungroup()

### Standard error of the mean (SEM)----
tb2.pang.mp <- tb1.pang.mp %>%
  group_by(Phylum) %>%
  mutate(SEM = SD / sqrt(n())) %>%
  ungroup()

### Convert mean, SD, and SEM to percentages----
tb3.pang.mp <- tb2.pang.mp %>%
  mutate(
    Mean_percent = Mean * 100,
    SD_percent = SD * 100,
    SEM_percent = SEM * 100
  ) %>%
  select(Phylum, Culture_system, Mean_percent, SD_percent, SEM_percent)

### Pairwise comparison (Kruskal-Wallis)of taxa----
# To perform the KW followed by dunn test, follow code from jamie, as the code form chatgpt seems giving very high p value
# Agglomerate at phylum level
V4.phy.pang.mp <- pseq_pang.mp %>% tax_glom("Phylum", NArm = FALSE) 

# Convert to relative abundance, (could use tidyPS(., TRUE) function from above)
V4.phy.ra.pang.mp <- transform_sample_counts(V4.phy.pang.mp, function(otu) otu/sum(otu))

# Keep only those groups with a mean relative abundance greater than 0.01%
V4.main.phy.pang.mp <- filter_taxa(V4.phy.ra.pang.mp, function(x) mean(x) > 0.0001, TRUE)

# Melt phyloseq object into a single tibble
V4.phy.melted.pang.mp <- as_tibble(psmelt(V4.main.phy.pang.mp), rownames = "Rows") %>%
  mutate(Abundance_adj = ifelse(Abundance == 0, NA, Abundance))   # Replace abundance zeros with NA

# Perform Kruskal-Wallis test for each phylum using the combined variable
v4.phy.kw.pang.mp <- V4.phy.melted.pang.mp %>%
  group_by(Phylum) %>%
  kruskal_test(Abundance_adj ~ Culture_system)

# Perform Dunn's post-hoc test for each phylum where Kruskal-Wallis was significant
v4.phy.dunn.pang.mp <- V4.phy.melted.pang.mp %>%
  group_by(Phylum) %>% #also group by sample type to
  dunn_test(Abundance_adj ~ Culture_system, p.adjust.method = "BH") #%>% 
#dplyr::mutate(FDR = p.adjust(.$p, method = "BH", 30)) #choosing number of comparisons as 30 as we have 10 phyla and three sample types from which comparisons are made
# p adjustment is needed based on how many times we are running the test.

### Combine relative abundance and test results together----
# Reshape tb3 to wide format, separating by Crop (Pangasius and Tilapia)
tb3_wide.pang.mp <- tb3.pang.mp %>%
  pivot_wider(
    names_from = Culture_system,
    values_from = c(Mean_percent, SD_percent, SEM_percent)
  )

# Merge tb3_wide with v4.phy.dunn by Phylum
v4.phy.combined_data.pang.mp <- v4.phy.dunn.pang.mp %>%
  left_join(tb3_wide.pang.mp, by = "Phylum")

# Save the combined data as an Excel file
#writexl::write_xlsx(v4.phy.combined_data.pang.mp, "combined RA and KW results for pp pm_20241001.xlsx")
# also save tb3 for RA results
#writexl::write_xlsx(tb3.pang.mp, "RA of phyla pp pm_20241001.xlsx")

## Family level----
# Aggregate at family level
pang.mp.fam <- tax_glom(pseq_pang.mp, taxrank = "Family")
pang.mp.fam # 549 taxa and 227 samples

# calculate relative abundance
pang.mp.fam.rel <- transform_sample_counts(pang.mp.fam, function(x) x / sum(x))# do not do this if want to see absulate abundance

# Convert to tibble
tb.pang.mp.fam <- psmelt(pang.mp.fam.rel) %>%
  as_tibble()

### Mean and standard deviation----
tb1.pang.mp.fam <- tb.pang.mp.fam %>%
  group_by(Family, Culture_system) %>%
  summarize(Mean = mean(Abundance), SD = sd(Abundance)) %>%
  ungroup()

### Standard error of the mean (SEM)----
tb2.pang.mp.fam <- tb1.pang.mp.fam %>%
  group_by(Family) %>%
  mutate(SEM = SD / sqrt(n())) %>%
  ungroup()

### Convert mean, SD, and SEM to percentages----
tb3.pang.mp.fam <- tb2.pang.mp.fam %>%
  mutate(
    Mean_percent = Mean * 100,
    SD_percent = SD * 100,
    SEM_percent = SEM * 100
  ) %>%
  select(Family, Culture_system, Mean_percent, SD_percent, SEM_percent)

### Pairwise comparison (Kruskal-Wallis)of taxa----
# To perform the KW followed by dunn test, follow code from jamie, as the code form chatgpt seems giving very high p value
## Agglomerate at Family level
#V4.fam.pang.mp <- pseq_pang.mp %>% tax_glom("Family", NArm = FALSE) 
#
## Convert to relative abundance, (could use tidyPS(., TRUE) function from above)
#V4.fam.ra.pang.mp <- transform_sample_counts(V4.fam.pang.mp, function(otu) otu/sum(otu))
#
# Keep only those groups with a mean relative abundance greater than 0.1%
V4.main.fam.pang.mp <- filter_taxa(pang.mp.fam.rel, function(x) mean(x) > 0.001, TRUE)

# Melt phyloseq object into a single tibble
V4.fam.melted.pang.mp <- as_tibble(psmelt(V4.main.fam.pang.mp), rownames = "Rows") %>%
  mutate(Abundance_adj = ifelse(Abundance == 0, NA, Abundance))   # Replace abundance zeros with NA

# Perform Kruskal-Wallis test for each Family using the combined variable
v4.fam.kw.pang.mp <- V4.fam.melted.pang.mp %>%
  group_by(Family) %>%
  kruskal_test(Abundance_adj ~ Culture_system)


# Perform Dunn's post-hoc test for each Family where Kruskal-Wallis was significant
v4.fam.dunn.pang.mp <- V4.fam.melted.pang.mp %>%
  group_by(Family) %>% #also group by sample type to
  dunn_test(Abundance_adj ~ Culture_system, p.adjust.method = "BH") #%>% 
#dplyr::mutate(FDR = p.adjust(.$p, method = "BH", 30)) #choosing number of comparisons as 30 as we have 10 phyla and three sample types from which comparisons are made
# p adjustment is needed based on how many times we are running the test.

### Combine relative abundance and test results together----
# Reshape tb3 to wide format, separating by Crop (Pangasius and Tilapia)
tb3_wide.pang.mp.fam <- tb3.pang.mp.fam %>%
  pivot_wider(
    names_from = Culture_system,
    values_from = c(Mean_percent, SD_percent, SEM_percent)
  )

# Merge tb3_wide with v4.fam.dunn by Family
v4.fam.combined_data.pang.mp <- v4.fam.dunn.pang.mp %>%
  left_join(tb3_wide.pang.mp.fam, by = "Family")

# Save the combined data as an Excel file
#writexl::write_xlsx(v4.fam.combined_data.pang.mp, "combined RA and KW results for pp pm family_20241005.xlsx")
# also save tb3 for RA results
#writexl::write_xlsx(tb3.pang.mp.fam, "RA of phyla pp pm family_20241005.xlsx")

# 2. Tilapia monoculture vs polyculture----
pseq_tila.mp <- ps %>% 
  ps_filter(
    Sampling_point > 8, 
    Crop == "Tilapia",
    Crop_species != "T.no.fish",
    Crop_species != "Gulsha.Carp" & Crop_species != "Gulsha.Pabda"
  ) %>% 
  ps_mutate(Culture_system2 = ifelse(grepl("\\.", Crop_species), "Polyculture", "Monoculture")) %>%  #add a column in metadata named Culture_system2 where it will look for any patern with "." in the Crop_species column and if it's true, it will mark as Polyculture, if false, Monoculture
  ps_mutate(Culture_system3 = ifelse(grepl("\\.", Crop_species), "01. Polyculture", Crop_species)) %>% #add a column in metadata named Culture.system where it will look for any patern with "." in the Crop_species column and if it's true, it will mark as Poly
  tax_fix()
pseq_tila.mp # 7069 taxa and 227 samples 

## Phylum----
# Aggregate at phyla level
tila.mp.phyla <- tax_glom(pseq_tila.mp, taxrank = "Phylum")
tila.mp.phyla # 52 taxa and 227 samples

# calculate relative abundance
tila.mp.phy.rel <- transform_sample_counts(tila.mp.phyla, function(x) x / sum(x))# do not do this if want to see absulate abundance

# Convert to tibble
tb.tila.mp <- psmelt(tila.mp.phy.rel) %>%
  as_tibble()

### Mean and standard deviation----
tb1.tila.mp <- tb.tila.mp %>%
  group_by(Phylum, Culture_system) %>%
  summarize(Mean = mean(Abundance), SD = sd(Abundance)) %>%
  ungroup()

### Standard error of the mean (SEM)----
tb2.tila.mp <- tb1.tila.mp %>%
  group_by(Phylum) %>%
  mutate(SEM = SD / sqrt(n())) %>%
  ungroup()

### Convert mean, SD, and SEM to percentages----
tb3.tila.mp <- tb2.tila.mp %>%
  mutate(
    Mean_percent = Mean * 100,
    SD_percent = SD * 100,
    SEM_percent = SEM * 100
  ) %>%
  select(Phylum, Culture_system, Mean_percent, SD_percent, SEM_percent)

### Pairwise comparison (Kruskal-Wallis)of taxa----
# To perform the KW followed by dunn test, follow code from jamie, as the code form chatgpt seems giving very high p value
# Agglomerate at phylum level
V4.phy.tila.mp <- pseq_tila.mp %>% tax_glom("Phylum", NArm = FALSE) 

# Convert to relative abundance, (could use tidyPS(., TRUE) function from above)
V4.phy.ra.tila.mp <- transform_sample_counts(V4.phy.tila.mp, function(otu) otu/sum(otu))

# Keep only those groups with a mean relative abundance greater than 0.01%
V4.main.phy.tila.mp <- filter_taxa(V4.phy.ra.tila.mp, function(x) mean(x) > 0.0001, TRUE)

# Melt phyloseq object into a single tibble
V4.phy.melted.tila.mp <- as_tibble(psmelt(V4.main.phy.tila.mp), rownames = "Rows") %>%
  mutate(Abundance_adj = ifelse(Abundance == 0, NA, Abundance))   # Replace abundance zeros with NA

# Perform Kruskal-Wallis test for each phylum using the combined variable
v4.phy.kw.tila.mp <- V4.phy.melted.tila.mp %>%
  group_by(Phylum) %>%
  kruskal_test(Abundance_adj ~ Culture_system)

# Perform Dunn's post-hoc test for each phylum where Kruskal-Wallis was significant
v4.phy.dunn.tila.mp <- V4.phy.melted.tila.mp %>%
  group_by(Phylum) %>% #also group by sample type to
  dunn_test(Abundance_adj ~ Culture_system, p.adjust.method = "BH") #%>% 
#dplyr::mutate(FDR = p.adjust(.$p, method = "BH", 30)) #choosing number of comparisons as 30 as we have 10 phyla and three sample types from which comparisons are made
# p adjustment is needed based on how many times we are running the test.

### Combine relative abundance and test results together----
# Reshape tb3 to wide format, separating by Crop (tilaasius and Tilapia)
tb3_wide.tila.mp <- tb3.tila.mp %>%
  pivot_wider(
    names_from = Culture_system,
    values_from = c(Mean_percent, SD_percent, SEM_percent)
  )

# Merge tb3_wide with v4.phy.dunn by Phylum
v4.phy.combined_data.tila.mp <- v4.phy.dunn.tila.mp %>%
  left_join(tb3_wide.tila.mp, by = "Phylum")

# Save the combined data as an Excel file
#writexl::write_xlsx(v4.phy.combined_data.tila.mp, "combined RA and KW results for TP TM_20241001.xlsx")
# also save tb3 for RA results
#writexl::write_xlsx(tb3.tila.mp, "RA of phyla TP TM_20241001.xlsx")

## Family----
# Aggregate at phyla level
tila.mp.fam <- tax_glom(pseq_tila.mp, taxrank = "Family")
tila.mp.fam # 546 taxa and 227 samples

# calculate relative abundance
tila.mp.fam.rel <- transform_sample_counts(tila.mp.fam, function(x) x / sum(x))# do not do this if want to see absulate abundance

# Convert to tibble
tb.tila.mp.fam <- psmelt(tila.mp.fam.rel) %>%
  as_tibble()

### Mean and standard deviation----
tb1.tila.mp.fam <- tb.tila.mp.fam %>%
  group_by(Family, Culture_system) %>%
  summarize(Mean = mean(Abundance), SD = sd(Abundance)) %>%
  ungroup()

### Standard error of the mean (SEM)----
tb2.tila.mp.fam <- tb1.tila.mp.fam %>%
  group_by(Family) %>%
  mutate(SEM = SD / sqrt(n())) %>%
  ungroup()

### Convert mean, SD, and SEM to percentages----
tb3.tila.mp.fam <- tb2.tila.mp.fam %>%
  mutate(
    Mean_percent = Mean * 100,
    SD_percent = SD * 100,
    SEM_percent = SEM * 100
  ) %>%
  select(Family, Culture_system, Mean_percent, SD_percent, SEM_percent)

### Pairwise comparison (Kruskal-Wallis)of taxa----
# To perform the KW followed by dunn test, follow code from jamie, as the code form chatgpt seems giving very high p value
## Agglomerate at Family level
#V4.phy.tila.mp <- pseq_tila.mp %>% tax_glom("Family", NArm = FALSE) 
#
## Convert to relative abundance, (could use tidyPS(., TRUE) function from above)
#V4.phy.ra.tila.mp <- transform_sample_counts(V4.phy.tila.mp, function(otu) otu/sum(otu))

# Keep only those groups with a mean relative abundance greater than 0.1%
V4.main.fam.tila.mp <- filter_taxa(tila.mp.fam.rel, function(x) mean(x) > 0.001, TRUE)

# Melt phyloseq object into a single tibble
V4.fam.melted.tila.mp <- as_tibble(psmelt(V4.main.fam.tila.mp), rownames = "Rows") %>%
  mutate(Abundance_adj = ifelse(Abundance == 0, NA, Abundance))   # Replace abundance zeros with NA

# Perform Kruskal-Wallis test for each Family using the combined variable
v4.fam.kw.tila.mp <- V4.fam.melted.tila.mp %>%
  group_by(Family) %>%
  kruskal_test(Abundance_adj ~ Culture_system)

# Perform Dunn's post-hoc test for each Family where Kruskal-Wallis was significant
v4.fam.dunn.tila.mp <- V4.fam.melted.tila.mp %>%
  group_by(Family) %>% #also group by sample type to
  dunn_test(Abundance_adj ~ Culture_system, p.adjust.method = "BH") #%>% 
#dplyr::mutate(FDR = p.adjust(.$p, method = "BH", 30)) #choosing number of comparisons as 30 as we have 10 phyla and three sample types from which comparisons are made
# p adjustment is needed based on how many times we are running the test.

### Combine relative abundance and test results together----
# Reshape tb3 to wide format, separating by Crop (tilaasius and Tilapia)
tb3_wide.tila.mp.fam <- tb3.tila.mp.fam %>%
  pivot_wider(
    names_from = Culture_system,
    values_from = c(Mean_percent, SD_percent, SEM_percent)
  )

# Merge tb3_wide with v4.phy.dunn by Family
v4.fam.combined_data.tila.mp <- v4.fam.dunn.tila.mp %>%
  left_join(tb3_wide.tila.mp.fam, by = "Family")

# Save the combined data as an Excel file
#writexl::write_xlsx(v4.fam.combined_data.tila.mp, "combined RA and KW results for TP TM family_20241003.xlsx")
# also save tb3 for RA results
#writexl::write_xlsx(tb3.tila.mp.fam, "RA of family TP TM_20241003.xlsx")

# 3. Pangasius and Tilapia monoculture----
pseq_pt_mono <- ps %>% 
  ps_filter(Culture_system %in% c("Pangasius-monoculture", "Tilapia-monoculture"))
pseq_pt_mono # 8892 taxa and 485 samples 

## Phylum----
# Aggregate at phyla level
pt.phyla <- tax_glom(pseq_pt_mono, taxrank = "Phylum")
pt.phyla # 53 taxa and 592 samples

# calculate relative abundance
pt.phy.rel <- transform_sample_counts(pt.phyla, function(x) x / sum(x))# do not do this if want to see absulate abundance

# Convert to tibble
tb.pt <- psmelt(pt.phy.rel) %>%
  as_tibble()

### Mean and standard deviation----
tb1.pt <- tb.pt %>%
  group_by(Phylum, Crop) %>%
  summarize(Mean = mean(Abundance), SD = sd(Abundance)) %>%
  ungroup()

### Standard error of the mean (SEM)----
tb2.pt <- tb1.pt %>%
  group_by(Phylum) %>%
  mutate(SEM = SD / sqrt(n())) %>%
  ungroup()

### Convert mean, SD, and SEM to percentages----
tb3.pt <- tb2.pt %>%
  mutate(
    Mean_percent = Mean * 100,
    SD_percent = SD * 100,
    SEM_percent = SEM * 100
  ) %>%
  select(Phylum, Crop, Mean_percent, SD_percent, SEM_percent)

### Pairwise comparison (Kruskal-Wallis)of taxa----
# To perform the KW followed by dunn test, follow code from jamie, as the code form chatgpt seems giving very high p value
# Agglomerate at phylum level
V4.phy.pt <- pseq_pt_mono %>% tax_glom("Phylum", NArm = FALSE) 

# Convert to relative abundance, (could use tidyPS(., TRUE) function from above)
V4.phy.ra.pt <- transform_sample_counts(V4.phy.pt, function(otu) otu/sum(otu))

# Keep only those grodiss with a mean relative abundance greater than 0.01%
V4.main.phy.pt <- filter_taxa(V4.phy.ra.pt, function(x) mean(x) > 0.0001, TRUE)

# Melt phyloseq object into a single tibble
V4.phy.melted.pt <- as_tibble(psmelt(V4.main.phy.pt), rownames = "Rows") %>%
  mutate(Abundance_adj = ifelse(Abundance == 0, NA, Abundance))   # Replace abundance zeros with NA

# Perform Kruskal-Wallis test for each phylum using the combined variable
v4.phy.kw.pt <- V4.phy.melted.pt %>%
  group_by(Phylum) %>%
  kruskal_test(Abundance_adj ~ Crop)

# Perform Dunn's post-hoc test for each phylum where Kruskal-Wallis was significant
v4.phy.dunn.pt2 <- V4.phy.melted.pt %>%
  group_by(Phylum) %>% #also grodis by sample type to
  dunn_test(Abundance_adj ~ Crop, p.adjust.method = "BH") %>% 
  dplyr::mutate(FDR = p.adjust(.$p, method = "BH", 27)) #choosing number of comparisons as 30 as we have 10 phyla and three sample types from which comparisons are made
# p adjustment is needed based on how many times we are running the test.

### Combine relative abundance and test results together----
# Reshape tb3 to wide format, separating by Crop (tilaasius and Tilapia)
tb3_wide.pt <- tb3.pt %>%
  pivot_wider(
    names_from = Crop,
    values_from = c(Mean_percent, SD_percent, SEM_percent)
  )

# Merge tb3_wide with v4.phy.dunn by Phylum
v4.phy.combined_data.pt <- v4.phy.dunn.pt %>%
  left_join(tb3_wide.pt, by = "Phylum")

# Save the combined data as an Excel file
#writexl::write_xlsx(v4.phy.combined_data.pt, "combined RA and KW results for mono crop_20241001.xlsx")
# also save tb3 for RA results
#writexl::write_xlsx(tb3.pt, "RA of phyla mono crop_20241001.xlsx")

## Family----
# Aggregate at fam level
pt.fam <- tax_glom(pseq_pt_mono, taxrank = "Family")
pt.fam # 361 taxa and 592 samples

# calculate relative abundance
pt.fam.rel <- transform_sample_counts(pt.fam, function(x) x / sum(x))# do not do this if want to see absulate abundance

# Convert to tibble
tb.pt.fam <- psmelt(pt.fam.rel) %>%
  as_tibble()

### Mean and standard deviation----
tb1.pt.fam <- tb.pt.fam %>%
  group_by(Family, Crop) %>%
  summarize(Mean = mean(Abundance), SD = sd(Abundance)) %>%
  ungroup()

### Standard error of the mean (SEM)----
tb2.pt.fam <- tb1.pt.fam %>%
  group_by(Family) %>%
  mutate(SEM = SD / sqrt(n())) %>%
  ungroup()

### Convert mean, SD, and SEM to percentages----
tb3.pt.fam <- tb2.pt.fam %>%
  mutate(
    Mean_percent = Mean * 100,
    SD_percent = SD * 100,
    SEM_percent = SEM * 100
  ) %>%
  select(Family, Crop, Mean_percent, SD_percent, SEM_percent)

### Pairwise comparison (Kruskal-Wallis)of taxa----
# To perform the KW followed by dunn test, follow code from jamie, as the code form chatgpt seems giving very high p value
## Agglomerate at Family level
#V4.fam.pt <- pseq_pt_mono %>% tax_glom("Family", NArm = FALSE) 
#
## Convert to relative abundance, (could use tidyPS(., TRUE) function from above)
#V4.fam.ra.pt <- transform_sample_counts(V4.fam.pt, function(otu) otu/sum(otu))
#
# Keep only those grodiss with a mean relative abundance greater than 0.1%
V4.main.fam.pt <- filter_taxa(pt.fam.rel, function(x) mean(x) > 0.001, TRUE)

# Melt phyloseq object into a single tibble
V4.fam.melted.pt <- as_tibble(psmelt(V4.main.fam.pt), rownames = "Rows") %>%
  mutate(Abundance_adj = ifelse(Abundance == 0, NA, Abundance))   # Replace abundance zeros with NA

# Perform Kruskal-Wallis test for each Family using the combined variable
v4.fam.kw.pt <- V4.fam.melted.pt %>%
  group_by(Family) %>%
  kruskal_test(Abundance_adj ~ Culture_system)

# Perform Dunn's post-hoc test for each Family where Kruskal-Wallis was significant
v4.fam.dunn.pt <- V4.fam.melted.pt %>%
  group_by(Family) %>% #also grodis by sample type to
  dunn_test(Abundance_adj ~ Crop, p.adjust.method = "BH") #%>% 
#dplyr::mutate(FDR = p.adjust(.$p, method = "BH", 30)) #choosing number of comparisons as 30 as we have 10 fam and three sample types from which comparisons are made
# p adjustment is needed based on how many times we are running the test.

### Combine relative abundance and test results together----
# Reshape tb3 to wide format, separating by Crop (tilaasius and Tilapia)
tb3.pt.fam_wide <- tb3.pt.fam %>%
  pivot_wider(
    names_from = Crop,
    values_from = c(Mean_percent, SD_percent, SEM_percent)
  )

# Merge tb3_wide with v4.fam.dunn by Family
v4.fam.combined_data.pt <- v4.fam.dunn.pt %>%
  left_join(tb3.pt.fam_wide, by = "Family")

# Save the combined data as an Excel file
#writexl::write_xlsx(v4.fam.combined_data.pt, "combined RA and KW results for mono crop family_20241005.xlsx")
# also save tb3 for RA results
#writexl::write_xlsx(tb3.pt.fam, "RA of fam mono crop family_20241005.xlsx")

# Combine all these results to make Table S6
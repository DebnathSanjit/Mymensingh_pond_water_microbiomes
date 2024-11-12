# 21/09/2024
# Sanjit Debnath
# This script is to do beta diversity to see how season affect microbial communities

# Load Libraries
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(microViz); packageVersion("microViz")
library(vegan); packageVersion("vegan") # needed for PERMANOVA test
library(pairwiseAdonis); packageVersion("pairwiseAdonis") # pairwise comparison
library(tidyverse); packageVersion("tidyverse")
library(ggpubr); packageVersion("ggpubr")
library(cowplot); packageVersion("cowplot")
library(stringr); packageVersion("stringr") # to wrap text
# Libraries for tests
library(rstatix); packageVersion("rstatix")
#library(dunn.test); packageVersion("dunn.test")
library(openxlsx) # to write workbook in excel

## Setup working dictionary first
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/statistical_analysis")

# set seed
set.seed(1234)
# Set theme for ggplot
theme_set(theme_bw())

# Customised colours
season.colors <- c("Monsoon" = "#EC823C","Winter" ="#1B9E77","Pre-monsoon" = "#100AFF")
season.colors <- c("Monsoon1" = "#EC823C","Winter1" ="#1B9E77","Pre_monsoon1" = "#100AFF", "Monsoon2" = "#EC823C","Winter2" ="#1B9E77","Pre_monsoon2" = "#100AFF")

# Prokaryotic----
# Load the phyloseq object
ps <- readRDS("phyloseq_metadata_6_v3_20240419.rds")
ps # 10523 taxa and 891 samples

# Pangasius ponds only----
pseq_pang <- ps %>% 
  ps_filter(
    Crop == "Pangasius" & Crop_species != "Shing") %>% 
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
      TRUE ~ as.character(Season)))  # handle any other cases, if they exist
pseq_pang # 8091 taxa and 421 samples

## B1.Beta diversity (weighted unifrac)----
#wunifrac:weighted UniFrac took about #01:59 - 02:12 = 13 minutes to calculate the distance
#p_wuni <- pseq_pang %>%
#  tax_transform(rank = "unique", trans = "compositional") %>%
#  dist_calc(dist = "wunifrac") 
#saveRDS(p_wuni, "weighted_unifrac_beta_diversity_pseq_pang.rds")
p_wuni <- readRDS("weighted_unifrac_beta_diversity_pseq_pang.rds")

# beta diversity: Season
p.s.we <- p_wuni %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Season3", fill = "Season3",
    shape = 16, 
    alpha = 1,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(colour = Season3)
  ) +
  #scale_colour_manual(values = season.colors) +
  #scale_fill_manual(values = season.colors) +
  theme(
    plot.title = element_text(hjust = 0, size = 18, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 18, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank(),
    #panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  labs(colour = "Season",
       fill = "Season",
       title = "Prokaryotes: Pangasius pond") # Add custom plot label
  #xlim(-1.5, 2.2)
p.s.we

#### Significance Test----
#Plot PERMANOVA with phyloseq
metadata_pang <- sample_data(pseq_pang) %>%
  data.frame() %>%
  tibble()
wuni_dist.p = phyloseq::distance(pseq_pang, method = "wunifrac")
PERM_s_wuni.p <- adonis2(wuni_dist.p ~ Season3, data = metadata_pang)

# to add test result on the plot
p.s.we2 <- p.s.we + annotate(geom = "label",
                             label = paste("PERMANOVA: R² = ", round(PERM_s_wuni.p["Season3","R2"], 3), 
                                           ", p = ", PERM_s_wuni.p["Season3", "Pr(>F)"], sep = ""),
                             x=Inf, y=Inf, hjust = 1, vjust = 1)
p.s.we2

#### Perform pairwise PERMANOVA----
#pairwise_perm <- pairwise.adonis(wuni_dist ~ Season2, data = metadata_pang)
pairwise.perm.season <- pairwise.adonis(wuni_dist.p, metadata_pang$Season3)
pairwise.perm.season

# Write the results to a CSV file
#write.csv(pairwise.perm.season, file = "pairwise permanova for season for pangasius_20240516.csv", row.names = TRUE)

# Tilapia pond water----
pseq_tila <- ps %>% 
  ps_filter(Crop == "Tilapia",
            Crop_species != "Gulsha.Carp" & Crop_species != "Gulsha.Pabda") %>% 
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
pseq_tila #  8380 taxa and 391 samples

## B2. Beta diversity (weighted unifrac)----
# wunifrac:weighted UniFrac took about #02:58 - 03:09 = 11 minutes to calculate the distance
#t_wuni <- pseq_tila %>%
#  tax_transform(rank = "unique", trans = "compositional") %>%
#  dist_calc(dist = "wunifrac") 
#saveRDS(t_wuni, "weighted_unifrac_beta_diversity_pseq_tila.rds")
t_wuni <- readRDS("weighted_unifrac_beta_diversity_pseq_tila.rds")

# beta diversity: Season
t.s.we <- t_wuni %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Season3", fill = "Season3",
    shape = 16, 
    alpha = 1,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(colour = Season3)
  ) +
  #scale_colour_manual(values = season.colors) +
  #scale_fill_manual(values = season.colors) +
  theme(
    plot.title = element_text(hjust = 0, size = 18, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 18, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank(),
    #panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  labs(colour = "Season",
       fill = "Season",
       title = "Prokaryotes: Tilapia pond") # Add custom plot label
  #xlim(-1.5, 2)
t.s.we

### Significance Test----
#Plot PERMANOVA with phyloseq
metadata_tila <- sample_data(pseq_tila) %>%
  data.frame() %>%
  tibble()
wuni_dist.t = phyloseq::distance(pseq_tila, method = "wunifrac")
PERM_s_wuni.t <- adonis2(wuni_dist.t ~ Season3, data = metadata_tila)

# to add test result on the plot
t.s.we2 <- t.s.we + annotate(geom = "label",
                             label = paste("PERMANOVA: R² = ", round(PERM_s_wuni.t["Season3","R2"], 3), 
                                           ", p = ", PERM_s_wuni.t["Season3", "Pr(>F)"], sep = ""),
                             x=Inf, y=Inf, hjust = 1, vjust = 1)
t.s.we2

#### Perform pairwise PERMANOVA----
pairwise.perm.season.t <- pairwise.adonis(wuni_dist.t, metadata_tila$Season3)
pairwise.perm.season.t
# Write the results to a CSV file
#write.csv(pairwise.perm.season.t, file = "pairwise permanova results for season for tilapia_20240516.csv", row.names = TRUE)

# Microeukaryotic----
ps.18s <- readRDS("phyloseq_18S_filtered_with_tree_pr2_90-150bp_20240416.rds")
ps.18s # 5390 taxa and 872 samples

## Pangasius ponds----
pseq_pang.18s <- ps.18s %>% 
  ps_filter(
    Crop == "Pangasius" & Crop_species != "Shing") %>% 
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
pseq_pang.18s # 3947 taxa and 421 samples

## B3. Beta diversity (weighted unifrac)----
# wunifrac:weighted UniFrac took about #01:59 - 02:12 = 13 minutes to calculate the distance
#p_wuni.18s <- pseq_pang.18s %>%
#  tax_transform(rank = "unique", trans = "compositional") %>%
#  dist_calc(dist = "wunifrac") 
#saveRDS(p_wuni.18s, "weighted_unifrac_beta_diversity_pseq_pang_18s.rds")
p_wuni.18s <- readRDS("weighted_unifrac_beta_diversity_pseq_pang_18s.rds")

# beta diversity: Season
p.s.we.18s <- p_wuni.18s %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Season3", fill = "Season3",
    shape = 21, 
    alpha = 1,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(colour = Season3)
  ) +
  #scale_colour_manual(values = season.colors) +
  #scale_fill_manual(values = season.colors) +
  theme(
    plot.title = element_text(hjust = 0, size = 18, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 18, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank(),
    #panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  labs(colour = "Season",
       fill = "Season",
       title = "Microeukaryotes: Pangasius pond")  # Add custom plot label
p.s.we.18s

### Significance Test----
#Plot PERMANOVA with phyloseq
metadata_pang.18s <- sample_data(pseq_pang.18s) %>%
  data.frame() %>%
  tibble()
wuni_dist.p.18s = phyloseq::distance(pseq_pang.18s, method = "wunifrac")
PERM_s_wuni.p.18s <- adonis2(wuni_dist.p.18s ~ Season3, data = metadata_pang.18s)

# to add test result on the plot
p.s.we2.18s <- p.s.we.18s + annotate(geom = "label",
                                     label = paste("PERMANOVA: R² = ", round(PERM_s_wuni.p.18s["Season3","R2"], 3), 
                                                   ", p = ", PERM_s_wuni.p.18s["Season3", "Pr(>F)"], sep = ""),
                                     x=Inf, y=Inf, hjust = 1, vjust = 1)
p.s.we2.18s

#### Perform pairwise PERMANOVA----
pairwise.perm.season.18s <- pairwise.adonis(wuni_dist.p.18s, metadata_pang.18s$Season3)
pairwise.perm.season.18s
# Write the results to a CSV file
#write.csv(pairwise.perm.season.18s, file = "pairwise permanova for season for pangasius.18s_20240517.csv", row.names = TRUE)

## Tilapia pond water----
pseq_tila.18s <- ps.18s %>% 
  ps_filter(Crop == "Tilapia",
            Crop_species != "Gulsha.Carp" & Crop_species != "Gulsha.Pabda") %>% 
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
pseq_tila.18s #  8380 taxa and 391 samples

## B4. Beta diversity (weighted unifrac)----
# wunifrac:weighted UniFrac took about #02:58 - 03:09 = 11 minutes to calculate the distance
#t_wuni.18s <- pseq_tila.18s %>%
#  tax_transform(rank = "unique", trans = "compositional") %>%
#  dist_calc(dist = "wunifrac") 
#saveRDS(t_wuni.18s, "weighted_unifrac_beta_diversity_pseq_tila_18s.rds")
t_wuni.18s <- readRDS("weighted_unifrac_beta_diversity_pseq_tila_18s.rds")

# beta diversity: Season
t.s.we.18s <- t_wuni.18s %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Season3", fill = "Season3",
    shape = 21, 
    alpha = 1,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(colour = Season3)
  ) +
  #scale_colour_manual(values = season.colors) +
  #scale_fill_manual(values = season.colors) +
  theme(
    plot.title = element_text(hjust = 0, size = 18, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 18, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank(),
    #panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  labs(colour = "Season",
       fill = "Season",
       title = "Microeukaryotes: Tilapia pond") # Add custom plot label
t.s.we.18s

### Significance Test----
#Plot PERMANOVA with phyloseq
metadata_tila.18s <- sample_data(pseq_tila.18s) %>%
  data.frame() %>%
  tibble()
wuni_dist.t.18s = phyloseq::distance(pseq_tila.18s, method = "wunifrac")
PERM_s_wuni.t.18s <- adonis2(wuni_dist.t.18s ~ Season3, data = metadata_tila.18s)

# to add test result on the plot
t.s.we2.18s <- t.s.we.18s + annotate(geom = "label",
                                     label = paste("PERMANOVA: R² = ", round(PERM_s_wuni.t.18s["Season3","R2"], 3), 
                                                   ", p = ", PERM_s_wuni.t.18s["Season3", "Pr(>F)"], sep = ""),
                                     x=Inf, y=Inf, hjust = 1, vjust = 1)
t.s.we2.18s

#### Perform pairwise PERMANOVA----
pairwise.perm.season.t.18s <- pairwise.adonis(wuni_dist.t.18s, metadata_tila.18s$Season3)
pairwise.perm.season.t.18s
# Write the results to a CSV file
#write.csv(pairwise.perm.season.t.18s, file = "pairwise permanova results for season for tilapia.18s_20240517.csv", row.names = TRUE)

# combine beta diversity----
comb.month.beta <- cowplot::plot_grid(
  p.s.we2 + theme(legend.position = "none"), 
  t.s.we2 + theme(legend.position = "none"), 
  p.s.we2.18s + theme(legend.position = "none"), 
  t.s.we2.18s + theme(legend.position = "none"), 
  ncol = 2, nrow = 2, labels = "AUTO")
comb.month.beta

## get legend
legend.beta <- get_legend(p.s.we)

# Final plot----
comb.month.beta2 <- cowplot::plot_grid(comb.month.beta,
                                       legend.beta, 
                                       ncol = 1,
                                       rel_heights = c(2, 0.2))
comb.month.beta2
# Save as 1600*1200

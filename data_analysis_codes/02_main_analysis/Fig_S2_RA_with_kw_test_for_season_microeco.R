# Date: 08/10/2024 
# author: Sanjit Debnath
# This script to plot kruskal-wallis with post hoc Dunn test using microeco package to see how microbial phyla and family vary on different seasons.
# The advantages with this package is, we can adjust the p value. so I tried with tax_fix(). I will use with tax_fix, otherwise few taxa get missing

# Load Libraries
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(microViz); packageVersion("microViz")
library(tidyverse); packageVersion("tidyverse")
library(cowplot); packageVersion("cowplot")
library(patchwork); packageVersion("patchwork")
library(colorblindr); packageVersion("colorblindr")
library(paletteer); packageVersion("paletteer")
library(file2meco); packageVersion("file2meco") # to convent into microtable
library(microeco); packageVersion("microeco")

## Setup working dictionary first
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/statistical_analysis")

# set seed
set.seed(1234)

# Set theme for ggplot
theme_set(theme_bw())

# Prokaryotes----
# Load the phyloseq object
ps <- readRDS("phyloseq_metadata_6_v3_20240419.rds")

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

# Manually define the correct order of months
month_order <- c("Jun-16","Oct-16","Nov-16","Dec-16","Jan-17","Feb-17","Mar-17",
                 "Apr-17","May-17","Jun-17","Jul-17","Aug-17","Sep-17","Oct-17", 
                 "Nov-17","Dec-17","Jan-18","Feb-18","Mar-18","Apr-18","May-18")
season_order2 <- c("Monsoon", "Winter", "Pre-monsoon")
season_order3 <- c("Monsoon1", "Winter1", "Pre_monsoon1","Monsoon2", "Winter2", "Pre_monsoon2")

# Access the sample data from the phyloseq object and modify the Sampling_months column
sample_data(ps)$Sampling_months <- factor(sample_data(ps)$Sampling_months, levels = month_order)
sample_data(ps)$Season2 <- factor(sample_data(ps)$Season2, levels = season_order2)
sample_data(ps)$Season3 <- factor(sample_data(ps)$Season3, levels = season_order3)

# Pangasius ponds only----
pseq_pang <- ps %>% 
  ps_filter(
    Crop == "Pangasius" & Crop_species != "Shing") %>% 
  tax_fix()
pseq_pang # 8091 taxa and 421 samples

# Convert to microeco package
dataset.pang <- phyloseq2meco(pseq_pang)

### Run KW_dunn----
#### Phylum----
kwd.pang.phy <- trans_diff$new(dataset = dataset.pang, 
                               method = "KW_dunn", 
                               group = "Season3", 
                               taxa_level = "Phylum",
                               alpha = 0.05, # significance threshold or p value
                               p_adjust_method = "bh") 

##### Relative abundance plot
pang.phy.bar <- kwd.pang.phy$plot_diff_abund(use_number = 1:10, 
                                             coord_flip = FALSE,
                                             keep_prefix = FALSE,
                                             add_sig = TRUE, 
                                             text_y_size = 12) + 
                theme(
                plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
                axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5, size = 16, margin = margin(t = 5)),
                axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
                axis.title.x = element_blank(), #text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
                axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
                legend.position = "right",
                legend.box = "vertical",
                legend.title = element_text(size = 18, face = "bold"),
                legend.background = element_rect(fill = "white", color = "black"),
                legend.key = element_blank(),
                legend.spacing.y = unit(1, "cm"),
                legend.text = element_text(size = 16),
                panel.grid = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(color = "black")) +
  labs(fill = "Season")
pang.phy.bar

#### Family----
kwd.pang.fam <- trans_diff$new(dataset = dataset.pang, 
                               method = "KW_dunn", 
                               group = "Season3", 
                               taxa_level = "Family",
                               alpha = 0.05, # significance threshold or p value
                               p_adjust_method = "bh") 

##### Relative abundance plot
pang.fam.bar <- kwd.pang.fam$plot_diff_abund(use_number = 1:10, coord_flip = FALSE,
                                             keep_prefix = FALSE,
                                             add_sig = TRUE, 
                                             text_y_size = 12) + 
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5, size = 16, face = "italic", margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(), #text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 16),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  labs(fill = "Season")
pang.fam.bar

# Tilapia ponds only----
pseq_tila <- ps %>% 
  ps_filter(Crop == "Tilapia",
            Crop_species != "Gulsha.Carp" & Crop_species != "Gulsha.Pabda") %>% 
  tax_fix()
pseq_tila #  8380 taxa and 391 samples

## Convert to microeco
dataset.tila <- phyloseq2meco(pseq_tila)

### Run KW_dunn----
#### Phylum----
kwd.tila.phy <- trans_diff$new(dataset = dataset.tila, 
                               method = "KW_dunn", 
                               group = "Season3", 
                               taxa_level = "Phylum",
                               alpha = 0.05, # significance threshold or p value
                               p_adjust_method = "bh") 

##### Relative abundance plot
tila.phy.bar <- kwd.tila.phy$plot_diff_abund(use_number = 1:10, coord_flip = FALSE,
                                             keep_prefix = FALSE,
                                             add_sig = TRUE, 
                                             text_y_size = 12) + 
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5, size = 16, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(), #text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 16),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  labs(fill = "Season")
tila.phy.bar

#### Family----
kwd.tila.fam <- trans_diff$new(dataset = dataset.tila, 
                               method = "KW_dunn", 
                               group = "Season3", 
                               taxa_level = "Family",
                               alpha = 0.05, # significance threshold or p value
                               p_adjust_method = "bh") 

##### Relative abundance plot----
tila.fam.bar <- kwd.tila.fam$plot_diff_abund(use_number = 1:10, coord_flip = FALSE,
                                             keep_prefix = FALSE,
                                             add_sig = TRUE, 
                                             text_y_size = 12) + 
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5, size = 16, face = "italic", margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(), #text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 16),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  labs(fill = "Season")
tila.fam.bar

# Combine v4 plots----
season.bar.16s <- cowplot::plot_grid(pang.phy.bar + theme(legend.position = "none"),
                   tila.phy.bar + theme(legend.position = "none"),
                   pang.fam.bar + theme(legend.position = "none"),
                   tila.fam.bar + theme(legend.position = "none"),
                   labels = "AUTO")
season.bar.16s

# Microeukaryotes----
ps.18s <- readRDS("phyloseq_18S_filtered_with_tree_pr2_90-150bp_20240416.rds")
ps.18s <- ps.18s %>% 
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
ps.18s # 5390 taxa and 872 samples

# Manually define the correct order of months
month_order <- c("Jun-16","Oct-16","Nov-16","Dec-16","Jan-17","Feb-17","Mar-17",
                 "Apr-17","May-17","Jun-17","Jul-17","Aug-17","Sep-17","Oct-17", 
                 "Nov-17","Dec-17","Jan-18","Feb-18","Mar-18","Apr-18","May-18")

# Access the sample data from the phyloseq object and modify the Sampling_months column
sample_data(ps.18s)$Sampling_months <- factor(sample_data(ps.18s)$Sampling_months, levels = month_order)
sample_data(ps.18s)$Season2 <- factor(sample_data(ps.18s)$Season2, levels = season_order2)
sample_data(ps.18s)$Season3 <- factor(sample_data(ps.18s)$Season3, levels = season_order3)

# Pangasius ponds only----
pseq_pang.18s <- ps.18s %>% 
  ps_filter(
    Crop == "Pangasius" & Crop_species != "Shing") %>% 
  tax_fix()
pseq_pang.18s # 3947 taxa and 421 samples

# Convert to microtable
dataset.pang.18s <- phyloseq2meco(pseq_pang.18s)

### Run KW_dunn----
#### Phylum----
kwd.pang.div <- trans_diff$new(dataset = dataset.pang.18s, 
                               method = "KW_dunn", 
                               group = "Season3", 
                               taxa_level = "Division",
                               alpha = 0.05, # significance threshold or p value
                               p_adjust_method = "bh") 

##### Relative abundance plot
pang.div.bar <- kwd.pang.div$plot_diff_abund(use_number = 1:10, coord_flip = FALSE,
                                             keep_prefix = FALSE,
                                             add_sig = TRUE, 
                                             text_y_size = 12) + 
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5, size = 16, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(), #text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 16),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  labs(fill = "Season")
pang.div.bar

#### Family----
kwd.pang.fam.18s <- trans_diff$new(dataset = dataset.pang.18s, 
                                   method = "KW_dunn", 
                                   group = "Season3", 
                                   taxa_level = "Family",
                                   alpha = 0.05, # significance threshold or p value
                                   p_adjust_method = "bh") 

##### Relative abundance plot
pang.fam.18s.bar <- kwd.pang.fam.18s$plot_diff_abund(use_number = 1:10, coord_flip = FALSE,
                                                     keep_prefix = FALSE,
                                                     add_sig = TRUE, 
                                                     text_y_size = 12) + 
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5, size = 16, face = "italic", margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(), #text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 16),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  labs(fill = "Season")
pang.fam.18s.bar

# Tilapia ponds only----
pseq_tila.18s <- ps.18s %>% 
  ps_filter(Crop == "Tilapia",
            Crop_species != "Gulsha.Carp" & Crop_species != "Gulsha.Pabda") %>% 
  tax_fix() 
pseq_tila.18s #  4463 taxa and 376 samples

# Convert to microtable
dataset.tila.18s <- phyloseq2meco(pseq_tila.18s)

### Run KW_dunn----
#### Phylum----
kwd.tila.div <- trans_diff$new(dataset = dataset.tila.18s, 
                               method = "KW_dunn", 
                               group = "Season3", 
                               taxa_level = "Division",
                               alpha = 0.05, # significance threshold or p value
                               p_adjust_method = "bh") 

##### Relative abundance plot
tila.div.bar <- kwd.tila.div$plot_diff_abund(use_number = 1:10, coord_flip = FALSE,
                                             keep_prefix = FALSE,
                                             add_sig = TRUE, 
                                             text_y_size = 12) + 
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5, size = 16, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(), #text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 16),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  labs(fill = "Season")
tila.div.bar

#### Family----
kwd.tila.18s.fam <- trans_diff$new(dataset = dataset.tila.18s, 
                                   method = "KW_dunn", 
                                   group = "Season3", 
                                   taxa_level = "Family",
                                   alpha = 0.05, # significance threshold or p value
                                   p_adjust_method = "bh") 

##### Relative abundance plot----
tila.fam.18s.bar <- kwd.tila.18s.fam$plot_diff_abund(use_number = 1:10, coord_flip = FALSE,
                                                     keep_prefix = FALSE,
                                                     add_sig = TRUE, 
                                                     text_y_size = 12) + 
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5, size = 16, face = "italic", margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(), #text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 16),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  labs(fill = "Season") +
  guides(fill = guide_legend(nrow = 1))
tila.fam.18s.bar

# Combine v9 plot----
season.bar.18s <- cowplot::plot_grid(pang.div.bar + theme(legend.position = "none"),
                                    tila.div.bar + theme(legend.position = "none"),
                                    pang.fam.18s.bar + theme(legend.position = "none"),
                                    tila.fam.18s.bar + theme(legend.position = "none"),
                                    labels = c("E", "F", "G", "H"))
season.bar.18s

# Combine all
season.bar <- cowplot::plot_grid(season.bar.16s,
                   season.bar.18s,
                   ncol = 1)
legend.season <- get_legend(tila.fam.18s.bar)

season.bar2 <- cowplot::plot_grid(season.bar,
                                  legend.season,
                                  ncol = 1,
                                  rel_heights = c(2, 0.1))
season.bar2
# Save both as 2200*1800 



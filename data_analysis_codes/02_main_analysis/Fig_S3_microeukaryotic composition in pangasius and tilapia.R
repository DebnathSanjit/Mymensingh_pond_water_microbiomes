# 10/09/2024
# Sanjit Debnath
# This script is to plot microeukaryotic composition in pangasius and tilapia pond

# Load Libraries
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(microViz); packageVersion("microViz")
library(tidyverse); packageVersion("tidyverse")
library(cowplot); packageVersion("cowplot")
library(colorblindr); packageVersion("colorblindr")
library(paletteer); packageVersion("paletteer")

## Setup working dictionary first
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/statistical_analysis")

# Set theme for ggplot
theme_set(theme_bw())
set.seed(1234)

# Microeukaryotes----
ps.18s <- readRDS("phyloseq_18S_filtered_with_tree_pr2_90-150bp_20240416.rds")
ps.18s # 5390 taxa and 872 samples

# Manually define the correct order of months
month_order <- c("Jun-16","Oct-16","Nov-16","Dec-16","Jan-17","Feb-17","Mar-17",
                 "Apr-17","May-17","Jun-17","Jul-17","Aug-17","Sep-17","Oct-17", 
                 "Nov-17","Dec-17","Jan-18","Feb-18","Mar-18","Apr-18","May-18")

# Access the sample data from the phyloseq object and modify the Sampling_months column
sample_data(ps.18s)$Sampling_months <- factor(sample_data(ps.18s)$Sampling_months, levels = month_order)

## Subset pangasius (18s)----
pseq_pang.18s <- ps.18s %>% 
  ps_filter(
    Crop == "Pangasius",
    Crop_species != "Shing"
  ) %>% tax_fix()
pseq_pang.18s # 3947 taxa and 421 samples

## Division----
# my palette
myPal.div.18s <- tax_palette(
  data = pseq_pang.18s, rank = "Division", n = 20, pal = "kelly",
  add = c(Other = "lightgrey")
) 
# Override existing values
myPal.div.18s["Eukaryota Domain"] <- "#00D9FF"
myPal.div.18s["Chlorophyta"] <- "#ABE496"
tax_palette_plot(myPal.div.18s)

### Pangasius pond----
topTaxa.pang.div <- pseq_pang.18s %>%
  tax_top(n = 10, rank = "Division") %>% # >= 1%
  sort() # this makes them alphabetical

## plot with alphabetical sorting
p2.18s.div <- pseq_pang.18s %>%
  tax_sort(by = sum, at = "Division") %>% # this orders all genera by abundance
  ps_select(Crop_species, Sampling_months) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Sampling_months") %>%
  comp_barplot(tax_order = topTaxa.pang.div, # this brings the named taxa to the front
               tax_level = "Division", n_taxa = 10, # RA >= 1%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.div.18s) +
  #coord_flip() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, margin = margin(b = 15)), # Increase margin for plot title
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_blank(),  # to remove x axis title
        axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
        legend.position = "right",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18)) +
  labs(title = "Pangasius ponds",
       y = "Relative abundance") +
  guides(fill = guide_legend(ncol = 1))  # Adjust the number of rows in the legend
p2.18s.div

## Subset tilapia----
pseq_tila.18s <- ps.18s %>% 
  ps_filter(
    Crop == "Tilapia",
    Crop_species != "Gulsha.Carp", 
    Crop_species != "Gulsha.Pabda"
  ) %>% tax_fix()
pseq_tila.18s # 4463 taxa and 376 samples

### Tilapia pond_phylum----
# set up for alphabetical sorting
topTaxa.tila.div <- pseq_tila.18s %>%
  tax_top(n = 9, rank = "Division") %>% # >=1%
  sort() # this makes them alphabetical

## plot with alphabetical sorting
t2.18s.div <- pseq_tila.18s %>%
  tax_sort(by = sum, at = "Division") %>% # this orders all genera by abundance
  ps_select(Crop_species, Sampling_months) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Sampling_months") %>%
  comp_barplot(tax_order = topTaxa.tila.div, # this brings the named taxa to the front
               tax_level = "Division", n_taxa = 9, # RA >= 1%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.div.18s) +
  #coord_flip() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, margin = margin(b = 15)), # Increase margin for plot title
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_blank(),  # to remove x axis title
        axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
        legend.position = "right",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18)) +
  labs(title = "Tilapia ponds",
       y = "Relative abundance") +
  guides(fill = guide_legend(ncol = 1))  # Adjust the number of rows in the legend
t2.18s.div

## Family----
# my palette
myPal.fam.18s <- tax_palette(
  data = pseq_pang.18s, rank = "Family", n = 35, pal = "brewerPlus",
  add = c(Other = "lightgrey"))
myPal.fam.18s["Alveolata Division"] <- "#CC00A7"
myPal.fam.18s["Peridiniales Order"] <- "#EA70FF"
myPal.fam.18s["Saccharomycetales"] <- "#c29545" 
myPal.fam.18s["Tovelliaceae"] <- "#5e79b2"
tax_palette_plot(myPal.fam.18s)

### Pangasius pond----
# set up for alphabetical sorting
topTaxa.pang.fam.18s <- pseq_pang.18s %>%
  tax_top(n = 23, rank = "Family") %>% # >= 1%
  sort() # this makes them alphabetical

## plot with alphabetical sorting
p2.18s.fam <- pseq_pang.18s %>%
  tax_sort(by = sum, at = "Family") %>% # this orders all genera by abundance
  ps_select(Crop_species, Sampling_months) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Sampling_months") %>%
  comp_barplot(tax_order = topTaxa.pang.fam.18s, # this brings the named taxa to the front
               tax_level = "Family", n_taxa = 23, # RA >= 1%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.fam.18s) +
  #coord_flip() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, margin = margin(b = 15)), # Increase margin for plot title
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_blank(),  # to remove x axis title
        axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
        legend.position = "right",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18, face = "italic")) +
  labs(title = "Pangasius ponds",
       y = "Relative abundance") +
  guides(fill = guide_legend(ncol = 1))  # Adjust the number of rows in the legend
p2.18s.fam

### Tilapia pond----
# set up for alphabetical sorting
topTaxa.tila.fam.18s <- pseq_tila.18s %>%
  tax_top(n = 25, rank = "Family") %>% # >= 1%
  sort() # this makes them alphabetical

## plot with alphabetical sorting
t2.18s.fam <- pseq_tila.18s %>%
  tax_sort(by = sum, at = "Family") %>% # this orders all genera by abundance
  ps_select(Crop_species, Sampling_months) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Sampling_months") %>%
  comp_barplot(tax_order = topTaxa.tila.fam.18s, # this brings the named taxa to the front
               tax_level = "Family", n_taxa = 25, # RA >= 1%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.fam.18s) +
  #coord_flip() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, margin = margin(b = 15)), # Increase margin for plot title
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_blank(),  # to remove x axis title
        axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
        legend.position = "right",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18, face = "italic")) +
  labs(title = "Tilapia ponds",
       y = "Relative abundance") +
  guides(fill = guide_legend(ncol = 1))  # Adjust the number of rows in the legend
t2.18s.fam

#### Combine microeukaryotic----
comb.18s <- cowplot::plot_grid(p2.18s.div + theme(legend.position = "none"),#+ guides(fill = guide_legend(ncol = 1)),
                               t2.18s.div + theme(legend.position = "none"),#+ guides(fill = guide_legend(ncol = 1)),
                               p2.18s.fam + theme(legend.position = "none"),#+ guides(fill = guide_legend(ncol = 1)),
                               t2.18s.fam + theme(legend.position = "none"),#+ guides(fill = guide_legend(ncol = 1)),
                               labels = "AUTO",
                               ncol = 2)
comb.18s

# Add legend for now, need to update in inkscape
legend.fam.18s <- get_legend(t2.18s.fam)
comb.18s2 <- cowplot::plot_grid(comb.18s,
                                legend.fam.18s,
                                ncol = 2,
                                rel_widths = c(2,0.43))
comb.18s2 # Save as 2000*1250 as final figure

# Combinedly plot all legends
legend.p2.d <- get_legend(p2.18s.div)
legend.t2.d <- get_legend(t2.18s.div)
legend.p2.fam <- get_legend(p2.18s.fam)
legend.t2.fam <- get_legend(t2.18s.fam)

comb.legend.18s <- cowplot::plot_grid(legend.p2.d,# + theme(legend.position = "right")+ guides(fill = guide_legend(ncol = 1)),
                                  legend.t2.d,# + theme(legend.position = "right")+ guides(fill = guide_legend(ncol = 1)),
                                  legend.p2.fam,#p2.fam + theme(legend.position = "right")+ guides(fill = guide_legend(ncol = 1)),
                                  legend.t2.fam,#t2.fam + theme(legend.position = "right")+ guides(fill = guide_legend(ncol = 1)),
                                  #labels = "AUTO",
                                  ncol = 2)
comb.legend.18s # Save as 2000*1250 and then edit the legend to add all the missing together and add this legend to the previous plot




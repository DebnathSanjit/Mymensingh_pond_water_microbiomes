# Date: 17/09/2024 
# Sanjit Debnath
# This code is to plot alpha and beta diversity for pond with and without fish for both prokaryotes and microeukaryotes together

# Load Libraries
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(microViz); packageVersion("microViz")
library(vegan); packageVersion("vegan") # needed for PERMANOVA test
library(tidyverse); packageVersion("tidyverse")
library(ggpubr); packageVersion("ggpubr")
library(cowplot); packageVersion("cowplot")
library(stringr); packageVersion("stringr") # to wrap text
# Libraries for tests
library(ggpubr); packageVersion("ggpubr")
library(rstatix); packageVersion("rstatix")
#library(dunn.test); packageVersion("dunn.test")

## Setup working dictionary first
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/statistical_analysis")

# Set theme for ggplot
theme_set(theme_bw())
# set seed
set.seed(1234)

# Prokaryotes----
# Load the phyloseq object
ps <- readRDS("phyloseq_metadata_6_v3_20240419.rds")
ps #10523 taxa and 891 samples

# set colors
fish.colors <- c("With_fish" = "#654CFFFF", "Without_fish" = "#920000FF")

## 1 month without fish----
# subset samples there there was no fish for one month compare with before and after no fish
pseq_fish1 <- ps %>%   
  ps_filter(
    Pond_name == "PA" & Sampling_point %in% c(11,12)|
      Pond_name == "PB" & Sampling_point %in% c(7, 8)|
      Pond_name == "PC" & Sampling_point %in% c(13,14,15)|
      Pond_name == "PD" & Sampling_point %in% c(11,12)|
      Pond_name == "PF" & Sampling_point %in% c(2,3)|
      Pond_name == "PH" & Sampling_point %in% c(4,5,6,16,17)|
      Pond_name == "TA" & Sampling_point %in% c(1,2)|
      Pond_name == "TF" & Sampling_point %in% c(2:5)|
      Pond_name == "TG" & Sampling_point %in% c(2:6)
  ) %>% 
  #tax_fix() #%>% 
  ps_mutate(Pond_water_condition1 = ifelse(grepl("\\.no.", Crop_species), "Without_fish", "With_fish")) #add a column in metadata named Culture.system where it will look for any patern with "." in the Crop_species column and if it's true, it will mark as Poly
pseq_fish1 # 4245 taxa and 77 samples 

### Alpha diversity----
set.seed(1234)
ps_rarefy_fish1 <- rarefy_even_depth(pseq_fish1, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#687OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_fish1 <- estimate_richness(ps_rarefy_fish1, measures = c("Chao1", "Shannon"))
alpha_estimates_fish1 <- cbind(alpha_estimates_fish1, sample_data(pseq_fish1))

#### Richness----
# Chao1: An estimator of total species richness that accounts for the presence of rare species.
fish1.c <- ggplot(data = alpha_estimates_fish1, aes(y = Chao1, x = Pond_water_condition1)) +
  geom_boxplot(aes(color = Pond_water_condition1)) +
  geom_jitter(aes(color = Pond_water_condition1), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = fish.colors, labels = c("With fish", "Without fish")) + # Customize labels
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 18, margin = margin(r = 5)),
    axis.title.x = element_blank(),#text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank(),
  ) + guides(fill = guide_legend(nrow = 1)) + # Set number of columns in the legend
  labs (color = "Pond water condition",
        y = "Chao1", # Add y-axis title to the plot
        x = "Pond water condition",
        title = "Ponds with and without fish for one month") +
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) + # Insert line break in category names
  ylim(0, 800)
fish1.c

##### Perform Kw test----
kruskal.test(Chao1 ~ Pond_water_condition1, data = alpha_estimates_fish1)
#data:  Chao1 by Crop_species2
#Kruskal-Wallis chi-squared = 0.25266, df = 1, p-value = 0.6152

#### Shannon----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
fish1.s <- ggplot(data = alpha_estimates_fish1, aes(y = Shannon, x = Pond_water_condition1)) +
  geom_boxplot(aes(color = Pond_water_condition1)) +
  geom_jitter(aes(color = Pond_water_condition1), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = fish.colors, labels = c("With fish", "Without fish")) + # Customize labels
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 18, face = "bold", margin = margin(r = 5)),
    axis.title.x = element_blank(),#text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank()
  ) + guides(fill = guide_legend(ncol = 1)) + # Set number of columns in the legend
  labs (color = "Pond water condition",
        y = "Shannon", # Add y-axis title to the plot
        x = "Pond water condition",
        title = "Ponds with and without fish for one month") +
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) + # Insert line break in category names
  ylim(3, 6)
fish1.s

##### Perform Kw test----
kruskal.test(Shannon ~ Pond_water_condition1, data = alpha_estimates_fish1)
#Kruskal-Wallis rank sum test
# Kruskal-Wallis chi-squared = 0.037988, df = 1, p-value = 0.8455

### Beta diversity----
# wunifrac:weighted UniFrac took about #01:35- 04:45 = 12 minutes to calculate the distance
fish1_wuni <- pseq_fish1 %>%
  tax_transform(rank = "unique", trans = "compositional") %>%
  dist_calc(dist = "wunifrac") 
#saveRDS(fish1_wuni, "pond_without_fish_for_1month_weighted_unifrac_20240917.rds")
fish1_wuni <- readRDS("pond_without_fish_for_1month_weighted_unifrac_20240917.rds")

# Culture.system
fish1.we <- fish1_wuni %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Pond_water_condition1", 
    fill = "Pond_water_condition1",
    shape = 21, 
    alpha = 1,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(colour = Pond_water_condition1)
  ) +
  scale_colour_manual(values = fish.colors, labels = c("With fish", "Without fish")) +
  scale_fill_manual(values = fish.colors, labels = c("With fish", "Without fish")) +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
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
  labs(colour = "Pond water condition",
       fill = "Pond water condition",
    #x = "PCoA1 [25.7%]",
    #y = "PCoA2 [16.2%]",
    title = "Ponds with and without fish for one month")
fish1.we

##### Significance Test----
#Plot PERMANOVA with phyloseq
metadata_fish1 <- sample_data(pseq_fish1) %>%
  data.frame() %>%
  tibble()
wuni_dist.fish1 = phyloseq::distance(pseq_fish1, method = "wunifrac")
PERM_fish1_wuni <- adonis2(wuni_dist.fish1 ~ Pond_water_condition1, data = metadata_fish1)

# to add test result on the plot
fish1.we2 <- fish1.we + annotate(geom = "label",
                               label = paste("PERMANOVA: R² = ", round(PERM_fish1_wuni["Pond_water_condition1","R2"], 3), 
                                             ", p = ", PERM_fish1_wuni["Pond_water_condition1", "Pr(>F)"], sep = ""),
                               x=Inf, y=Inf, hjust = 1, vjust = 1)
fish1.we2

## 2 months without fish----
pseq_fish2 <- ps %>%   
  ps_filter(
    Sampling_point %in% c(3,4,5,6), 
    Pond_name %in% c("PD", "PE", "PF")
  ) %>% 
  #tax_fix() #%>% 
  ps_mutate(Pond_water_condition2 = ifelse(grepl("\\.no.", Crop_species), "Without_fish", "With_fish")) #add a column in metadata named Crop_species2 where it will look for any patern with ".no." in the Crop_species column and if it's true, it will mark as No_fish, if false, Fish
pseq_fish2 # 2356 taxa and 36 samples

### Alpha diversity----
set.seed(1234)
ps_rarefy_fish2 <- rarefy_even_depth(pseq_fish2, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#328OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_fish2 <- estimate_richness(ps_rarefy_fish2, measures = c("Chao1", "Shannon"))
alpha_estimates_fish2 <- cbind(alpha_estimates_fish2, sample_data(pseq_fish2))

#### Richness----
# Chao1: An estimator of total species richness that accounts for the presence of rare species.
fish2.c <- ggplot(data = alpha_estimates_fish2, aes(y = Chao1, x = Pond_water_condition2)) +
  geom_boxplot(aes(color = Pond_water_condition2)) +
  geom_jitter(aes(color = Pond_water_condition2), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = fish.colors, labels = c("With fish", "Without fish")) + # Customize labels
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 18, margin = margin(r = 5)),
    axis.title.x = element_blank(),#text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank()
  ) + guides(fill = guide_legend(nrow = 1)) + # Set number of columns in the legend
  labs (color = "Pond water condition",
        y = "Chao1", # Add y-axis title to the plot
        x = "Pond water condition",
        title = "Ponds with and without fish for two months") +
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) + # Insert line break in category names
  ylim(0, 800)
fish2.c

##### Perform Kw test----
kruskal.test(Chao1 ~ Pond_water_condition2, data = alpha_estimates_fish2)
#Kruskal-Wallis rank sum test
#data:  Chao1 by Pond_water_condition2  Kruskal-Wallis chi-squared = 9.6136, df = 1, p-value = 0.001931

#### Shannon----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
fish2.s <- ggplot(data = alpha_estimates_fish2, aes(y = Shannon, x = Pond_water_condition2)) +
  geom_boxplot(aes(color = Pond_water_condition2)) +
  geom_jitter(aes(color = Pond_water_condition2), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = fish.colors, labels = c("With fish", "Without fish")) + # Customize labels
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 18, margin = margin(r = 5)),
    axis.title.x = element_blank(),#text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank()
  ) + guides(fill = guide_legend(nrow = 1)) + # Set number of columns in the legend
  labs (color = "Pond water condition",
        y = "Shannon", # Add y-axis title to the plot
        x = "Pond water condition",
        title = "Ponds with and without fish for two months") +
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) + # Insert line break in category names
ylim(3, 6)
fish2.s

##### Perform Kw test----
kruskal.test(Shannon ~ Pond_water_condition2, data = alpha_estimates_fish2)
#Kruskal-Wallis rank sum test
#data:  Shannon by Pond_water_condition2 Kruskal-Wallis chi-squared = 3.1391, df = 1, p-value = 0.07643

# perform Dunn test (https://www.youtube.com/watch?v=pyLQmUfrel8)
stat.test.c <- dunn_test(Chao1 ~ Pond_water_condition2, data = alpha_estimates_fish2, p.adjust.method = "BH")

# add the value on the plotBox plot
stat.test.c2 <- stat.test.c %>% add_xy_position(x = "Pond_water_condition2")
fish2.c2 <- fish2.c + stat_pvalue_manual(stat.test.c2, label = "p.adj.signif", 
                                       step.increase = 0.05, tip.length = 0.005, 
                                       size = 8,  # Adjust the size of the asterisk here
                                       hide.ns = TRUE)
fish2.c2

### Beta diversity (weighted unifrac)----
# wunifrac:weighted UniFrac took about #01:35- 04:45 = 12 minutes to calculate the distance
fish2_wuni <- pseq_fish2 %>%
  tax_transform(rank = "unique", trans = "compositional") %>%
  dist_calc(dist = "wunifrac") 
#saveRDS(fish2_wuni, "pond_without_fish_for_2month_weighted_unifrac_20240917.rds")
fish2_wuni <- readRDS("pond_without_fish_for_2month_weighted_unifrac_20240917.rds")

# Culture.system
fish2.we <- fish2_wuni %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Pond_water_condition2", 
    fill = "Pond_water_condition2",
    shape = 21, 
    alpha = 1,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(colour = Pond_water_condition2)
  ) +
  scale_colour_manual(values = fish.colors, labels = c("With fish", "Without fish")) +
  scale_fill_manual(values = fish.colors, labels = c("With fish", "Without fish")) +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
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
  labs(colour = "Pond water condition",
       fill = "Pond water condition",
       #x = "PCoA1 [25.7%]",
       #y = "PCoA2 [16.2%]",
       title = "Ponds with and without fish for two months")
fish2.we

#### Significance Test----
#Plot PERMANOVA with phyloseq
metadata_fish2 <- sample_data(pseq_fish2) %>%
  data.frame() %>%
  tibble()
wuni_dist.fish2 = phyloseq::distance(pseq_fish2, method="wunifrac")
PERM_fish2_wuni <- adonis2(wuni_dist.fish2 ~ Pond_water_condition2, data = metadata_fish2)

# to add test result on the plot
fish2.we2 <- fish2.we + annotate(geom = "label",
                           label = paste("PERMANOVA: R² = ", round(PERM_fish2_wuni["Pond_water_condition2","R2"], 3), 
                                         ", p = ", PERM_fish2_wuni["Pond_water_condition2", "Pr(>F)"], sep = ""),
                           x=Inf, y=Inf, hjust = 1, vjust = 1)
fish2.we2

# Microeukaryotes----
# Load the phyloseq object
ps.18s <- readRDS("phyloseq_18S_filtered_with_tree_pr2_90-150bp_20240416.rds")
ps.18s # 5390 taxa and 872 samples

## 1 month without fish----
pseq_fish1.18s <- ps.18s  %>%   
  ps_filter(
    Pond_name == "PA" & Sampling_point %in% c(11,12)|
      Pond_name == "PB" & Sampling_point %in% c(7, 8)|
      Pond_name == "PC" & Sampling_point %in% c(13,14,15)|
      Pond_name == "PD" & Sampling_point %in% c(11,12)|
      Pond_name == "PF" & Sampling_point %in% c(2,3)|
      Pond_name == "PH" & Sampling_point %in% c(4,5,6,16,17)|
      Pond_name == "TA" & Sampling_point %in% c(1,2)|
      Pond_name == "TF" & Sampling_point %in% c(2:5)|
      Pond_name == "TG" & Sampling_point %in% c(2:6)
  ) %>% 
  #tax_fix() #%>% 
  ps_mutate(Pond_water_condition1 = ifelse(grepl("\\.no.", Crop_species), "Without_fish", "With_fish")) #add a column in metadata named Culture.system where it will look for any patern with "." in the Crop_species column and if it's true, it will mark as Poly
pseq_fish1.18s  # 4245 taxa and 77 samples 

### Alpha diversity----
set.seed(1234)
ps_rarefy_fish1.18s  <- rarefy_even_depth(pseq_fish1.18s , rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#687OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_fish1.18s <- estimate_richness(ps_rarefy_fish1.18s , measures = c("Chao1", "Shannon"))
alpha_estimates_fish1.18s <- cbind(alpha_estimates_fish1.18s , sample_data(pseq_fish1.18s))

#### Richness----
# Chao1: An estimator of total species richness that accounts for the presence of rare species.
fish1.c.18s <- ggplot(data = alpha_estimates_fish1.18s, aes(y = Chao1, x = Pond_water_condition1)) +
  geom_boxplot(aes(color = Pond_water_condition1)) +
  geom_jitter(aes(color = Pond_water_condition1), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = fish.colors, labels = c("With fish", "Without fish")) + # Customize labels
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 18, margin = margin(r = 5)),
    axis.title.x = element_blank(),#text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank()
  ) + guides(fill = guide_legend(nrow = 1)) + # Set number of columns in the legend
  labs (color = "Pond water condition",
        y = "Chao1", # Add y-axis title to the plot
        x = "Pond water condition",
        title = "Ponds with and without fish for one month") +
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) + # Insert line break in category names
  ylim(0, 400)
fish1.c.18s

##### Perform Kw test----
kruskal.test(Chao1 ~ Pond_water_condition1, data = alpha_estimates_fish1.18s)
#data:  Chao1 by Crop_species2
#Kruskal-Wallis chi-squared = 1.727, df = 1, p-value = 0.1888

#### Shannon----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
fish1.s.18s <- ggplot(data = alpha_estimates_fish1.18s, aes(y = Shannon, x = Pond_water_condition1)) +
  geom_boxplot(aes(color = Pond_water_condition1)) +
  geom_jitter(aes(color = Pond_water_condition1), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = fish.colors, labels = c("With fish", "Without fish")) + # Customize labels
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 18, margin = margin(r = 5)),
    axis.title.x = element_blank(),#text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank()
  ) + guides(fill = guide_legend(nrow = 1)) + # Set number of columns in the legend
  labs (color = "Pond water condition",
        y = "Shannon", # Add y-axis title to the plot
        x = "Pond water condition",
        title = "Ponds with and without fish for one month") +
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) + # Insert line break in category names
  ylim(0, 6)
fish1.s.18s

##### Perform Kw test----
kruskal.test(Shannon ~ Pond_water_condition1, data = alpha_estimates_fish1.18s)
#Kruskal-Wallis rank sum test
# Kruskal-Wallis chi-squared = 1.1485, df = 1, p-value = 0.2839

### Beta diversity----
# wunifrac:weighted UniFrac took about #01:35- 04:45 = 12 minutes to calculate the distance
fish1_wuni.18s <- pseq_fish1.18s %>%
  tax_transform(rank = "unique", trans = "compositional") %>%
  dist_calc(dist = "wunifrac") 
#saveRDS(fish1_wuni.18s, "pond_without_fish_for_1month_weighted_unifrac.18s_20240917.rds")
fish1_wuni.18s <- readRDS("pond_without_fish_for_1month_weighted_unifrac.18s_20240917.rds")

# Culture.system
fish1.we.18s <- fish1_wuni.18s %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Pond_water_condition1", 
    fill = "Pond_water_condition1",
    shape = 21, 
    alpha = 1,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(colour = Pond_water_condition1)
  ) +
  scale_colour_manual(values = fish.colors, labels = c("With fish", "Without fish")) +
  scale_fill_manual(values = fish.colors, labels = c("With fish", "Without fish")) +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
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
  labs(colour = "Pond water condition",
       fill = "Pond water condition",
       #x = "PCoA1 [25.7%]",
       #y = "PCoA2 [16.2%]",
       title = "Ponds with and without fish for one month")
fish1.we.18s

#### Significance Test----
#Plot PERMANOVA with phyloseq
metadata_fish1.18s <- sample_data(pseq_fish1.18s) %>%
  data.frame() %>%
  tibble()
wuni_dist.fish1.18s = phyloseq::distance(pseq_fish1.18s, method = "wunifrac")
PERM_fish1_wuni.18s <- adonis2(wuni_dist.fish1.18s ~ Pond_water_condition1, data = metadata_fish1.18s)

# to add test result on the plot
fish1.we2.18s <- fish1.we.18s + annotate(geom = "label",
                                         label = paste("PERMANOVA: R² = ", round(PERM_fish1_wuni.18s["Pond_water_condition1","R2"], 3), 
                                                       ", p = ", PERM_fish1_wuni.18s["Pond_water_condition1", "Pr(>F)"], sep = ""),
                                         x=Inf, y=Inf, hjust = 1, vjust = 1)
fish1.we2.18s

## 2 months without fish----
pseq_fish2.18s <- ps.18s %>%   
  ps_filter(
    Sampling_point %in% c(3,4,5,6), 
    Pond_name %in% c("PD", "PE", "PF")
  ) %>% 
  #tax_fix() #%>% 
  ps_mutate(Pond_water_condition2 = ifelse(grepl("\\.no.", Crop_species), "Without_fish", "With_fish")) #add a column in metadata named Crop_species2 where it will look for any patern with ".no." in the Crop_species column and if it's true, it will mark as No_fish, if false, Fish
pseq_fish2.18s # 2356 taxa and 36 samples

## Alpha diversity----
set.seed(1234)
ps_rarefy_fish2.18s <- rarefy_even_depth(pseq_fish2.18s, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#328OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_fish2.18s <- estimate_richness(ps_rarefy_fish2.18s, measures = c("Chao1", "Shannon"))
alpha_estimates_fish2.18s <- cbind(alpha_estimates_fish2.18s, sample_data(pseq_fish2.18s))

#### Richness----
# Chao1: An estimator of total species richness that accounts for the presence of rare species.
fish2.c.18s <- ggplot(data = alpha_estimates_fish2.18s, aes(y = Chao1, x = Pond_water_condition2)) +
  geom_boxplot(aes(color = Pond_water_condition2)) +
  geom_jitter(aes(color = Pond_water_condition2), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = fish.colors, labels = c("With fish", "Without fish")) + # Customize labels
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 18, margin = margin(r = 5)),
    axis.title.x = element_blank(),#text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank()
  ) + guides(fill = guide_legend(nrow = 1)) + # Set number of columns in the legend
  labs (color = "Pond water condition",
        y = "Chao1", # Add y-axis title to the plot
        x = "Pond water condition",
        title = "Ponds with and without fish for two months") +
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) + # Insert line break in category names
  ylim(0, 400)
fish2.c.18s

##### Perform Kw test----
kruskal.test(Chao1 ~ Pond_water_condition2, data = alpha_estimates_fish2.18s)
#Kruskal-Wallis rank sum test
#data:  Chao1 by Pond_water_condition2  Kruskal-Wallis chi-squared = 0.0012976, df = 1, p-value = 0.9713

### Shannon----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
fish2.s.18s <- ggplot(data = alpha_estimates_fish2.18s, aes(y = Shannon, x = Pond_water_condition2)) +
  geom_boxplot(aes(color = Pond_water_condition2)) +
  geom_jitter(aes(color = Pond_water_condition2), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = fish.colors, labels = c("With fish", "Without fish")) + # Customize labels
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 18, margin = margin(r = 5)),
    axis.title.x = element_blank(),#text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank()
  ) + guides(fill = guide_legend(nrow = 1)) + # Set number of columns in the legend
  labs (color = "Pond water condition",
        y = "Shannon", # Add y-axis title to the plot
        x = "Pond water condition",
        title = "Ponds with and without fish for two months") +
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) + # Insert line break in category names
  ylim(0, 6)
fish2.s.18s

##### Perform Kw test----
kruskal.test(Shannon ~ Pond_water_condition2, data = alpha_estimates_fish2.18s)
#Kruskal-Wallis rank sum test
#data:  Shannon by Pond_water_condition2 Kruskal-Wallis chi-squared = 0.063581, df = 1, p-value = 0.8009

### Beta diversity----
# wunifrac:weighted UniFrac took about #01:35- 04:45 = 12 minutes to calculate the distance
fish2_wuni.18s <- pseq_fish2.18s %>%
  tax_transform(rank = "unique", trans = "compositional") %>%
  dist_calc(dist = "wunifrac") 
#saveRDS(fish2_wuni.18s, "pond_without_fish_for_2month_weighted_unifrac.18s_20240917.rds")
fish2_wuni.18s <- readRDS("pond_without_fish_for_2month_weighted_unifrac.18s_20240917.rds")

# Culture.system
fish2.we.18s <- fish2_wuni.18s %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Pond_water_condition2", 
    fill = "Pond_water_condition2",
    shape = 21, 
    alpha = 1,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(colour = Pond_water_condition2)
  ) +
  scale_colour_manual(values = fish.colors, labels = c("With fish", "Without fish")) +
  scale_fill_manual(values = fish.colors, labels = c("With fish", "Without fish")) +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
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
  labs(colour = "Pond water condition",
       fill = "Pond water condition",
       #x = "PCoA1 [25.7%]",
       #y = "PCoA2 [16.2%]",
       title = "Ponds with and without fish for two months")
fish2.we.18s

#### Significance Test----
#Plot PERMANOVA with phyloseq
metadata_fish2.18s <- sample_data(pseq_fish2.18s) %>%
  data.frame() %>%
  tibble()
wuni_dist.fish2.18s = phyloseq::distance(pseq_fish2.18s, method="wunifrac")
PERM_fish2_wuni.18s <- adonis2(wuni_dist.fish2.18s ~ Pond_water_condition2, data = metadata_fish2.18s)

# to add test result on the plot
fish2.we2.18s <- fish2.we.18s + annotate(geom = "label",
                                         label = paste("PERMANOVA: R² = ", round(PERM_fish2_wuni.18s["Pond_water_condition2","R2"], 3), 
                                                       ", p = ", PERM_fish2_wuni.18s["Pond_water_condition2", "Pr(>F)"], sep = ""),
                                         x=Inf, y=Inf, hjust = 1, vjust = 1)
fish2.we2.18s

# Combine plots----
cowplot::plot_grid(fish1.c + theme(legend.position = "none", plot.title = element_blank(), axis.text.x = element_blank()),
                   fish1.s + theme(legend.position = "none", plot.title = element_blank(), axis.text.x = element_blank()),
                   fish1.we2 + theme(legend.position = "none", plot.title = element_blank()),
                   fish1.c.18s + theme(legend.position = "none", plot.title = element_blank(), axis.text.x = element_blank()),
                   fish1.s.18s + theme(legend.position = "none", plot.title = element_blank(), axis.text.x = element_blank()),
                   fish1.we2.18s + theme(legend.position = "none", plot.title = element_blank()),
                   fish2.c2 + theme(legend.position = "none", plot.title = element_blank(), axis.text.x = element_blank()),
                   fish2.s + theme(legend.position = "none", plot.title = element_blank(), axis.text.x = element_blank()),
                   fish2.we2 + theme(legend.position = "none", plot.title = element_blank()),
                   fish2.c.18s + theme(legend.position = "none", plot.title = element_blank()),
                   fish2.s.18s + theme(legend.position = "none", plot.title = element_blank()),
                   fish2.we2.18s + theme(legend.position = "none", plot.title = element_blank()),
                   rel_heights = c(0.9,0.9,0.9,1,0.9,0.9,0.9,1,0.9,0.9, 0.9,1),  # Adjust relative widths of the plots
                   labels = c("A", "", "", "B","", "", "C","", "", "D"),  # Adjust plot labels
                   align = "v",
                   nrow = 4)

# save as 2000*1500





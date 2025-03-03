#15/02/2025
# Sanjit Debnath
# This script is to plot water physicochemistry and bacterial composition in pangasius and tilapia ponds.
# Microbial composition with a relative abundance >= 1% was plotted
# If I use tax_filter(min_prevalence = 0.01), it's a bit difficult to make as I want. so from the excel file where I have the microbial composition, I checked how many taxa are above RA 1 and plotted them.

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
library(paletteer); packageVersion("paletteer")
library(writexl); packageVersion("writexl")
library(phylosmith);packageVersion("phylosmith")
library(dplyr); packageVersion("dplyr")

## Setup working dictionary first
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/statistical_analysis")

# Set theme for ggplot
theme_set(theme_bw())
set.seed(1234)

# Prokaryotes----
# Load the phyloseq object
ps <- readRDS("phyloseq_metadata_6_v3_20240419.rds")
ps # 10523 taxa and 891 samples

# Replace taxonomy names at the phylum level with new names
tax_table(ps)[, "Phylum"] <- gsub("Proteobacteria", "Pseudomonadota", tax_table(ps)[, "Phylum"])
tax_table(ps)[, "Phylum"] <- gsub("Firmicutes", "Bacillota", tax_table(ps)[, "Phylum"])
tax_table(ps)[, "Phylum"] <- gsub("Actinobacteriota", "Actinomycetota", tax_table(ps)[, "Phylum"])

# Check if changes were applied
table(tax_table(ps)[, "Phylum"])

## Convert ppm to ppt for salinity
# Extract sample metadata
sample_data(ps)$Salinity_ppt <- sample_data(ps)$Salinity_ppm / 1000

# Check if the new column has been added correctly
head(sample_data(ps))


# Set season color
#season.colors <- c("Monsoon" = "#EC823C","Winter" ="#1B9E77","Pre-monsoon" = "#100AFF")
#season.colors2 <- c("Monsoon" = "#199C75","Winter" ="#D95E00","Pre-monsoon" = "#736EB3")
#season.colors3 <- c("Monsoon1" = "#199C75","Winter1" ="#D95E00","Pre_monsoon1" = "#736EB3", "Monsoon2" = "#E6288A","Winter2" ="#66A61C","Pre_monsoon2" = "#E6AB00")

# 1. Plot the parameters----
# I think instead of plotting pond wise, if i plot mean of all pangasius and tilapia ponds and plot them, coloured by season, it makes the plot easier to understand and interpret

## Pangasius ponds----
pang.ep <- ps %>% 
  ps_filter(Crop == "Pangasius")
pang.ep #  8334 taxa and 451 samples

# Define colors for each upazila
upazila.colors <- c("Jamalpur_Sadar" = "#1B9E77FF", "Muktagacha" = "#E6AB02FF", "Tarakanda" = "#490092FF")
sample.colors <- c("Monsoon" = "#199C75","Winter" ="#D95E00","Pre-monsoon" = "#736EB3",
                   "Pangasius" = "#9DCC00", "Tilapia" = "#006DDBFF")

# Make a dataframe (to plot, I need my data as a dataframe and the variables for x and y axis, besides choose a varible to make fill or color of the plot. so need to make a dataframe from the phyloseq object with the sample data)
df.p <- sample_data(pang.ep) %>%     # %>% works as a pipe that allows the output of a previous command to be used as input to another command instead of using nested functions.
  data.frame() 

# Rename season name
df.p <- df.p %>%
  mutate(Season = recode(Season,
                         "01.Monsoon" = "Monsoon",
                         "02.Winter" = "Winter",
                         "03.Pre.Monsoon" = "Pre-monsoon",
                         "04.Monsoon" = "Monsoon",
                         "05.Winter" = "Winter",
                         "06.Pre.Monsoon" = "Pre-monsoon",
                         "06.Pre.monsoon" = "Pre-monsoon"))
# Manually define the correct order of months
month_order <- c("Jun-16","Oct-16","Nov-16","Dec-16","Jan-17","Feb-17","Mar-17",
                 "Apr-17","May-17","Jun-17","Jul-17","Aug-17","Sep-17","Oct-17", 
                 "Nov-17","Dec-17","Jan-18","Feb-18","Mar-18","Apr-18","May-18")

# Convert Sampling_months to a factor with the specified levels
df.p$Sampling_months <- factor(df.p$Sampling_months, levels = month_order)

# Calculate the mean of the parameters and Summarize the data
pang_df_summary <- df.p %>%
  group_by(Crop, Sampling_months, Season) %>%
  summarize(
    mean_salinity = mean(Salinity_ppt, na.rm = TRUE),
    sd_salinity = sd(Salinity_ppt, na.rm = TRUE),
    mean_temperature = mean(Temperature_C, na.rm = TRUE),
    sd_temperature = sd(Temperature_C, na.rm = TRUE),
    mean_pH = mean(pH, na.rm = TRUE),
    sd_pH = sd(pH, na.rm = TRUE),
    mean_DO = mean(DO_mg_l, na.rm = TRUE),
    sd_DO = sd(DO_mg_l, na.rm = TRUE)
  )
# Save the summary as an Excel file
#write_xlsx(pang_df_summary, "pangasius_parameter_summary_20250215.xlsx")

# Line plot for salinity
p1 <- ggplot(pang_df_summary, aes(x = Sampling_months, y = mean_salinity, 
                                  color = Season, group = Crop)) + 
  geom_point(size = 3) +
  geom_line() +
  scale_colour_manual(values = sample.colors) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_blank(),#text(angle = 90, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(),#text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank()) + 
  labs(color = "Season", # Custom legend name
       y = "Salinity (ppt)",  # Add y-axis label
       x = "Sampling month",   # Add x-axis label
       title = "Pangasius ponds") + # to add title
  ylim(0.09, 0.6)
p1

# Line plot for temperature
p2 <- ggplot(pang_df_summary, aes(x = Sampling_months, y = mean_temperature, 
                                  color = Season, group = Crop)) + 
  geom_point(size = 3) +
  geom_line() +
  scale_colour_manual(values = sample.colors) +
  theme(
    plot.title = element_blank(),#element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_blank(),#element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(),#element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank()) + 
  labs(color = "Season", # Custom legend name
       y = "Temp (°C)",  # Add y-axis label
       x = "Sampling month",   # Add x-axis label
       title = "Pangasius ponds") +# to add title
  ylim(20, 35)
p2

# Line plot for pH
p3 <- ggplot(pang_df_summary, aes(x = Sampling_months, y = mean_pH, 
                                  color = Season, group = Crop)) + 
  geom_point(size = 3) +
  geom_line() +
  scale_colour_manual(values = sample.colors) +
  theme(
    plot.title = element_blank(),#element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_blank(),#element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(),#element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank()) + 
  labs(color = "Season", # Custom legend name
       y = "pH",  # Add y-axis label
       x = "Sampling month",   # Add x-axis label
       title = "Pangasius ponds") +
  ylim(7,9)
p3

# Line plot for DO
p4 <- ggplot(pang_df_summary, aes(x = Sampling_months, y = mean_DO, 
                                  color = Season, group = Crop)) + 
  geom_point(size = 3) +
  geom_line() +
  scale_colour_manual(values = sample.colors) +
  theme(
    plot.title = element_blank(),#element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_blank(),#element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(),#element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank()) + 
  labs(color = "Season", # Custom legend name
       y = "DO (mg/L)",  # Add y-axis label
       x = "Sampling month",   # Add x-axis label
       title = "Pangasius ponds") +
  ylim(4, 15)
p4

### Combine plots with all parameters----
comb.p.ep <- cowplot::plot_grid(p1 + theme(legend.position = "none"),# axis.title.x = element_blank()), 
                                p2 + theme(legend.position = "none"),# axis.title.x = element_blank()), 
                                p3 + theme(legend.position = "none"),# axis.title.x = element_blank()), 
                                p4 + theme(legend.position = "none"), #axis.title.x = element_blank()), 
                                rel_heights = c(0.7, 0.6, 0.6, 0.6),
                                labels = "A",
                                align = "v",
                                ncol = 1)
comb.p.ep


## Tilapia ponds----
tila.ep <- ps %>% 
  ps_filter(Crop == "Tilapia")
tila.ep # 8624 taxa and 440 samples

# Make a dataframe
df.t <- sample_data(tila.ep) %>%     # %>% works as a pipe that allows the output of a previous command to be used as input to another command instead of using nested functions.
  data.frame() 

# Rename season name
df.t <- df.t %>%
  mutate(Season = recode(Season,
                         "01.Monsoon" = "Monsoon",
                         "02.Winter" = "Winter",
                         "03.Pre.Monsoon" = "Pre-monsoon",
                         "04.Monsoon" = "Monsoon",
                         "05.Winter" = "Winter",
                         "06.Pre.Monsoon" = "Pre-monsoon",
                         "06.Pre.monsoon" = "Pre-monsoon"))

# Manually define the correct order of months
month_order <- c("Jun-16","Oct-16","Nov-16","Dec-16","Jan-17","Feb-17","Mar-17",
                 "Apr-17","May-17","Jun-17","Jul-17","Aug-17","Sep-17","Oct-17", 
                 "Nov-17","Dec-17","Jan-18","Feb-18","Mar-18","Apr-18","May-18")

# Convert Sampling_months to a factor with the specified levels
df.t$Sampling_months <- factor(df.t$Sampling_months, levels = month_order)

# Calculate the mean of the parameters and Summarize the data
tila_df_summary <- df.t %>%
  group_by(Crop, Sampling_months, Season) %>%
  summarize(
    mean_salinity = mean(Salinity_ppt, na.rm = TRUE),
    sd_salinity = sd(Salinity_ppt, na.rm = TRUE),
    mean_temperature = mean(Temperature_C, na.rm = TRUE),
    sd_temperature = sd(Temperature_C, na.rm = TRUE),
    mean_pH = mean(pH, na.rm = TRUE),
    sd_pH = sd(pH, na.rm = TRUE),
    mean_DO = mean(DO_mg_l, na.rm = TRUE),
    sd_DO = sd(DO_mg_l, na.rm = TRUE))
# Save the summary as an Excel file
#write_xlsx(tila_df_summary, "tilapia_parameter_summary_20250215.xlsx")

# Line plot for salinity
t1 <- ggplot(tila_df_summary, aes(x = Sampling_months, y = mean_salinity, 
                                  color = Season, group = Crop)) + 
  geom_point(size = 3) +
  geom_line() +
  scale_colour_manual(values = sample.colors) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_blank(),#element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(),#element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank()) + 
  labs(color = "Season", # Custom legend name
       y = "Salinity (ppt)",  # Add y-axis label
       x = "Sampling month",   # Add x-axis label
       title = "Tilapia ponds") +# to add title
  ylim(0.09, 0.6)
t1

# Line plot for temperature
t2 <- ggplot(tila_df_summary, aes(x = Sampling_months, y = mean_temperature, 
                                  color = Season, group = Crop)) + 
  geom_point(size = 3) +
  geom_line() +
  scale_colour_manual(values = sample.colors) +
  theme(
    plot.title = element_blank(),#element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_blank(),#element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(),#element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank()) + 
  labs(color = "Season", # Custom legend name
       y = "Temp (°C)",  # Add y-axis label
       x = "Sampling month",   # Add x-axis label
       title = "Tilapia ponds") +# to add title
  ylim(20, 35)
t2

# Line plot for pH
t3 <- ggplot(tila_df_summary, aes(x = Sampling_months, y = mean_pH, 
                                  color = Season, group = Crop)) + 
  geom_point(size = 3) +
  geom_line() +
  scale_colour_manual(values = sample.colors) +
  theme(
    plot.title = element_blank(),#element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_blank(),#element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(),#element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank()) + 
  labs(color = "Season", # Custom legend name
       y = "pH",  # Add y-axis label
       x = "Sampling month",   # Add x-axis label
       title = "Tilapia ponds") +
  ylim(7,9)
t3

# Line plot for DO
t4 <- ggplot(tila_df_summary, aes(x = Sampling_months, y = mean_DO, 
                                  color = Season, group = Crop)) + 
  geom_point(size = 3) +
  geom_line() +
  scale_colour_manual(values = sample.colors) +
  theme(
    plot.title = element_blank(),#element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_blank(),#element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(),#element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank()) + 
  labs(color = "Season", # Custom legend name
       y = "DO (mg/L)",  # Add y-axis label
       x = "Sampling month",   # Add x-axis label
       title = "Tilapia ponds") +
  ylim(4, 15)
t4

### Combine water physicochemistry plots----
comb.t.ep <- cowplot::plot_grid(t1 + theme(legend.position = "none"), #axis.title.x = element_blank()), 
                                t2 + theme(legend.position = "none"), #axis.title.x = element_blank()), 
                                t3 + theme(legend.position = "none"),# axis.title.x = element_blank()), 
                                t4 + theme(legend.position = "none"), #axis.title.x = element_blank()), 
                                rel_heights = c(0.7, 0.6, 0.6, 0.6),
                                labels = "B",
                                align = "v",
                                ncol = 1)
comb.t.ep

# 2. Microbial composition barplots----
# Manually define the correct order of months
month_order <- c("Jun-16","Oct-16","Nov-16","Dec-16","Jan-17","Feb-17","Mar-17",
                 "Apr-17","May-17","Jun-17","Jul-17","Aug-17","Sep-17","Oct-17", 
                 "Nov-17","Dec-17","Jan-18","Feb-18","Mar-18","Apr-18","May-18")

# Access the sample data from the phyloseq object and modify the Sampling_months column
sample_data(ps)$Sampling_months <- factor(sample_data(ps)$Sampling_months, levels = month_order)

## Subset pangasius----
pseq_pang <- ps %>% 
  ps_filter(
    Crop == "Pangasius",
    Crop_species != "Shing"
  ) %>% tax_fix()
pseq_pang # 8091 taxa and 421 samples

### Phylum (my palette)----
myPal.phy.16s2 <- tax_palette(
  data = pseq_pang, rank = "Phylum", n = 20, pal = "kelly",
  add = c(Other = "lightgrey"))
myPal.phy.16s2["Bacillota"] <- "#00D9FF" # Override existing color if any color is not good
myPal.phy.16s2["Cyanobacteria"] <- "#ABE496"
#myPal.phy.16s2["Cyanobacteria"] <- "#ABE496"
tax_palette_plot(myPal.phy.16s2)

### Pangasius pond_phylum----
topTaxa.pang.phy <- pseq_pang %>%
  tax_top(n = 10, rank = "Phylum") %>% # >=1%
  sort() # this makes them alphabetical

## plot with alphabetical sorting
p2.phy <- pseq_pang %>%
  tax_sort(by = sum, at = "Phylum") %>% # this orders all genera by abundance
  ps_select(Crop_species, Sampling_months) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Sampling_months") %>%
  comp_barplot(tax_order = topTaxa.pang.phy, # this brings the named taxa to the front
               tax_level = "Phylum", n_taxa = 10, # RA >= 1%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.phy.16s2) +
  #coord_flip() +
  theme(plot.title = element_blank(),#text(hjust = 0.5, size = 18, margin = margin(b = 15)), # Increase margin for plot title
        axis.text.x = element_blank(),#text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_blank(),  # to remove x axis title
        axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
        legend.position = "right",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18)) +
  labs(title = "Pangasius ponds",
       y = "Relative abundance") +
  guides(fill = guide_legend(ncol = 1))  # Adjust the number of rows in the legend
p2.phy

## Subset tilapia----
pseq_tila <- ps %>% 
  ps_filter(
    Crop == "Tilapia",
    Crop_species != "Gulsha.Carp", 
    Crop_species != "Gulsha.Pabda"
  ) %>% tax_fix()
pseq_tila # 8380 taxa and 391 samples

### Tilapia pond_phylum----
topTaxa.tila.phy <- pseq_tila %>%
  tax_top(n = 10, rank = "Phylum") %>% # >=1%
  sort() # this makes them alphabetical

## plot with alphabetical sorting
t2.phy <- pseq_tila %>%
  tax_sort(by = sum, at = "Phylum") %>% # this orders all genera by abundance
  ps_select(Crop_species, Sampling_months) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Sampling_months") %>%
  comp_barplot(tax_order = topTaxa.tila.phy, # this brings the named taxa to the front
               tax_level = "Phylum", n_taxa = 10, # RA > 1%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.phy.16s2) +
  #coord_flip() +
  theme(plot.title = element_blank(),#text(hjust = 0.5, size = 18, margin = margin(b = 15)), # Increase margin for plot title
        axis.text.x = element_blank(),#text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_blank(),  # to remove x axis title
        axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
        legend.position = "right",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18)) +
  labs(title = "Tilapia ponds",
       y = "Relative abundance") +
  guides(fill = guide_legend(ncol = 1))  # Adjust the number of rows in the legend
t2.phy

## Family----
# my palette
myPal.fam.16s <- tax_palette(
  data = pseq_pang, rank = "Family", n = 35, pal = "brewerPlus",
  add = c(Other = "lightgrey"))
#myPal.fam.16s["Sporichthyaceae"] <- "#00D9FF" # Override existing color if any color is not good
#myPal.fam.16s["Ilumatobacteraceae"] <- "#ABE496"
myPal.fam.16s["Microbacteriaceae"] <- "#EA70FF" 
myPal.fam.16s["Mycobacteriaceae"] <-  "#5f7b35"
myPal.fam.16s["Chthoniobacteraceae"] <- "#CC00A7"
tax_palette_plot(myPal.fam.16s)

### Pangasius pond----
topTaxa.pang.fam <- pseq_pang %>%
  tax_top(n = 22, rank = "Family") %>% # >= 1%
  sort() # this makes them alphabetical

## plot with alphabetical sorting
p2.fam <- pseq_pang %>%
  tax_sort(by = sum, at = "Family") %>% # this orders all genera by abundance
  ps_select(Crop_species, Sampling_months) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Sampling_months") %>%
  comp_barplot(tax_order = topTaxa.pang.fam, # this brings the named taxa to the front
               tax_level = "Family", n_taxa = 22, # RA >= 1%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.fam.16s) +
  #coord_flip() +
  theme(plot.title = element_blank(),#text(hjust = 0.5, size = 18, margin = margin(b = 15)), # Increase margin for plot title
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
p2.fam

### Tilapia pond----
topTaxa.tila.fam <- pseq_tila %>%
  tax_top(n = 25, rank = "Family") %>% # >=1%
  sort() # this makes them alphabetical

## plot with alphabetical sorting
t2.fam <- pseq_tila %>%
  tax_sort(by = sum, at = "Family") %>% # this orders all genera by abundance
  ps_select(Crop_species, Sampling_months) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Sampling_months") %>%
  comp_barplot(tax_order = topTaxa.tila.fam, # this brings the named taxa to the front
               tax_level = "Family", n_taxa = 25, # RA >= 1%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.fam.16s
  ) +
  #coord_flip() +
  theme(plot.title = element_blank(),#text(hjust = 0.5, size = 18, margin = margin(b = 15)), # Increase margin for plot title
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
t2.fam


# Add season variables as require
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
                                             plot_type = "barerrorbar", errorbar_addpoint = FALSE, 
                                             errorbar_color_black = TRUE, plot_SE = TRUE,
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
  labs(fill = "Season",
       colour = "Season") 
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
                                             plot_type = "barerrorbar", errorbar_addpoint = FALSE, 
                                             errorbar_color_black = TRUE, plot_SE = TRUE,
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
  labs(fill = "Season", colour = "Season")
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
                                             plot_type = "barerrorbar", errorbar_addpoint = FALSE, 
                                             errorbar_color_black = TRUE, plot_SE = TRUE,
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
  labs(fill = "Season", colour = "Season")
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
                                             plot_type = "barerrorbar", errorbar_addpoint = FALSE, 
                                             errorbar_color_black = TRUE, plot_SE = TRUE,
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
  labs(fill = "Season", colour = "Season")
tila.fam.bar

#### Combine prokaryotic barplot----
comb.bars <- cowplot::plot_grid(p2.phy + theme(legend.position = "none"),#+ guides(fill = guide_legend(ncol = 1)),
                                t2.phy + theme(legend.position = "none"),#+ guides(fill = guide_legend(ncol = 1)),
                                p2.fam + theme(legend.position = "none"),#+ guides(fill = guide_legend(ncol = 1)),
                                t2.fam + theme(legend.position = "none"),#+ guides(fill = guide_legend(ncol = 1)),
                                rel_heights = c(0.6, 0.7, 0.6, 0.7),
                                labels = c("C", "D", "E", "F"),
                                ncol = 2)
comb.bars

# Combine environmental parameters and micorbial composition together----
# Combine environmental parameter
ep.pt <- cowplot::plot_grid(comb.p.ep,
                            comb.t.ep,
                            ncol = 2)
ep.pt

# Combine prokaryotic barplots plots----
season.bar.16s <- cowplot::plot_grid(pang.phy.bar + theme(legend.position = "none"),
                                     tila.phy.bar + theme(legend.position = "none"),
                                     pang.fam.bar + theme(legend.position = "none"),
                                     tila.fam.bar + theme(legend.position = "none"),
                                     labels = c("G", "H", "I", "J"))
season.bar.16s

# Combine all
comb1 <- cowplot::plot_grid(ep.pt,
                            comb.bars,
                            #rel_heights = c(0.7, 1),
                            align = "v",
                            ncol = 1)
comb1

# Add legend for now, need to update in inkscape
legend.fam <- get_legend(p2.fam)
comb1.leg <- cowplot::plot_grid(comb1,
                                legend.fam,
                                ncol = 2,
                                rel_widths = c(2,0.43))
comb1.leg 

# Final figure
fig2 <- cowplot::plot_grid(comb1.leg,
                           season.bar.16s,
                           ncol = 1,
                           rel_heights = c(6,4))
fig2 # Save as 2200*1800 as final figure 


# Combinedly plot all legends
legend.p.ep <- get_legend(p1)
legend.t.ep <- get_legend(t1)
legend.p2.p <- get_legend(p2.phy)
legend.t2.p <- get_legend(t2.phy)
legend.p2.f <- get_legend(p2.fam)
legend.t2.f <- get_legend(t2.fam)

comb.legend <- cowplot::plot_grid(legend.p.ep,
                                  legend.t.ep,
                                  legend.p2.p,# + theme(legend.position = "right")+ guides(fill = guide_legend(ncol = 1)),
                                  legend.t2.p,# + theme(legend.position = "right")+ guides(fill = guide_legend(ncol = 1)),
                                  legend.p2.f,#p2.fam + theme(legend.position = "right")+ guides(fill = guide_legend(ncol = 1)),
                                  legend.t2.f,#t2.fam + theme(legend.position = "right")+ guides(fill = guide_legend(ncol = 1)),
                                  #labels = "AUTO",
                                  ncol = 2)
comb.legend # Save as 2100*2500 and then edit the legend to add all the missing together and add this legend to the previous plot




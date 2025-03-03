# 12/09/2024
# Sanjit Debnath
# This script is to plot correlation between water physicochemistry and alpha diversity for microeukaryotes

# Load Libraries
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(microViz); packageVersion("microViz")
library(vegan); packageVersion("vegan") # needed for PERMANOVA test
library(tidyverse); packageVersion("tidyverse")
library(ggpubr); packageVersion("ggpubr")

## Setup working dictionary first
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/statistical_analysis")

# Set theme for ggplot
theme_set(theme_bw())
set.seed(1234)

# Load the phyloseq object
ps.18s <- readRDS("phyloseq_18S_filtered_with_tree_pr2_90-150bp_20240416.rds")
ps.18s # 5390 taxa and 872 samples

## Convert ppm to ppt for salinity
# Extract sample metadata
sample_data(ps.18s)$Salinity_ppt_avg <- sample_data(ps.18s)$Salinity_ppm_avg / 1000

# Check if the new column has been added correctly
head(sample_data(ps.18s))

# Pangasius----
### For the period where the phys chem data were collected only – plot the microbiomes (alpha diversity) against the change in the different phys chem  parameters. 
# first 6 months, no physoochemical parameters available
# subset pangasius pond
ps.pang.18s <- ps.18s %>% 
  ps_filter(
    Sampling_point > 6,
    Crop != "Tilapia" & Crop_species != "Shing")
ps.pang.18s # 7363 taxa and 291 samples

## 1. Correlation with alpha diversity----
### Estimate alpha diversity----
ps_rarefy_7.18s <- rarefy_even_depth(ps.pang.18s, rngseed = 1234) #subsample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#263OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_7.18s <- estimate_richness(ps_rarefy_7.18s, measures = c("Chao1", "Shannon"))
alpha_estimates_7.18s <- cbind(alpha_estimates_7.18s, sample_data(ps.pang.18s))

### Test to decide preason or spearman----
# Check normality for Shannon diversity
shapiro.test(alpha_estimates_7.18s$Shannon) # p-value  = 6.347e-09
# Check normality for Salinity
shapiro.test(alpha_estimates_7.18s$Salinity_ppt_avg) #  p-value < 2.2e-16
shapiro.test(alpha_estimates_7.18s$Temperature_C_avg) # p-value = 1.814e-06
shapiro.test(alpha_estimates_7.18s$pH_avg) # p-value = 3.466e-09
shapiro.test(alpha_estimates_7.18s$DO_mg_l_avg) # p-value = 2.003e-06
# as p < 0.05, these are not normally distributed. So shouldn't do preason
# Spearman’s correlation is a non-parametric test that does not require normality for either variable. So doing spearman is not any problem

# Linear regression (or Correlation?) to examine the relationship between environmental factors and microbial diversity?
### Plot alpha diversity----
#### Salinity----
p1_p.18s <- alpha_estimates_7.18s %>% 
  ggscatter(y = "Shannon", x = "Salinity_ppt_avg",
            ylab = "Shannon", xlab = "Salinity (ppt)",
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, 
            legend = "none",
            #conf.int.level = 0.95,
            cor.method = "spearman",
            label.x.npc = "left", 
            label.y.npc = "top", # Adjust the position of the labels here
            size = 2) +
  geom_point(aes(color = Salinity_ppm_avg), alpha = 1) +
  #gradient_color(c("blue", "red")) +
  scale_color_viridis_c() +  # Use viridis color scale
  ggtitle("A Pangasius ponds") +
  ylim(0,6) +  theme(plot.title = element_blank(),#element_text(size = 18, hjust = 0, face = "bold"),
                     axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
                     axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
                     axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
                     axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)))  # Increase margin for y-axis title)#, face = "bold") # Adjust text size, horizontal alignment, and style
p1_p.18s

#### Temperature----
p2_p.18s <- alpha_estimates_7.18s %>% 
  ggscatter(y = "Shannon", x = "Temperature_C_avg",
            ylab = "Shannon", xlab = "Temperature (°C)",
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE,
            legend = "none",
            cor.method = "spearman",
            label.x.npc = "left", 
            label.y.npc = "top", # Adjust the position of the labels here
            size = 2) +
  geom_point(aes(color = Temperature_C_avg)) +
  #scale_color_gradient(low = "blue", high = "red") +
  scale_color_viridis_c() +  # Use viridis color scale
  ggtitle("A Pangasius ponds") +
  ylim(0,6)+ 
  theme(plot.title = element_blank(),#element_text(size = 18, hjust = 0, face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
        axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)))  # Increase margin for y-axis title)#, face = "bold") # Adjust text size, horizontal alignment, and style
p2_p.18s
     
#### pH----
p3_p.18s <- alpha_estimates_7.18s %>% 
  ggscatter(y = "Shannon", x = "pH_avg",
            add = "reg.line", conf.int = TRUE,
            ylab = "Shannon", xlab = "pH",
            cor.coef = TRUE, 
            legend = "none",
            cor.method = "spearman",
            size = 2) +
  geom_point(aes(color = pH_avg)) +
  #scale_color_gradient(low = "blue", high = "red") +
  scale_color_viridis_c() +  # Use viridis color scale
  ggtitle("A Pangasius pond") +
  ylim(0, 6) + 
  theme(plot.title = element_blank(),#element_text(size = 18, hjust = 0, face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
        axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)))  # Increase margin for y-axis title)#, face = "bold") # Adjust text size, horizontal alignment, and style
p3_p.18s

#### DO----
p4_p.18s <- alpha_estimates_7.18s %>% 
  ggscatter(y = "Shannon", x = "DO_mg_l_avg",
            ylab = "Shannon", xlab = " DO (mg/L)",
            legend = "none",
            add = "reg.line", 
            conf.int = TRUE, 
            cor.coef = TRUE, 
            cor.method = "spearman",
            size = 2) +
  geom_point(aes(color = DO_mg_l_avg)) +
  #scale_color_gradient(low = "blue", high = "red") +
  scale_color_viridis_c() +  # Use viridis color scale
  ggtitle("A Pangasius ponds") +
  ylim(0,6) + 
  theme(plot.title = element_blank(),#element_text(size = 18, hjust = 0, face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
        axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)))  # Increase margin for y-axis title)#, face = "bold") # Adjust text size, horizontal alignment, and style
p4_p.18s

## 2.RDA plot----
## rank = Family
rda.pang.f.18s <- ps.pang.18s %>%
  tax_fix() %>% 
  tax_transform(rank = "Family", trans = "compositional") %>%
  ord_calc(
    constraints = c("Salinity_ppt_avg", "Temperature_C_avg", "pH_avg", "DO_mg_l_avg"),
    method = "RDA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    plot_taxa = 1:10,
    colour = "Culture_system", 
    #fill = "Culture_system",
    #shape = "Village", 
    alpha = 1,
    size = 3,
    tax_lab_style = tax_lab_style(size = 6, alpha = 0.5, fontface = "bold.italic"), # create a list of options to tweak the taxa labels' default style
    constraint_lab_style = constraint_lab_style(size = 6, alpha = 0.8) # this styles the constraint labels
  ) + 
  #scale_colour_brewer(palette = "Set1") + # can pick a different colour scale, such as a color_brewer palette
  #scale_fill_manual(values = c("Pangasius-monoculture" = "#9DCC00", "Pangasius-polyculture" = "#E7298A", "P.no.fish" = "#619CFF")) +  # Set specific fill colors for each Crop
  scale_colour_manual(values = c("Pangasius-monoculture" = "#9DCC00", "Pangasius-polyculture" = "#E7298A", "P.no.fish" = "#619CFF"),
                      labels = c("Pangasius-monoculture" = "PM", "Pangasius-polyculture" = "PP", "P.no.fish" = "Without fish")) +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = c(0.9, 0.95),
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
  labs(color = "Cultured species", # Custom legend name
       title = "RDA plot of pangasius ponds",
       subtitle = "RDA with compositional-transformed family: constraints in red, taxa in black") #+
#coord_fixed(ratio = 1, clip = "off", xlim = c(-1, 2), ylim = c(-1, 1.5)) #to set the aspect ratio of the plot
rda.pang.f.18s

#### Combine ----
pang_para.18s <- cowplot::plot_grid(p1_p.18s,
                                p2_p.18s,
                                p3_p.18s,
                                p4_p.18s,
                                ncol = 2,
                                align = "v",
                                labels = "AUTO")
pang_para.18s
# Combine spearman and RDA
comb_pang.18s <- cowplot::plot_grid(pang_para.18s,
                                rda.pang.f.18s,
                                ncol = 2,
                                align = "v",
                                labels = c("", "E"))
comb_pang.18s

# Tilapia ponds----
ps.tila.18s <- ps.18s %>% 
  ps_filter(
    Sampling_point > 6,
    Crop == "Tilapia",
    Crop_species != "Gulsha.Carp" & Crop_species != "Gulsha.Pabda")
ps.tila.18s #  3828 taxa and 255 samples

## 1. Correlation with alpha diversity----
### Estimate alpha diversity----
ps_rarefy_7t.18s <- rarefy_even_depth(ps.tila.18s, rngseed = 1234) #subsample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#281OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_7t.18s <- estimate_richness(ps_rarefy_7t.18s, measures = c("Chao1", "Shannon"))
alpha_estimates_7t.18s <- cbind(alpha_estimates_7t.18s, sample_data(ps.tila.18s))

### Plot alpha diversity----
#plotting estimated richness
#### Salinity----
p1_t.18s <- alpha_estimates_7t.18s %>% 
  ggscatter(y = "Shannon", x = "Salinity_ppt_avg",
            ylab = "Shannon", xlab = "Salinity (ppt)",
            legend = "none",
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, 
            #conf.int.level = 0.95,
            cor.method = "spearman",
            label.x.npc = "left", 
            label.y.npc = "top", # Adjust the position of the labels here
            size =2) +
  geom_point(aes(color = Salinity_ppm_avg), alpha = 1) +
  #gradient_color(c("blue", "red")) +
  scale_color_viridis_c() +  # Use viridis color scale
  labs(color = "ppm") +  # Change the legend title here
  ggtitle("B Tilapia ponds") +
  ylim(0,6) + 
  theme(plot.title = element_blank(),#element_text(size = 18, hjust = 0, face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
        axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)))  # Increase margin for y-axis title)#, face = "bold") # Adjust text size, horizontal alignment, and style
p1_t.18s

#### Temperature----
p2_t.18s <- alpha_estimates_7t.18s %>% 
  ggscatter(y = "Shannon", x = "Temperature_C_avg",
            ylab = "Shannon", xlab = "Temperature (°C)",
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE,
            legend = "none",
            cor.method = "spearman",
            label.x.npc = "left", 
            label.y.npc = "top", # Adjust the position of the labels here
            size = 2) +
  geom_point(aes(color = Temperature_C_avg)) +
  #scale_color_gradient(low = "blue", high = "red") +
  scale_color_viridis_c() +  # Use viridis color scale
  labs(color = "°C") +  # Change the legend title here
  ggtitle("B Tilapia ponds") +
  ylim(0,6) + 
  theme(plot.title = element_blank(),#element_text(size = 18, hjust = 0, face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
        axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)))  # Increase margin for y-axis title)#, face = "bold") # Adjust text size, horizontal alignment, and style
p2_t.18s

#### pH----
p3_t.18s <- alpha_estimates_7t.18s %>% 
  ggscatter(y = "Shannon", x = "pH_avg",
            add = "reg.line", conf.int = TRUE,
            ylab = "Shannon", xlab = "pH",
            cor.coef = TRUE, 
            legend = "none",
            cor.method = "spearman",
            size = 2) +
  geom_point(aes(color = pH_avg)) +
  #scale_color_gradient(low = "blue", high = "red") +
  scale_color_viridis_c() +  # Use viridis color scale
  labs(color = "pH") +  # Change the legend title here
  ggtitle("B Tilapia pond") +
  ylim(0,6) + 
  theme(plot.title = element_blank(),#element_text(size = 18, hjust = 0, face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
        axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)))  # Increase margin for y-axis title)#, face = "bold") # Adjust text size, horizontal alignment, and style
p3_t.18s

#### DO----
p4_t.18s <- alpha_estimates_7t.18s %>% 
  ggscatter(y = "Shannon", x = "DO_mg_l_avg",
            ylab = "Shannon", xlab = "DO (mg/L)",
            legend = "none",
            add = "reg.line", 
            conf.int = TRUE, 
            cor.coef = TRUE, 
            cor.method = "spearman",
            size = 2) +
  geom_point(aes(color = DO_mg_l_avg)) +
  #scale_color_gradient(low = "blue", high = "red") +
  scale_color_viridis_c() +  # Use viridis color scale
  labs(color = "mg/L") +  # Change the legend title here
  ggtitle("B Tilapia ponds") +
  ylim(0,6) + 
  theme(plot.title = element_blank(),#element_text(size = 18, hjust = 0, face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
        axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)))  # Increase margin for y-axis title)#, face = "bold") # Adjust text size, horizontal alignment, and style
p4_t.18s

## 2. RDA plot----
## Rank = Family
rda.tila.f.18s <- ps.tila.18s %>%
  tax_fix() %>%
  tax_transform(rank = "Family", trans = "compositional") %>%
  ord_calc(
    constraints = c("Salinity_ppt_avg", "Temperature_C_avg", "pH_avg", "DO_mg_l_avg"),
    method = "RDA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    plot_taxa = 1:10,
    colour = "Culture_system", 
    #fill = "Culture_system",
    #shape = "Village", 
    alpha = 1,
    size = 3,
    tax_lab_style = tax_lab_style(size = 6, alpha = 0.5, fontface = "bold.italic"), # create a list of options to tweak the taxa labels' default style
    constraint_lab_style = constraint_lab_style(size = 6, alpha = 0.8) # this styles the constraint labels
  ) + 
  #scale_colour_brewer(palette = "Set1") + # can pick a different colour scale, such as a color_brewer palette
  #scale_fill_manual(values = c("Pangasius-monoculture" = "#9DCC00", "Pangasius-polyculture" = "#E7298A", "P.no.fish" = "#619CFF")) +  # Set specific fill colors for each Crop
  scale_colour_manual(values = c("Tilapia-monoculture" = "#7cb9cb",  "Tilapia-polyculture" = "#6A3D9A"),
                      labels = c("Tilapia-monoculture" = "TM", "Tilapia-polyculture" = "TP")) +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = c(0.9, 0.95),
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
  labs(color = "Cultured species", # Custom legend name
       title = "RDA plot of tilapia ponds",
       subtitle = "RDA with compositional-transformed family: constraints in red, taxa in black") #+
#coord_fixed(ratio = 1, clip = "off"), xlim = c(-0.8, 0.8), ylim = c(-0.5, 1)) #to set the aspect ratio of the plot
rda.tila.f.18s

#### Combine ----
tila_para.18s <- cowplot::plot_grid(p1_t.18s,
                                p2_t.18s,
                                p3_t.18s,
                                p4_t.18s,
                                ncol = 2,
                                align = "v",
                                labels = c("F", "G", "H", "I"))
tila_para.18s
# Combine spearman and RDA
comb_tila.18s <- cowplot::plot_grid(tila_para.18s,
                                rda.tila.f.18s,
                                ncol = 2,
                                align = "v",
                                labels = c("", "J"))
comb_tila.18s

# Final plot----
comb_pt.18s <- cowplot::plot_grid(comb_pang.18s,
                              comb_tila.18s,
                              align = "v",
                              ncol = 1)
comb_pt.18s
# save as 2000*1600

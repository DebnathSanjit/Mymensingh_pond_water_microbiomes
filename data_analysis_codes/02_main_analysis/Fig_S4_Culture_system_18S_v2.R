# 19/09/2024 
# Sanjit Debnath
# This script is to plot the microbial composition, alpha, beta diversity and shared ASVs between pangasius and tilapia mono and polyculture system, as well as crop species (pangasius vs tilapia) 

# Load Libraries
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(microViz); packageVersion("microViz")
library(vegan); packageVersion("vegan") # needed for PERMANOVA test
library(tidyverse); packageVersion("tidyverse")
library(ggpubr); packageVersion("ggpubr")
library(cowplot); packageVersion("cowplot")
library(patchwork); packageVersion("patchwork")
library(stringr); packageVersion("stringr") # to wrap text
# Libraries for tests
library(ggpubr); packageVersion("ggpubr")
library(rstatix); packageVersion("rstatix")
#library(dunn.test); packageVersion("dunn.test")
library(VennDiagram); packageVersion("VennDiagram") # to plot venn diagram
library(grid); packageVersion("grid") # to show venn plot
library(gridExtra); packageVersion("gridExtra") # to combine venn diagram with other ggplot
library(phylosmith); packageVersion("phylosmith")

## Setup working dictionary first
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/statistical_analysis")

# Set theme for ggplot
theme_set(theme_bw())
#theme_set(theme(panel.background = element_blank(), axis.line = element_line(colour = "black")))
#when plotting with ggplot, removes grid lines etc for a clean plot N>B double check if needed
# set seed
set.seed(1234)

# Microeukaryotes (18s)----
# Load the phyloseq object
ps.18s <- readRDS("phyloseq_18S_filtered_with_tree_pr2_90-150bp_20240416.rds")
ps.18s # 5390 taxa and 872 samples

# 1. Crop species (Pangasius-Tilapia)----
crop_species <- c("Pangasius" = "#9DCC00", "Tilapia" = "#006DDBFF")

# only pangasius-tilapia
pseq_mono.18s <- ps.18s %>% 
  ps_filter(Culture_system %in% c("Pangasius-monoculture", "Tilapia-monoculture"))
pseq_mono.18s # 4844 taxa and 587 samples

ps_rarefy.mono.18s <- rarefy_even_depth(pseq_mono.18s, rngseed = 1234)
#263OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha.estimate.mono.18s <- estimate_richness(ps_rarefy.mono.18s, measures = c("Chao1", "Shannon"))
alpha.estimate.mono.18s <- cbind(alpha.estimate.mono.18s, sample_data(pseq_mono.18s))

## A. Estimate alpha diversity----
### Shannon----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
crop.s.18s <- ggplot(data = alpha.estimate.mono.18s, aes(y = Shannon, x = Crop)) +
  geom_boxplot(aes(color = Crop)) +
  geom_jitter(aes(color = Crop), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = crop_species) +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 22, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 20, margin = margin(r = 5)),
    axis.title.x = element_blank(), #text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 22, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "top",
    legend.box = "vertical",
    legend.title = element_text(size = 22, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 20),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")
  ) + guides(fill = guide_legend(ncol = 1)) + # Set number of columns in the legend
  labs (color = "Crop species",  # Customize legend title
        y = "Shannon") + # Add y-axis title to the plot
  #x = "Crop species") +
  ylim(0,6)
crop.s.18s

#### Perform Kw test----
kruskal.test(Shannon ~ Crop, data = alpha.estimate.mono.18s)
# data:  Shannon by Crop
# Kruskal-Wallis chi-squared = 13.726, df = 1, p-value = 0.0002115

##### perform Dunn test---- 
#(https://www.youtube.com/watch?v=pyLQmUfrel8)
# Dunn test is not necessary but I'm doing it to add the p value on plot. Furthermore, it doesn't change the p value
stat.test.crop.s.18s <- dunn_test(Shannon ~ Crop, data = alpha.estimate.mono.18s, p.adjust.method = "BH")

# add the value on the plotBox plot
stat.test.crop.s.18s <- stat.test.crop.s.18s %>% add_xy_position(x = "Crop")
crop.s2.18s <- crop.s.18s + stat_pvalue_manual(stat.test.crop.s.18s, label = "p.adj.signif", 
                                               step.increase = 0.05, tip.length = 0.005, 
                                               size = 8)  # Adjust the size of the asterisk here
#hide.ns = TRUE)
crop.s2.18s

## B. Beta diversity (weighted unifrac)----
#crop_wuni.18s <- pseq_mono.18s %>%
#  tax_transform(rank = "unique", trans = "compositional") %>%
#  dist_calc(dist = "wunifrac") 
#saveRDS(crop_wuni.18s, "weighted_unifrac_beta_pseq_mono.18s.rds")
crop_wuni.18s <- readRDS("weighted_unifrac_beta_pseq_mono.18s.rds")

# Culture.system
crop.we.18s <- crop_wuni.18s %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Crop", fill = "Crop",
    shape = 21, 
    alpha = 1,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(colour = Culture_system)
  ) +
  scale_colour_manual(values = crop_species) +
  scale_fill_manual(values = crop_species) +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 22, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 20, margin = margin(r = 5)),
    axis.title.x = element_text(size = 22, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 22, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "top",
    legend.box = "vertical",
    legend.title = element_text(size = 22, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 20),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  labs(colour = "Crop species",
       fill = "Crop species")
crop.we.18s

### Significance Test----
#Plot PERMANOVA with phyloseq
metadata_crop.18s <- sample_data(pseq_mono.18s) %>%
  data.frame() #%>%
#  tibble()
wuni_dist.crop.18s = phyloseq::distance(pseq_mono.18s, method = "wunifrac")
PERM_crop_wuni.18s <- adonis2(wuni_dist.crop.18s ~ Crop, data = metadata_crop.18s)

# to add test result on the plot
crop.we2.18s <- crop.we.18s + annotate(geom = "label",
                                       label = paste("PERMANOVA: R² = ", round(PERM_crop_wuni.18s["Crop","R2"], 3), 
                                                     ", p = ", PERM_crop_wuni.18s["Crop", "Pr(>F)"], sep = ""),
                                       x=Inf, y=Inf, hjust = 1, vjust = 1)
crop.we2.18s

## C. Shared taxa (Venn diagram)----
## pangasius-tilapia
# Subset tilapia
tila.18s <- pseq_mono.18s %>% 
  ps_filter(
    Crop == "Tilapia") %>% 
  tax_fix()
tila.18s # 3920 taxa and 234 samples

# Subset pangasius
pang.18s <- pseq_mono.18s %>% 
  ps_filter(
    Crop == "Pangasius") %>% 
  tax_fix()
pang.18s # 3520 taxa and 353 samples

tila.18s_asvs <- tila.18s %>% otu_table() %>%  rownames() # already did before
pang.18s_asvs <- pang.18s %>% otu_table() %>%  rownames()

venn_list_pang.tila.18s <- list("Tilapia" = tila.18s_asvs, "Pangasius" = pang.18s_asvs)

# Plot
venn_plot_pang.tila.18s <- venn.diagram(venn_list_pang.tila.18s, filename = NULL,
                                        #circles
                                        lwd = 2, lty = 'blank', fill = c("#006DDBFF", "#9DCC00"),
                                        #numbers
                                        cex = 2, # To adjust font size inside the diagram
                                        fontface = "italic", fontfamily = "sans",
                                        #set names
                                        cat.cex = 2, # To adjust font size of the variables
                                        cat.fontface = "bold",
                                        cat.default.pos = "outer",
                                        cat.fontfamily = "sans",
                                        print.mode = c("raw","percent"), # to add percentages of the shared taxa
                                        output = NULL, #to suppress the creation of text file in the working dictionary
                                        auto_scale = FALSE) # it resize the circle
#main = "d") # main title set to empty string
# Create a new plot
plot.new()

# Draw the Venn diagram
pt <- grid.arrange(venn_plot_pang.tila.18s)
pt <- grid.draw(venn_plot_pang.tila.18s)

## D. Compositional barplot----
pseq_mono2.18s <- pseq_mono.18s %>% tax_fix()
pseq_mono2.18s # 4844 taxa and 587 samples

# my palette
myPal.fam.pt.18s <- tax_palette(
  data = pseq_mono2.18s, rank = "Family", n = 25, pal = "greenArmytage", #pal = "brewerPlus",
  add = c(Other = "lightgrey")) 
myPal.fam.pt.18s["Alveolata Division"] <- "orange"
myPal.fam.pt.18s["Euglenaceae"] <- "#FFE100"  
myPal.fam.pt.18s["Raphidophyceae_XX"] <- "#5EF1F2"
myPal.fam.pt.18s["Stephanodiscaceae"] <- "#00998F"
myPal.fam.pt.18s["Gyrista Subdivision"] <- "#F0A3FF"
tax_palette_plot(myPal.fam.pt.18s)

### Family proportion-----
#### Between upazila-----
pt.fam.18s <- tax_glom(pseq_mono2.18s, taxrank = "Family")
pt.fam.18s.ra <- relative_abundance(pt.fam.18s)
pt.fam.18s.pro <- taxa_proportions(pt.fam.18s.ra, 'Family') # Overall to check who are dominant 
pt.fam.18s.pro2 <- taxa_proportions(pt.fam.18s.ra, 'Family', treatment = "Crop")

# set up for alphabetical sorting
topTaxa.crop.18s.fam <- pseq_mono2.18s %>%
  tax_top(n = 10, rank = "Family") %>%
  sort() # this makes them alphabetical

## plot with alphabetical sorting
crop.18s.fam <- pseq_mono2.18s %>%
  tax_sort(by = sum, at = "Family") %>% # this orders all genera by abundance
  ps_select(Crop_species, Crop) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Crop") %>%
  comp_barplot(tax_order = topTaxa.crop.18s.fam, # this brings the named taxa to the front
               tax_level = "Family", n_taxa = 10,
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.fam.pt.18s) +
  coord_flip() +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 22, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 24, margin = margin(r = 5)),
    axis.title.x = element_text(size = 22, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 22, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18, face = "italic")) +
  labs(title = "Pangasius vs Tilapia",
       y = "Relative abundance") +
  guides(fill = guide_legend(ncol = 1))  # Adjust the number of rows in the legend
crop.18s.fam

### Compare taxa between groups----
### Make a dataframe with dominant phyla
pseq_mono3.18s.fam <- pseq_mono2.18s %>%
  tax_transform(trans = "compositional", rank = "Family") %>%
  ps_get() %>%
  ps_otu2samdat(taxa = c("Alveolata Division","Chlamydomonadales_X","Cryptophyceae Class", 
                         "Dinophyceae Class","Euglenaceae","Eukaryota Domain",  
                         "Fungi Subdivision", "Gyrista Subdivision","Raphidophyceae_XX",
                         "Stephanodiscaceae")) %>% 
  samdat_tbl()

# run kw and dunn test for dominant taxa
kruskal.test(Stephanodiscaceae ~ Crop, data = pseq_mono3.18s.fam) # p-value  = 9.891e-08
kruskal.test(Raphidophyceae_XX ~ Crop, data = pseq_mono3.18s.fam) # p-value < 2.2e-16
kruskal.test(Euglenaceae ~ Crop, data = pseq_mono3.18s.fam) # p-value < 2.2e-16

## Culture system----
# set colors
culture_system <- c("Pangasius-monoculture" = "#9DCC00", "Pangasius-polyculture" ="#E7298A",  
                    "Tilapia-monoculture" = "#006DDBFF",  "Tilapia-polyculture" = "#490092FF")

# culture system
# to see how microbiome differ between crop species, we can take all the samples
# on those timepoints where they changed from mono to polyculture. This way we can
# avoid the effect of time. So, month 9 onward, compare pangasius and tilapia separately 
# if i mix all together, it will not give the effect of species
# 2. Pangasius samples----
pseq_pang.mp.18s <- ps.18s %>% 
  ps_filter(
    Sampling_point > 8, 
    Crop_species != "Shing" & Crop_species != "P.no.fish",
    Crop == "Pangasius",
  ) %>% 
  ps_mutate(Culture_system2 = ifelse(grepl("\\.", Crop_species), "Polyculture", "Monoculture")) %>%  #add a column in metadata named Culture_system2 where it will look for any patern with "." in the Crop_species column and if it's true, it will mark as Polyculture, if false, Monoculture
  ps_mutate(Culture_system3 = ifelse(grepl("\\.", Crop_species), "01. Polyculture", Crop_species)) #add a column in metadata named Culture.system where it will look for any patern with "." in the Crop_species column and if it's true, it will mark as Poly
pseq_pang.mp.18s # 2980 taxa and 233 samples

## A. Estimate alpha diversity (pangasius pond)----
ps_rarefy_2.3.p.18s <- rarefy_even_depth(pseq_pang.mp.18s, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#252OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_2.3.p.18s <- estimate_richness(ps_rarefy_2.3.p.18s, measures = c("Chao1", "Shannon"))
alpha_estimates_2.3.p.18s <- cbind(alpha_estimates_2.3.p.18s, sample_data(pseq_pang.mp.18s))

### Shannon----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
cs.p.s.18s <- ggplot(data = alpha_estimates_2.3.p.18s, aes(y = Shannon, x = Culture_system)) +
  geom_boxplot(aes(color = Culture_system)) +
  geom_jitter(aes(color = Culture_system), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = culture_system,
                     labels = c("Pangasius-monoculture" = "PM", 
                                "Pangasius-polyculture" = "PP" 
                                #"Tilapia-monoculture" = "TM",
                                #"Tilapia-polyculture" = "TP"
                     )) +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 22, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 20, margin = margin(r = 5)),
    axis.title.x = element_blank(), #text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 22, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "top",
    legend.box = "vertical",
    legend.title = element_text(size = 22, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 20),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  guides(fill = guide_legend(ncol = 1)) + # Set number of columns in the legend
  labs (color = "Crop species",  # Customize legend title
        y = "Chao1") +
  ylim(0, 6) +
  scale_x_discrete(labels = c("Pangasius-monoculture" = "PM", 
                              "Pangasius-polyculture" = "PP")) # Insert line break in category names
cs.p.s.18s

#### Perform Kw test----
kruskal.test(Shannon ~ Culture_system, data = alpha_estimates_2.3.p.18s)
# data:  Shannon by Culture_system
# Kruskal-Wallis chi-squared = 7.5946, df = 1, p-value = 0.005854

##### perform Dunn test---- 
#(https://www.youtube.com/watch?v=pyLQmUfrel8)
# Dunn test is not necessary but I'm doing it to add the p value on plot. Furthermore, it doesn't change the p value
stat.test.cs.p.s.18s <- dunn_test(Shannon ~ Culture_system, data = alpha_estimates_2.3.p.18s, p.adjust.method = "BH")

# add the value on the plotBox plot
stat.test.cs.p.s.18s <- stat.test.cs.p.s.18s %>% add_xy_position(x = "Culture_system")
cs.p.s2.18s <- cs.p.s.18s + stat_pvalue_manual(stat.test.cs.p.s.18s, label = "p.adj.signif", 
                                               step.increase = 0.05, tip.length = 0.005, 
                                               size = 8)  # Adjust the size of the asterisk here
#hide.ns = TRUE)
cs.p.s2.18s

## B. Beta diversity (weighted unifrac)----
#cs.p_wuni.18s <- pseq_pang.mp.18s %>%
#  tax_transform(rank = "unique", trans = "compositional") %>%
#  dist_calc(dist = "wunifrac")
#saveRDS(cs.p_wuni.18s, "weighted_unifrac_beta_diversity_pseq_pang.mp.18s.rds")
cs.p_wuni.18s <- readRDS("weighted_unifrac_beta_diversity_pseq_pang.mp.18s.rds")

# Plot
cs.p.we.18s <- cs.p_wuni.18s %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Culture_system", fill = "Culture_system",
    shape = 21, 
    alpha = 1,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(colour = Culture_system)
  ) +
  scale_colour_manual(values = culture_system,
                      labels = c("Pangasius-monoculture" = "PM", 
                                 "Pangasius-polyculture" = "PP")) +
  scale_fill_manual(values = culture_system,
                    labels = c("Pangasius-monoculture" = "PM", 
                               "Pangasius-polyculture" = "PP")) +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 22, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 20, margin = margin(r = 5)),
    axis.title.x = element_text(size = 22, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 22, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "top",
    legend.box = "vertical",
    legend.title = element_text(size = 22, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 20),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  labs(colour = "Culture system",
       fill = "Culture system")
cs.p.we.18s

### Significance Test----
#Plot PERMANOVA with phyloseq
metadata_cs.p.18s <- sample_data(pseq_pang.mp.18s) %>%
  data.frame() #%>%
#  tibble()
wuni_dist.cs.p.18s = phyloseq::distance(pseq_pang.mp.18s, method = "wunifrac")
PERM_cs.p_wuni.18s <- adonis2(wuni_dist.cs.p.18s ~ Culture_system, data = metadata_cs.p.18s)

# to add test result on the plot
cs.p.we2.18s <- cs.p.we.18s + annotate(geom = "label",
                                       label = paste("PERMANOVA: R² = ", round(PERM_cs.p_wuni.18s["Culture_system","R2"], 3), 
                                                     ", p = ", PERM_cs.p_wuni.18s["Culture_system", "Pr(>F)"], sep = ""),
                                       x=Inf, y=Inf, hjust = 1, vjust = 1)
cs.p.we2.18s

## C. Shared taxa (pangasius ponds)----
# Subset Pangasius-monoculture
pang.m.18s <- pseq_pang.mp.18s %>% 
  ps_filter(
    Culture_system == "Pangasius-monoculture") %>% 
tax_fix()
pang.m.18s # 2612 taxa and 202 samples

# Subset Pangasius-polyculture
pang.p.18s <- pseq_pang.mp.18s %>% 
  ps_filter(
    Culture_system == "Pangasius-polyculture") %>% 
tax_fix()
pang.p.18s #  1499 taxa and 31 samples

pang.m.18s_asvs <- pang.m.18s %>% otu_table() %>%  rownames() # already did before
pang.p.18s_asvs <- pang.p.18s %>% otu_table() %>%  rownames()

venn_list_pang.mp.18s <- list("PM" = pang.m.18s_asvs, "PP" = pang.p.18s_asvs)

# Plot
venn_plot_pang.mp.18s <- venn.diagram(venn_list_pang.mp.18s, filename = NULL,
                                      #circles
                                      lwd = 2, lty = 'blank', fill = c("#9DCC00", "#E7298A"),
                                      #numbers
                                      cex = 2, # To adjust font size inside the diagram
                                      fontface = "italic", fontfamily = "sans",
                                      #set names
                                      cat.cex = 2, # To adjust font size of the variables
                                      cat.fontface = "bold",
                                      cat.default.pos = "outer",
                                      cat.fontfamily = "sans",
                                      print.mode = c("raw","percent"), # to add percentages of the shared taxa
                                      output = NULL, #to suppress the creation of text file in the working dictionary
                                      auto_scale = FALSE) # it resize the circle
#main = "d") # main title set to empty string
# Create a new plot
plot.new()

# Draw the Venn diagram
pa.mp <- grid.arrange(venn_plot_pang.mp.18s)
pa.mp <- grid.draw(venn_plot_pang.mp.18s)

## D. Compositional barplot----
## subset only pangasius
pseq_pang.mp2.18s <- pseq_pang.mp.18s %>% 
 tax_fix()
pseq_pang.mp2.18s # 2980 taxa and 233 samples

# Check the relative abundance of taxa, not going to use, just to mention in which group higher or lower
pang.fam.18s <- tax_glom(pseq_pang.mp2.18s, taxrank = "Family")
pang.fam.18s.ra <- relative_abundance(pang.fam.18s)
pang.fam.18s.pro <- taxa_proportions(pang.fam.18s.ra, classification = "Family")
pang.fam.18s.pro2 <- taxa_proportions(pang.fam.18s.ra, classification = "Family", treatment = "Culture_system")

# set up for alphabetical sorting
topTaxa.p.18s.fam <- pseq_pang.mp2.18s %>%
  tax_top(n = 10, rank = "Family") %>%
  sort() # this makes them alphabetical

## plot with alphabetical sorting
p.18s.fam <- pseq_pang.mp2.18s %>%
  tax_sort(by = sum, at = "Family") %>% # this orders all genera by abundance
  ps_select(Crop_species, Culture_system) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Culture_system") %>%
  comp_barplot(tax_order = topTaxa.p.18s.fam, # this brings the named taxa to the front
               tax_level = "Family", n_taxa = 10, 
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.fam.pt.18s) +
  coord_flip() +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 22, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 24, margin = margin(r = 5)),
    axis.title.x = element_blank(), #text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 22, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18, face = "italic")) +
  labs(title = "Pangasius monoculture vs polyculture") +
  scale_x_discrete(labels = c("Pangasius-monoculture" = "PM", 
                              "Pangasius-polyculture" = "PP")) + # Insert line breaks in category names
  guides(fill = guide_legend(ncol = 1))  # Adjust the number of rows in the legend
p.18s.fam

### Compare taxa between groups----
### Make a dataframe with dominant phyla
pseq_pang.mp3.18s <- pseq_pang.mp2.18s %>%
  tax_transform(trans = "compositional", rank = "Family") %>%
  ps_get() %>%
  ps_otu2samdat(taxa = c("Alveolata Division","Dinophyceae Class","Euglenozoa Subdivision",
                        "Eukaryota Domain", "Fungi Subdivision","Gyrista Subdivision",   
                        "Kinetoplastea Class","Oxytrichidae","Stephanodiscaceae",     
                         "Thalassiosirales Order")) %>% 
  samdat_tbl()

# run kw and dunn test for dominant taxa
kruskal.test(Stephanodiscaceae ~ Culture_system, data = pseq_pang.mp3.18s) # p-value = 0.839

# 3. Tilapia samples----
pseq_tila.mp.18s <- ps.18s %>% 
  ps_filter(
    Sampling_point > 8, 
    Crop == "Tilapia",
    Crop_species != "T.no.fish",
    Crop_species != "Gulsha.Carp" & Crop_species != "Gulsha.Pabda"
  ) %>% 
  ps_mutate(Culture_system2 = ifelse(grepl("\\.", Crop_species), "Polyculture", "Monoculture")) %>%  #add a column in metadata named Culture_system2 where it will look for any patern with "." in the Crop_species column and if it's true, it will mark as Polyculture, if false, Monoculture
  ps_mutate(Culture_system3 = ifelse(grepl("\\.", Crop_species), "01. Polyculture", Crop_species)) #add a column in metadata named Culture.system where it will look for any patern with "." in the Crop_species column and if it's true, it will mark as Poly
pseq_tila.mp.18s # 3632 taxa and 215 samples 

## A. Estimate alpha diversity (tilapia ponds)----
ps_rarefy_2.3.t.18s <- rarefy_even_depth(pseq_tila.mp.18s, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#275OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_2.3.t.18s <- estimate_richness(ps_rarefy_2.3.t.18s, measures = c("Chao1", "Shannon"))
alpha_estimates_2.3.t.18s <- cbind(alpha_estimates_2.3.t.18s, sample_data(pseq_tila.mp.18s))

### Shannon----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
cs.t.s.18s <- ggplot(data = alpha_estimates_2.3.t.18s, aes(y = Shannon, x = Culture_system)) +
  geom_boxplot(aes(color = Culture_system)) +
  geom_jitter(aes(color = Culture_system), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = culture_system,
                     labels = c(#"Pangasius-monoculture" = "PM", 
                       #"Pangasius-polyculture" = "PP" 
                       "Tilapia-monoculture" = "TM",
                       "Tilapia-polyculture" = "TP"
                     )) +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 22, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 20, margin = margin(r = 5)),
    axis.title.x = element_blank(), #text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 22, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "top",
    legend.box = "vertical",
    legend.title = element_text(size = 22, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 20),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  guides(fill = guide_legend(ncol = 1)) + # Set number of columns in the legend
  labs (color = "Culture system",  # Customize legend title
        y = "Shannon", # Add y-axis title to the plot
        x = "Crop species") +
  #title = "a") +
  ylim(0,6) +
  scale_x_discrete(labels = c(#"Pangasius-monoculture" = "PM", 
    #"Pangasius-polyculture" = "PP" 
    "Tilapia-monoculture" = "TM",
    "Tilapia-polyculture" = "TP"
  )) # Insert line break in category names
cs.t.s.18s

#### Perform Kw test----
kruskal.test(Shannon ~ Culture_system, data = alpha_estimates_2.3.t.18s)
# data:  Shannon by Culture_system
# Kruskal-Wallis chi-squared = 0.77167, df = 1, p-value = 0.3797

## B. Beta diversity----
#cs.t_wuni.18s <- pseq_tila.mp.18s %>%
#  tax_transform(rank = "unique", trans = "compositional") %>%
#  dist_calc(dist = "wunifrac") 
#saveRDS(cs.t_wuni.18s, "weighted_unifrac_beta_pseq_tila_mp.18s.rds")
cs.t_wuni.18s <- readRDS("weighted_unifrac_beta_pseq_tila_mp.18s.rds")

# Culture.system
cs.t.we.18s <- cs.t_wuni.18s %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Culture_system", fill = "Culture_system",
    shape = 21, 
    alpha = 1,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(colour = Culture_system)
  ) +
  scale_colour_manual(values = culture_system,
                      labels = c(#"Pangasius-monoculture" = "PM", 
                        #"Pangasius-polyculture" = "PP" 
                        "Tilapia-monoculture" = "TM",
                        "Tilapia-polyculture" = "TP"
                      )) +
  scale_fill_manual(values = culture_system,
                    labels = c(#"Pangasius-monoculture" = "PM", 
                      #"Pangasius-polyculture" = "PP" 
                      "Tilapia-monoculture" = "TM",
                      "Tilapia-polyculture" = "TP"
                    )) +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 22, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 20, margin = margin(r = 5)),
    axis.title.x = element_text(size = 22, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 22, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "top",
    legend.box = "vertical",
    legend.title = element_text(size = 22, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 20),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  labs(colour = "Culture system",
       fill = "Culture system")
cs.t.we.18s

### Significance Test----
#Plot PERMANOVA with phyloseq
metadata_cs.t.18s <- sample_data(pseq_tila.mp.18s) %>%
  data.frame() #%>%
#  tibble()
wuni_dist.cs.t.18s = phyloseq::distance(pseq_tila.mp.18s, method = "wunifrac")
PERM_cs.t_wuni.18s <- adonis2(wuni_dist.cs.t.18s ~ Culture_system, data = metadata_cs.t.18s)

# to add test result on the plot
cs.t.we2.18s <- cs.t.we.18s + annotate(geom = "label",
                               label = paste("PERMANOVA: R² = ", round(PERM_cs.t_wuni.18s["Culture_system","R2"], 3), 
                                             ", p = ", PERM_cs.t_wuni.18s["Culture_system", "Pr(>F)"], sep = ""),
                               x=Inf, y=Inf, hjust = 1, vjust = 1)
cs.t.we2.18s

## C. Shared taxa----
# Subset Tilapia-monoculture
tila.m.18s <- pseq_tila.mp.18s %>% 
  ps_filter(
    Culture_system == "Tilapia-monoculture") %>% 
tax_fix()
tila.m.18s # 2600 taxa and 88 samples

# Subset Tilapia-polyculture
tila.p.18s <- pseq_tila.mp.18s %>% 
  ps_filter(
    Culture_system == "Tilapia-polyculture") %>% 
tax_fix()
tila.p.18s # 2795 taxa and 127 samples

tila.m.18s_asvs <- tila.m.18s %>% otu_table() %>%  rownames() # already did before
tila.p.18s_asvs <- tila.p.18s %>% otu_table() %>%  rownames()

venn_list_tila.mp.18s <- list("TM" = tila.m.18s_asvs, "TP" = tila.p.18s_asvs)

# Plot
venn_plot_tila.mp.18s <- venn.diagram(venn_list_tila.mp.18s, filename = NULL,
                                  #circles
                                  lwd = 2, lty = 'blank', fill = c("#006DDBFF", "#490092FF"),
                                  #numbers
                                  cex = 2, # To adjust font size inside the diagram
                                  fontface = "italic", fontfamily = "sans",
                                  #set names
                                  cat.cex = 2, # To adjust font size of the variables
                                  cat.fontface = "bold",
                                  cat.default.pos = "outer",
                                  cat.fontfamily = "sans",
                                  print.mode = c("raw","percent"), # to add percentages of the shared taxa
                                  output = NULL, #to suppress the creation of text file in the working dictionary
                                  auto_scale = FALSE) # it resize the circle
#main = "d") # main title set to empty string
# Create a new plot
plot.new()

# Draw the Venn diagram
ti.mp <- grid.arrange(venn_plot_tila.mp.18s)
ti.mp <- grid.draw(venn_plot_tila.mp.18s)

## D. Compositional barplot----
pseq_tila.mp2.18s <- pseq_tila.mp.18s %>% tax_fix()
pseq_tila.mp2.18s #3632 taxa and 215 samples 

# Check the relative abundace of taxa, not going to use, just to mention in which group higher or lower
tila.18s.fam <- tax_glom(pseq_tila.mp2.18s, taxrank = "Family")
tila.18s.fam.ra <- relative_abundance(tila.18s.fam)
tila.18s.fam.pro <- taxa_proportions(tila.18s.fam.ra, classification = "Family")
tila.18s.fam.pro2 <- taxa_proportions(tila.18s.fam.ra, classification = "Family", treatment = "Culture_system")

# set up for alphabetical sorting
topTaxa.t.18s.fam <- pseq_tila.mp2.18s %>%
  tax_top(n = 10, rank = "Family") %>%
  sort() # this makes them alphabetical

## plot with alphabetical sorting
t.18s.fam <- pseq_tila.mp2.18s %>%
  tax_sort(by = sum, at = "Family") %>% # this orders all genera by abundance
  ps_select(Crop_species, Culture_system) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Culture_system") %>%
  comp_barplot(tax_order = topTaxa.t.18s.fam, # this brings the named taxa to the front
               tax_level = "Family", n_taxa = 10, 
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.fam.pt.18s) +
  coord_flip() +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 22, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 20, margin = margin(r = 5)),
    axis.title.x = element_blank(), #text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 22, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18, face = "italic")) +
  labs(title = "Tilapia monoculture vs polyculture") +
  scale_x_discrete(labels = c("Tilapia-monoculture" = "TM",
                              "Tilapia-polyculture" = "TP")) +# Insert line breaks in category names
  guides(fill = guide_legend(ncol = 1))  # Adjust the number of rows in the legend
t.18s.fam

### Compare taxa between groups----
### Make a dataframe with dominant phyla
pseq_tila.mp3.18s <- pseq_tila.mp2.18s %>%
  tax_transform(trans = "compositional", rank = "Family") %>%
  ps_get() %>%
  ps_otu2samdat(taxa = c("Alveolata Division","Chlamydomonadales_X", "Cryptomonadales_X",
                         "Cryptophyceae Class","Dinophyceae Class","Eukaryota Domain",
                         "Fungi Subdivision","Gyrista Subdivision","Raphidophyceae_XX",
                         "Stephanodiscaceae")) %>% 
  samdat_tbl()

# run kw and dunn test for dominant taxa
kruskal.test(Stephanodiscaceae ~ Culture_system, data = pseq_tila.mp3.18s) # p-value = 0.4605


# Combine plots----
## Combine barplot----
cs.barplot.18s <- cowplot::plot_grid(p.18s.fam + theme(plot.title = element_blank(),axis.text.x = element_blank(), panel.border = element_blank()), 
                                     t.18s.fam + theme(plot.title = element_blank(),axis.text.x = element_blank(), panel.border = element_blank()), 
                                     crop.18s.fam + theme(plot.title = element_blank()), 
                                     nrow = 3, rel_heights = c(0.7, 0.7, 1),
                                     labels = "A",
                                     align = "v")
cs.barplot.18s # first save this with legends, edit the legend into one, then plot again without legend and otehr formating

## Combine pangasius diversity plots----
pang.cs.18s <- cowplot::plot_grid(cs.p.s2.18s,
                                  cs.p.we2.18s,
                                  align = "v",
                                  labels = c("B", "E H"),
                                  ncol = 1)
## Combina tilapia diversity plots----
tila.cs.18s <- cowplot::plot_grid(cs.t.s.18s,
                                  cs.t.we2.18s,
                                  align = "v",
                                  labels = c("C", "F I"),
                                  ncol = 1)
## Combine pangasius tilapia diversity plots----
pt.cs.18s <- cowplot::plot_grid(crop.s2.18s,
                                crop.we2.18s,
                                align = "v",
                                labels = c("D", "G J"),
                                ncol = 1)
## Combine all diversity plots----
cs.all.18s <- cowplot::plot_grid(pang.cs.18s,
                                 tila.cs.18s,
                                 pt.cs.18s, 
                                 ncol = 3) 
cs.all.18s 

## Combine venn diagram----
venn.18s <- grid.arrange(venn_plot_pang.mp.18s,
                         venn_plot_tila.mp.18s, # + theme(legend.position = "none"),
                         venn_plot_pang.tila.18s, # + theme(legend.position = "none"),
                         ncol = 3)  # Set the number of columns as per your preference

## Final plot----
fig_s4 <- grid.arrange(cs.barplot.18s,
                       cs.all.18s,
                       venn.18s,
                       heights = c(0.25, 0.5, 0.25))
# If i save this A4 size 2100*2970, only then the bars in the barplot semmes ok, otherwithe, the lower was is squeezed
# need to increase the text size also for pangasius barplot, one division is missing in legend, add that one


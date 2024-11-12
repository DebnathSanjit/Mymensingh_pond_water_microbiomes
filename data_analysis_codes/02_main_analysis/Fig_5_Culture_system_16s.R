# 02/10/2024
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

# Prokaryotes (16s)----
# Load the phyloseq object
ps <- readRDS("phyloseq_metadata_6_v3_20240419.rds")
ps #10523 taxa and 891 samples

# 1. Crop species (Pangasius-Tilapia)----
crop_species <- c("Pangasius" = "#9DCC00", "Tilapia" = "#006DDBFF")

# Only Pangasius-tilapia (all monocltured)
pseq_pt_mono <- ps %>% 
  ps_filter(Culture_system %in% c("Pangasius-monoculture", "Tilapia-monoculture"))
pseq_pt_mono # 9594 taxa and 592 samples 

## A. Estimate alpha diversity----
ps_rarefy.mono <- rarefy_even_depth(pseq_pt_mono, rngseed = 1234)
#15457TUs were removed because they are no longer 
#present in any sample after random subsampling
alpha.estimate.mono <- estimate_richness(ps_rarefy.mono, measures = c("Chao1", "Shannon"))
alpha.estimate.mono <- cbind(alpha.estimate.mono, sample_data(pseq_pt_mono))

### 1. Richness----
# Chao1: An estimator of total species richness that accounts for the presence of rare species.
crop.c <- ggplot(data = alpha.estimate.mono, aes(y = Chao1, x = Crop)) +
  geom_boxplot(aes(color = Crop)) +
  geom_jitter(aes(color = Crop), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = crop_species) +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(), #text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "top",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")
  ) + guides(fill = guide_legend(ncol = 1)) + # Set number of columns in the legend
  labs (color = "Crop species",  # Customize legend title
        y = "Chao1") +
  ylim(0,1000)
crop.c

#### Perform Kw test----
kruskal.test(Chao1 ~ Crop, data = alpha.estimate.mono)
# data:  Chao1 by Crop
# Kruskal-Wallis chi-squared = 4.1078, df = 1, p-value = 0.04269

##### perform Dunn test---- 
#(https://www.youtube.com/watch?v=pyLQmUfrel8)
# Dunn test is not necessary but I'm doing it to add the p value on plot. Furthermore, it doesn't change the p value
stat.test.crop.c <- dunn_test(Chao1 ~ Crop, data = alpha.estimate.mono, p.adjust.method = "BH")


# add the value on the plotBox plot
stat.test.crop.c <- stat.test.crop.c %>% add_xy_position(x = "Crop")
crop.c2 <- crop.c + stat_pvalue_manual(stat.test.crop.c, label = "p.adj.signif", 
                                       step.increase = 0.05, tip.length = 0.005, 
                                       size = 8)  # Adjust the size of the asterisk here
#hide.ns = TRUE)
crop.c2

### 2. Shannon----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
crop.s <- ggplot(data = alpha.estimate.mono, aes(y = Shannon, x = Crop)) +
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
  ylim(2.5,6)
crop.s

#### Perform Kw test----
kruskal.test(Shannon ~ Crop, data = alpha.estimate.mono)
# data:  Shannon by Crop
# Kruskal-Wallis chi-squared =  17.019, df = 1, p-value = 3.701e-05

##### perform Dunn test---- 
#(https://www.youtube.com/watch?v=pyLQmUfrel8)
# Dunn test is not necessary but I'm doing it to add the p value on plot. Furthermore, it doesn't change the p value
stat.test.crop.s <- dunn_test(Shannon ~ Crop, data = alpha.estimate.mono, p.adjust.method = "BH")

# add the value on the plotBox plot
stat.test.crop.s <- stat.test.crop.s %>% add_xy_position(x = "Crop")
crop.s2 <- crop.s + stat_pvalue_manual(stat.test.crop.s, label = "p.adj.signif", 
                                       step.increase = 0.05, tip.length = 0.005, 
                                       size = 8)  # Adjust the size of the asterisk here
#hide.ns = TRUE)
crop.s2

## B. Beta diversity (weighted unifrac)----
#crop_wuni <- pseq_pt_mono %>%
#  tax_transform(rank = "unique", trans = "compositional") %>%
#  dist_calc(dist = "wunifrac") 
## save this
#saveRDS(crop_wuni, "weighted_unifrac_beta_pseq_pt_mono.rds")
crop_wuni <- readRDS("weighted_unifrac_beta_pseq_pt_mono.rds")

# Culture.system
crop.we <- crop_wuni %>%
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
crop.we

### Significance Test----
#Plot PERMANOVA with phyloseq
metadata_crop <- sample_data(pseq_pt_mono) %>%
  data.frame() #%>%
#  tibble()
wuni_dist.crop = phyloseq::distance(pseq_pt_mono, method = "wunifrac")
PERM_crop_wuni <- adonis2(wuni_dist.crop ~ Crop, data = metadata_crop)

# to add test result on the plot
crop.we2 <- crop.we + annotate(geom = "label",
                               label = paste("PERMANOVA: R² = ", round(PERM_crop_wuni["Crop","R2"], 3), 
                                             ", p = ", PERM_crop_wuni["Crop", "Pr(>F)"], sep = ""),
                               x=Inf, y=Inf, hjust = 1, vjust = 1)
crop.we2

## C. Shared taxa (Venn diagram)----
## pangasius-tilapia
# for tilapia
tila <- pseq_pt_mono %>% 
  ps_filter(
    Crop == "Tilapia") %>% 
  tax_fix()
tila # 5545 taxa and 135 samples

pang <- pseq_pt_mono %>% 
  ps_filter(
    Crop == "Pangasius") %>% 
  tax_fix()
pang #  7333 taxa and 350 samples

tila_asvs <- tila %>% otu_table() %>%  rownames() # already did before
pang_asvs <- pang %>% otu_table() %>%  rownames()

venn_list_pang.tila <- list("Tilapia" = tila_asvs, "Pangasius" = pang_asvs)

venn_plot_pang.tila <- venn.diagram(venn_list_pang.tila, filename = NULL,
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
pt <- grid.arrange(venn_plot_pang.tila)

## D. Compositional barplot----
pseq_pt_mono2 <- pseq_pt_mono %>% tax_fix()
pseq_pt_mono2 # 9594 taxa and 592 samples 

# my palette
myPal.fam.pt.16s <- tax_palette(
  data = pseq_pt_mono2, rank = "Family", n = 35, pal = "brewerPlus",
  add = c(Other = "lightgrey"))
#myPal.fam.16s["Sporichthyaceae"] <- "#00D9FF" # Override existing color if any color is not good
#myPal.fam.16s["Ilumatobacteraceae"] <- "#ABE496"
#myPal.fam.16s["Microbacteriaceae"] <- "#EA70FF" 
#myPal.fam.16s["Sporichthyaceae"] <-  "#CC00A7"
myPal.fam.pt.16s["Clostridiaceae"] <- "#CC00A7"
tax_palette_plot(myPal.fam.pt.16s)

# to check the relative abundace of taxa, not going to use, just to mention in which group higher or lower
pt.fam <- tax_glom(pseq_pt_mono2, taxrank = "Family")
pt.fam.ra <- relative_abundance(pt.fam)
pt.fam.pro <- taxa_proportions(pt.fam.ra, classification = "Family")
pt.fam.pro2 <- taxa_proportions(pt.fam.ra, classification = "Family", treatment = "Crop")

# set up for alphabetical sorting
topTaxa.crop.fam <- pseq_pt_mono2 %>%
  tax_top(n = 10, rank = "Family") %>%
  sort() # this makes them alphabetical

## plot with alphabetical sorting
crop.fam <- pseq_pt_mono2 %>%
  tax_sort(by = sum, at = "Family") %>% # this orders all genera by abundance
  ps_select(Crop_species, Crop) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Crop") %>%
  comp_barplot(tax_order = topTaxa.crop.fam, # this brings the named taxa to the front
               tax_level = "Family", n_taxa = 10, 
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.fam.pt.16s) +
  coord_flip() +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 22, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 20, margin = margin(r = 5)),
    axis.title.x = element_text(size = 22, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 22, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18, face = "italic")) +
  labs(title = "Pangasius vs Tilapia",
       y = "Relative abundance") +
  guides(fill = guide_legend(ncol = 1))  # Adjust the number of rows in the legend
crop.fam

### Compare taxa between groups----
### Make a dataframe with dominant phyla
pseq_pt_mono2.fam <- pseq_pt_mono2 %>%
  tax_transform(trans = "compositional", rank = "Family") %>%
  ps_get() %>%
  ps_otu2samdat(taxa = c("Pirellulaceae", "Phycisphaeraceae", "Cyanobiaceae", # RA > 4.49%         
                         "Chthoniobacteraceae", "Clostridiaceae","Sporichthyaceae", # RA > 3.99%     
                         "Isosphaeraceae","Kapabacteriales Order","Pedosphaeraceae",        
                         "PeM15 Order")) %>% 
  samdat_tbl()

# run kw and dunn test for dominant taxa
kruskal.test(Pirellulaceae ~ Crop, data = pseq_pt_mono2.fam) # p-value = 0.001104
kruskal.test(Phycisphaeraceae ~ Crop, data = pseq_pt_mono2.fam) # p-value = 4.525e-14
kruskal.test(Cyanobiaceae ~ Crop, data = pseq_pt_mono2.fam) # p-value  = 1.508e-14
kruskal.test(Chthoniobacteraceae ~ Crop, data = pseq_pt_mono2.fam) # p-value = 0.01503
kruskal.test(Clostridiaceae ~ Crop, data = pseq_pt_mono2.fam) # p-value < 2.2e-16
kruskal.test(Sporichthyaceae ~ Crop, data = pseq_pt_mono2.fam) # p-value = 0.6188

# # Culture system----
# set colors
culture_system <- c("Pangasius-monoculture" = "#9DCC00", "Pangasius-polyculture" ="#E7298A",  
                    "Tilapia-monoculture" = "#006DDBFF",  "Tilapia-polyculture" = "#490092FF")

# to see how microbiome differ between crop species, we can take all the samples
# on those timepoints where they changed from mono to polyculture. This way we can
# avoid the effect of time. So, month 9 onward, compare pangasius and tilapia separately 
# if i mix all together, it will not give the effect of species
# 2. Pangasius samples----
pseq_pang.mp <- ps %>% 
  ps_filter(
    Sampling_point > 8, 
    Crop_species != "Shing" & Crop_species != "P.no.fish",
    Crop == "Pangasius",
   ) %>% 
  ps_mutate(Culture_system2 = ifelse(grepl("\\.", Crop_species), "Polyculture", "Monoculture")) %>%  #add a column in metadata named Culture_system2 where it will look for any patern with "." in the Crop_species column and if it's true, it will mark as Polyculture, if false, Monoculture
  ps_mutate(Culture_system3 = ifelse(grepl("\\.", Crop_species), "01. Polyculture", Crop_species)) #add a column in metadata named Culture.system where it will look for any patern with "." in the Crop_species column and if it's true, it will mark as Poly
pseq_pang.mp # 6518 taxa and 227 samples

## A. Estimate alpha diversity (pangasius pond)----
ps_rarefy_2.3.p <- rarefy_even_depth(pseq_pang.mp, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#1214OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_2.3.p <- estimate_richness(ps_rarefy_2.3.p, measures = c("Chao1", "Shannon"))
alpha_estimates_2.3.p <- cbind(alpha_estimates_2.3.p, sample_data(pseq_pang.mp))

### 1. Richness----
# Chao1: An estimator of total species richness that accounts for the presence of rare species.
cs.p.c <- ggplot(data = alpha_estimates_2.3.p, aes(y = Chao1, x = Culture_system)) +
  geom_boxplot(aes(color = Culture_system)) +
  geom_jitter(aes(color = Culture_system), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = culture_system,
                     labels = c("Pangasius-monoculture" = "PM", 
                                "Pangasius-polyculture" = "PP" 
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
    axis.line = element_line(color = "black")
  ) + guides(fill = guide_legend(ncol = 1)) + # Set number of columns in the legend
  labs (color = "Culture system",  # Customize legend title
        y = "Chao1", # Add y-axis title to the plot
        x = "Culture system")+
  ylim(0,1000)+
  #title = "a) Chao1 index of different crop species") +
  scale_x_discrete(labels = c("Pangasius-monoculture" = "PM", 
                              "Pangasius-polyculture" = "PP" 
                              #"Tilapia-monoculture" = "TM",
                              #"Tilapia-polyculture" = "TP"
                              )) # Insert line breaks in category names
cs.p.c

#### Perform Kw test----
kruskal.test(Chao1 ~ Culture_system, data = alpha_estimates_2.3.p)
# data:  Chao1 by Culture_system
# Kruskal-Wallis chi-squared = 9.2024, df = 1, p-value = 0.002417

##### perform Dunn test---- 
#(https://www.youtube.com/watch?v=pyLQmUfrel8)
# Dunn test is not necessary but I'm doing it to add the p value on plot. Furthermore, it doesn't change the p value
stat.test.cs.p.c <- dunn_test(Chao1 ~ Culture_system, data = alpha_estimates_2.3.p, p.adjust.method = "BH")

# add the value on the plotBox plot
stat.test.cs.p.c <- stat.test.cs.p.c %>% add_xy_position(x = "Culture_system")
cs.p.c2 <- cs.p.c + stat_pvalue_manual(stat.test.cs.p.c, label = "p.adj.signif", 
                                       step.increase = 0.05, tip.length = 0.005, 
                                       size = 8)  # Adjust the size of the asterisk here
#hide.ns = TRUE)
cs.p.c2

### 2. Shannon----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
cs.p.s <- ggplot(data = alpha_estimates_2.3.p, aes(y = Shannon, x = Culture_system)) +
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
  labs (color = "Culture system",  # Customize legend title
        y = "Shannon" # Add y-axis title to the plot
        #x = "Crop species"
        ) +
  #title = "a") +
  ylim(2.5, 6) +
  scale_x_discrete(labels = c("Pangasius-monoculture" = "PM", 
                              "Pangasius-polyculture" = "PP" 
                              #"Tilapia-monoculture" = "TM",
                              #"Tilapia-polyculture" = "TP"
                              )) # Insert line break in category names
cs.p.s

#### Perform Kw test----
kruskal.test(Shannon ~ Culture_system, data = alpha_estimates_2.3.p)
# data:  Shannon by Culture_system
# Kruskal-Wallis chi-squared = 18.695, df = 1, p-value = 1.534e-05

##### perform Dunn test---- 
#(https://www.youtube.com/watch?v=pyLQmUfrel8)
# Dunn test is not necessary but I'm doing it to add the p value on plot. Furthermore, it doesn't change the p value
stat.test.cs.p.s <- dunn_test(Shannon ~ Culture_system, data = alpha_estimates_2.3.p, p.adjust.method = "BH")

# add the value on the plotBox plot
stat.test.cs.p.s <- stat.test.cs.p.s %>% add_xy_position(x = "Culture_system")
cs.p.s2 <- cs.p.s + stat_pvalue_manual(stat.test.cs.p.s, label = "p.adj.signif", 
                                      step.increase = 0.05, tip.length = 0.005, 
                                      size = 8)  # Adjust the size of the asterisk here
#hide.ns = TRUE)
cs.p.s2

## B. Beta diversity (weighted unifrac)----
#Un-weighted UniFrac, dist_calc(dist = "unifrac"), does not consider the relative abundance of taxa, only their presence (detection) or absence, which can make it (overly) sensitive to rare taxa, sequencing artefacts, and abundance filtering choices. Conversely, weighted UniFrac, "wunifrac", does put (perhaps too much) more importance on highly abundant taxa, when determining dissimilarities. 

#cs.p_wuni <- pseq_pang.mp %>%
#  tax_transform(rank = "unique", trans = "compositional") %>%
#  dist_calc(dist = "wunifrac")
# save this
#saveRDS(cs.p_wuni, "weighted_unifrac_beta_diversity_pseq_pang.mp.rds")
cs.p_wuni <- readRDS("weighted_unifrac_beta_diversity_pseq_pang.mp.rds")

cs.p.we <- cs.p_wuni %>%
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
                                 "Pangasius-polyculture" = "PP" 
                                 #"Tilapia-monoculture" = "TM",
                                 #"Tilapia-polyculture" = "TP"
                                 )) +
  scale_fill_manual(values = culture_system,
                    labels = c("Pangasius-monoculture" = "PM", 
                               "Pangasius-polyculture" = "PP" 
                               #"Tilapia-monoculture" = "TM",
                               #"Tilapia-polyculture" = "TP"
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
    axis.line = element_line(color = "black")
    ) +
  labs(colour = "Culture system",
       fill = "Culture system")
cs.p.we

### Significance Test----
#Plot PERMANOVA with phyloseq
metadata_cs.p <- sample_data(pseq_pang.mp) %>%
  data.frame() #%>%
#  tibble()
wuni_dist.cs.p = phyloseq::distance(pseq_pang.mp, method = "wunifrac")
PERM_cs.p_wuni4 <- adonis2(wuni_dist.cs.p ~ Culture_system, data = metadata_cs.p)

# to add test result on the plot
cs.p.we2 <- cs.p.we + annotate(geom = "label",
                              label = paste("PERMANOVA: R² = ", round(PERM_cs.p_wuni4["Culture_system","R2"], 3), 
                                            ", p = ", PERM_cs.p_wuni4["Culture_system", "Pr(>F)"], sep = ""),
                              x=Inf, y=Inf, hjust = 1, vjust = 1)
cs.p.we2

## C. Shared taxa (pangasius ponds)----
pang.m <- pseq_pang.mp %>% 
  ps_filter(
    Culture_system == "Pangasius-monoculture") #%>% 
#tax_fix()
pang.m # 5757 taxa and 194 samples

pang.p <- pseq_pang.mp %>% 
  ps_filter(
    Culture_system == "Pangasius-polyculture") #%>% 
#tax_fix()
pang.p #  2972 taxa and 33 samples

pang.m_asvs <- pang.m %>% otu_table() %>%  rownames() # already did before
pang.p_asvs <- pang.p %>% otu_table() %>%  rownames()

venn_list_pang.mp <- list("PM" = pang.m_asvs, "PP" = pang.p_asvs)

venn_plot_pang.mp <- venn.diagram(venn_list_pang.mp, filename = NULL,
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
pa.mp <- grid.arrange(venn_plot_pang.mp)

## D. Compositional barplot----
## 1. Culture system
pseq_pang.mp #  6518 taxa and 227 samples
pseq_pang.mp2 <- pseq_pang.mp %>% tax_fix()
pseq_pang.mp2 # 6518 taxa and 227 samples

# Check the relative abundace of taxa, not going to use, just to mention in which group higher or lower
pang.fam <- tax_glom(pseq_pang.mp2, taxrank = "Family")
pang.fam.ra <- relative_abundance(pang.fam)
pang.fam.pro <- taxa_proportions(pang.fam.ra, classification = "Family")
pang.fam.pro2 <- taxa_proportions(pang.fam.ra, classification = "Family", treatment = "Culture_system")

# set up for alphabetical sorting
topTaxa.p.fam <- pseq_pang.mp2 %>%
  tax_top(n = 10, rank = "Family") %>%
  sort() # this makes them alphabetical

## plot with alphabetical sorting
p.fam <- pseq_pang.mp2 %>%
  tax_sort(by = sum, at = "Family") %>% # this orders all genera by abundance
  ps_select(Crop_species, Culture_system) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Culture_system") %>%
  comp_barplot(tax_order = topTaxa.p.fam, # this brings the named taxa to the front
               tax_level = "Family", n_taxa = 10, 
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.fam.pt.16s) +
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
  labs(title = "Pangasius monoculture vs polyculture") +
  scale_x_discrete(labels = c("Pangasius-monoculture" = "PM", 
                              "Pangasius-polyculture" = "PP")) + # Insert line breaks in category names
  guides(fill = guide_legend(ncol = 1))  # Adjust the number of rows in the legend
p.fam

### Compare taxa between groups----
### Make a dataframe with dominant family
pseq_pang.mp2.fam <- pseq_pang.mp2 %>%
  tax_transform(trans = "compositional", rank = "Family") %>%
  ps_get() %>%
  ps_otu2samdat(taxa = c("Pirellulaceae","Clostridiaceae","Phycisphaeraceae", # RA > 5%
                         "Sporichthyaceae", "Chthoniobacteraceae", "Isosphaeraceae", # 3-4%         
                         "Ilumatobacteraceae","Peptostreptococcaceae","Kapabacteriales Order","PeM15 Order"  # > 2.27%  
  )) %>% 
  samdat_tbl()

# run kw and dunn test for dominant taxa
kruskal.test(Pirellulaceae ~ Culture_system, data = pseq_pang.mp2.fam) # p-value = 4.067e-08
kruskal.test(Clostridiaceae ~ Culture_system, data = pseq_pang.mp2.fam) # p-value = 4.481e-10
kruskal.test(Phycisphaeraceae ~ Culture_system, data = pseq_pang.mp2.fam) # p-value = 0.0003132
kruskal.test(Sporichthyaceae ~ Culture_system, data = pseq_pang.mp2.fam) # p-value = 1.965e-07
kruskal.test(Chthoniobacteraceae ~ Culture_system, data = pseq_pang.mp2.fam) # p-value = 6.087e-08
kruskal.test(Isosphaeraceae ~ Culture_system, data = pseq_pang.mp2.fam) # p-value = 0.006452

# 3. Tilapia samples----
pseq_tila.mp <- ps %>% 
  ps_filter(
    Sampling_point > 8, 
    Crop == "Tilapia",
    Crop_species != "T.no.fish",
    Crop_species != "Gulsha.Carp" & Crop_species != "Gulsha.Pabda"
  ) %>% 
  ps_mutate(Culture_system2 = ifelse(grepl("\\.", Crop_species), "Polyculture", "Monoculture")) %>%  #add a column in metadata named Culture_system2 where it will look for any patern with "." in the Crop_species column and if it's true, it will mark as Polyculture, if false, Monoculture
  ps_mutate(Culture_system3 = ifelse(grepl("\\.", Crop_species), "01. Polyculture", Crop_species)) #add a column in metadata named Culture.system where it will look for any patern with "." in the Crop_species column and if it's true, it will mark as Poly
pseq_tila.mp #   7069 taxa and 227 samples 

## A. Estimate alpha diversity (tilapia ponds)----
ps_rarefy_2.3.t <- rarefy_even_depth(pseq_tila.mp, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#1716OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_2.3.t <- estimate_richness(ps_rarefy_2.3.t, measures = c("Chao1", "Shannon"))
alpha_estimates_2.3.t <- cbind(alpha_estimates_2.3.t, sample_data(pseq_tila.mp))

### 1. Richness----
# Chao1: An estimator of total species richness that accounts for the presence of rare species.
cs.t.c <- ggplot(data = alpha_estimates_2.3.t, aes(y = Chao1, x = Culture_system)) +
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
        y = "Chao1", # Add y-axis title to the plot
        x = "Crop species") +
  ylim(0,1000)+
  #title = "a) Chao1 index of different crop species") +
  scale_x_discrete(labels = c(#"Pangasius-monoculture" = "PM", 
    #"Pangasius-polyculture" = "PP" 
    "Tilapia-monoculture" = "TM",
    "Tilapia-polyculture" = "TP"
  )) # Insert line breaks in category names
cs.t.c

#### Perform Kw test----
kruskal.test(Chao1 ~ Culture_system, data = alpha_estimates_2.3.t)
# data:  Chao1 by Culture_system
# Kruskal-Wallis chi-squared = 1.085, df = 1, p-value = 0.2976

### 2. Shannon----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
cs.t.s <- ggplot(data = alpha_estimates_2.3.t, aes(y = Shannon, x = Culture_system)) +
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
  ylim(2.5,6) +
  scale_x_discrete(labels = c(#"Pangasius-monoculture" = "PM", 
    #"Pangasius-polyculture" = "PP" 
    "Tilapia-monoculture" = "TM",
    "Tilapia-polyculture" = "TP"
  )) # Insert line break in category names
cs.t.s

#### Perform Kw test----
kruskal.test(Shannon ~ Culture_system, data = alpha_estimates_2.3.t)
# data:  Shannon by Culture_system
# Kruskal-Wallis chi-squared =  0.24613, df = 1, p-value = 0.6198

## B. Beta diversity----
#cs.t_wuni <- pseq_tila.mp %>%
#  tax_transform(rank = "unique", trans = "compositional") %>%
#  dist_calc(dist = "wunifrac") 
## save this
#saveRDS(cs.t_wuni, "weighted_unifrac_beta_pseq_tila_mp.rds")
cs.t_wuni <- readRDS("weighted_unifrac_beta_pseq_tila_mp.rds")

# Culture.system
cs.t.we <- cs.t_wuni %>%
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
cs.t.we

### Significance Test----
#Plot PERMANOVA with phyloseq
metadata_cs.t <- sample_data(pseq_tila.mp) %>%
  data.frame() #%>%
#  tibble()
wuni_dist.cs.t = phyloseq::distance(pseq_tila.mp, method = "wunifrac")
PERM_cs.t_wuni4 <- adonis2(wuni_dist.cs.t ~ Culture_system, data = metadata_cs.t)

# to add test result on the plot
cs.t.we2 <- cs.t.we + annotate(geom = "label",
                                  label = paste("PERMANOVA: R² = ", round(PERM_cs.t_wuni4["Culture_system","R2"], 3), 
                                                ", p = ", PERM_cs.t_wuni4["Culture_system", "Pr(>F)"], sep = ""),
                                  x=Inf, y=Inf, hjust = 1, vjust = 1)
cs.t.we2

## C. Shared taxa----
# for tilapia
tila.m <- pseq_tila.mp %>% 
  ps_filter(
    Culture_system == "Tilapia-monoculture") %>% 
tax_fix()
tila.m # 4652 taxa and 92 samples

tila.p <- pseq_tila.mp %>% 
  ps_filter(
    Culture_system == "Tilapia-polyculture") %>% 
tax_fix()
tila.p # 5545 taxa and 135 samples

tila.m_asvs <- tila.m %>% otu_table() %>%  rownames() # already did before
tila.p_asvs <- tila.p %>% otu_table() %>%  rownames()

venn_list_tila.mp <- list("TM" = tila.m_asvs, "TP" = tila.p_asvs)

venn_plot_tila.mp <- venn.diagram(venn_list_tila.mp, filename = NULL,
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
ti.mp <- grid.arrange(venn_plot_tila.mp)

## D. Compositional barplot----
pseq_tila.mp2 <- pseq_tila.mp %>% tax_fix()
pseq_tila.mp2 # 7069 taxa and 227 samples

# to check the relative abundace of taxa, not going to use, just to mention in which group higher or lower
tila.fam <- tax_glom(pseq_tila.mp2, taxrank = "Family")
tila.fam.ra <- relative_abundance(tila.fam)
tila.fam.pro <- taxa_proportions(tila.fam.ra, classification = "Family")
tila.fam.pro2 <- taxa_proportions(tila.fam.ra, classification = "Family", treatment = "Culture_system")

# set up for alphabetical sorting
topTaxa.t.fam <- pseq_tila.mp2 %>%
  tax_top(n = 10, rank = "Family") %>%
  sort() # this makes them alphabetical

## plot with alphabetical sorting
t.fam <- pseq_tila.mp2 %>%
  tax_sort(by = sum, at = "Family") %>% # this orders all genera by abundance
  ps_select(Crop_species, Culture_system) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Culture_system") %>%
  comp_barplot(tax_order = topTaxa.t.fam, # this brings the named taxa to the front
               tax_level = "Family", n_taxa = 10, 
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.fam.pt.16s) +
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
                              "Tilapia-polyculture" = "TP")) + # Insert line breaks in category names
  guides(fill = guide_legend(ncol = 1))  # Adjust the number of rows in the legend
t.fam

### Compare taxa between groups----
### Make a dataframe with dominant phyla
pseq_tila.mp2.fam <- pseq_tila.mp2 %>%
  tax_transform(trans = "compositional", rank = "Family") %>%
  ps_get() %>%
  ps_otu2samdat(taxa = c("Pirellulaceae", "Sporichthyaceae", "Cyanobiaceae", # RA > 5.68%         
                         "Chthoniobacteraceae", "Phycisphaeraceae","Clostridiaceae",  # RA > 2.63%    
                         "Isosphaeraceae","Kapabacteriales Order","Peptostreptococcaceae",        
                         "PeM15 Order")) %>% 
  samdat_tbl()

# run kw and dunn test for dominant taxa
kruskal.test(Pirellulaceae ~ Culture_system, data = pseq_tila.mp2.fam) # p-value = 0.3511
kruskal.test(Sporichthyaceae ~ Culture_system, data = pseq_tila.mp2.fam) # p-value = 9.405e-11
kruskal.test(Cyanobiaceae ~ Culture_system, data = pseq_tila.mp2.fam) # p-value  = 0.1502
kruskal.test(Chthoniobacteraceae ~ Culture_system, data = pseq_tila.mp2.fam) # p-value = 0.03179
kruskal.test(Phycisphaeraceae ~ Culture_system, data = pseq_tila.mp2.fam) # p-value = 0.002688
kruskal.test(Clostridiaceae ~ Culture_system, data = pseq_tila.mp2.fam) # p-value = 0.01382

# Combine plots----
## combine barplot
cs.barplot <- cowplot::plot_grid(p.fam + theme( plot.title = element_blank(),axis.text.x = element_blank(), panel.border = element_blank()), 
                   t.fam + theme(plot.title = element_blank(),axis.text.x = element_blank(), panel.border = element_blank()), 
                   crop.fam + theme( plot.title = element_blank()), 
                   nrow = 3, rel_heights = c(0.7, 0.7, 1),
                   labels = "A",
                   align = "v")
cs.barplot # first save this with legends, edit the legend into one, then plot again without legend and other formating

## Combine venn diagrams
venn <- grid.arrange(venn_plot_pang.mp,
                     venn_plot_tila.mp, # + theme(legend.position = "none"),
                     venn_plot_pang.tila, # + theme(legend.position = "none"),
                     ncol = 3)  # Set the number of columns as per your preference

# Another way to plot with better allignment
# Combine pangasius diversity plots
p1 <- cowplot::plot_grid(cs.p.s2,
                          cs.p.we2,
                          align = "v",
                         labels = c("B", "E"),
                          ncol = 1)
# Combine tilapia diversity plots
p2 <- cowplot::plot_grid(cs.t.s,
                        cs.t.we2,
                        align = "v",
                        labels = c("C", "F"),
                        ncol = 1)
# Combine pangasius tilapia diversity plots
p3 <- cowplot::plot_grid(crop.s2,
                        crop.we2,
                        align = "v",
                        labels = c("D", "G"),
                        ncol = 1)
# Combine all three diversity plots
p <- cowplot::plot_grid(p1,p2,p3, ncol = 3) 
p 

# Final plot
grid.arrange(cs.barplot,
             p,
             venn,
             heights = c(0.25, 0.5, 0.25))
# Save as 2100*2970 (A4 size, so that each plot keeps the ratio properly)


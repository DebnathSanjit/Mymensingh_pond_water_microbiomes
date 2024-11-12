# Date: 20240612
# Author: Sanjit Debnath
# This code is to see how disease and algal bloom effem the microbial composition, diveristy

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
library(gridExtra); packageVersion("gridExtra")
# Libraries for tests
library(ggpubr); packageVersion("ggpubr")
library(rstatix); packageVersion("rstatix")
library(dunn.test); packageVersion("dunn.test")

## Setup working dictionary first
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/statistical_analysis")

# Set theme for ggplot
theme_set(theme_bw())
set.seed(1234)

##set variables desired color: all variables
dis.colors <- c("Non-diseased"= "#4363d8", "Diseased" = "#cc79a7", "Algal.bloom"= "#5EC747")

# Load the phyloseq object----
ps <- readRDS("phyloseq_metadata_6_v3_20240419.rds")
ps #10523 taxa and 891 samples

# Prokaryotes----
pseq_pang <- ps %>%
  ps_filter(
    Crop == "Pangasius",
    Crop_species != "Shing") #%>% 
#tax_fix() #%>% 
pseq_pang # 8091 taxa and 421 samples

## A. Estimate Alpha diversity----
set.seed(1234)#The set. seed() function in R is used to create reproducible results when writing code that involves creating variables that take on random values. By using the set. seed() function, you guarantee that the same random values are produced each time you run the code
ps_rarefy_dis <- rarefy_even_depth(pseq_pang, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#1346OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_dis <- estimate_richness(ps_rarefy_dis, measures = c("Chao1", "Shannon"))
alpha_estimates_dis <- cbind(alpha_estimates_dis, sample_data(pseq_pang))

# Reorder variables
alpha_estimates_dis$Reported_disease <- factor(alpha_estimates_dis$Reported_disease, levels = c("Non-diseased", "Diseased", "Algal.bloom"))

### Shannon----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
p.dis.s <- ggplot(data = alpha_estimates_dis, aes(y = Shannon, x = Reported_disease)) +
  geom_boxplot(aes(color = Reported_disease), fill = NA) +  # Outline color, no fill for box
  geom_jitter(aes(color = Reported_disease), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  scale_color_manual(values = dis.colors, labels = c("ND", "DS", "BL")) +  # Ensure points and box outlines use the same colors
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
    panel.grid = element_blank()) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_x_discrete(labels = c("Non-diseased" = "ND", 
                              "Diseased" = "DS", 
                              "Algal.bloom" = "BL"))  + # Insert line breaks in category names
  labs(color = "Reported disease",  # Customize legend title
       y = "Shannon",  # Add y-axis title to the plot
       x = "Reported disease",
       title = "Alpha diversity") +
  guides(color = guide_legend(nrow = 1))  # Set number of columns in the legend
p.dis.s

#### Perform Kw test----
kruskal.test(Shannon ~ Reported_disease, data = alpha_estimates_dis)
#Kruskal-Wallis chi-squared = 12.892, df = 2, p-value = 0.001587

# perform Dunn test (https://www.youtube.com/watch?v=pyLQmUfrel8)
stat.test.d.s <- dunn_test(Shannon ~ Reported_disease, data = alpha_estimates_dis, p.adjust.method = "BH")

# add the value on the plotBox plot
stat.test.d.s2 <- stat.test.d.s %>% add_xy_position(x = "Reported_disease")
p.dis.s2 <- p.dis.s + stat_pvalue_manual(stat.test.d.s2, label = "p.adj.signif", 
                                         step.increase = 0.05, tip.length = 0.005, 
                                         size = 8)#,  # Adjust the size of the asterisk here
#hide.ns = TRUE)
p.dis.s2

## B. Beta diversity (weighted unifrac)----
# wunifrac:weighted UniFrac took about #01:59 - 02:12 = 13 minutes to calculate the distance
#dis_wuni <- pseq_pang %>%
#  tax_transform(rank = "unique", trans = "compositional") %>%
#  dist_calc(dist = "wunifrac") 
dis_wuni <- readRDS("weighted_unifrac_beta_diversity_pseq_pang.rds")

# Plot
p.dis.we <- dis_wuni %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Reported_disease", 
    #fill = "Reported_disease",
    shape = 16, 
    alpha = 1,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(colour = Reported_disease)
  ) +
  scale_color_manual(values = dis.colors, labels = c("ND", "DS", "BL")) +  # Ensure points and box outlines use the same colors
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
    #panel.border = element_blank(),
    #panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  scale_x_discrete(labels = c("Non-diseased" = "ND", 
                              "Diseased" = "DS", 
                              "Algal.bloom" = "BL"))  + # Insert line breaks in category names
  labs(color = "Reported disease",  # Customize legend title
       #y = "Chao1",  # Add y-axis title to the plot
       #x = "Reported disease",
       title = "PCoA") +
  guides(color = guide_legend(nrow = 1)) + # Set number of columns in the legend
  xlim(-1.5, 3.5) +
  ylim(-1.5, 1.75)
p.dis.we

### PERMANOVA Test----
#Plot PERMANOVA with phyloseq
metadata_pang.dis <- sample_data(pseq_pang) %>%
  data.frame() %>%
  tibble()
wuni_dist.dis = phyloseq::distance(pseq_pang, method="wunifrac")
PERM_dis_wuni <- adonis2(wuni_dist.dis ~ Reported_disease, data = metadata_pang.dis)

# to add test result on the plot
p.dis.we2 <- p.dis.we + annotate(geom = "label",
                                 label = paste("PERMANOVA: R² = ", round(PERM_dis_wuni["Reported_disease","R2"], 3), 
                                               ", p = ", PERM_dis_wuni["Reported_disease", "Pr(>F)"], sep = ""),
                                 x=Inf, y=Inf, hjust = 1, vjust = 1)
p.dis.we2

#### Pairwise adonis----
pairwise_adonis_dis <- pairwise.adonis(wuni_dist.dis, metadata_pang.dis$Reported_disease, p.adjust.m = "BH")
pairwise_adonis_dis

# Add both PERMANOVA and pairwise adonis results on the plot
p.dis.we3 <- p.dis.we + 
  annotate(
    geom = "label",
    label = paste(
      "PERMANOVA: R² = ", round(PERM_dis_wuni["Reported_disease","R2"], 3), 
      ", p =", PERM_dis_wuni["Reported_disease", "Pr(>F)"], "\n",
      "DS vs ND: R² =", round(pairwise_adonis_dis[1, "R2"], 3),
      ", FDR =", pairwise_adonis_dis[1, "p.adjusted"], "\n",
      "DS vs BL: R² =", round(pairwise_adonis_dis[2, "R2"], 3),
      ", FDR =", pairwise_adonis_dis[2, "p.adjusted"], "\n",
      "ND vs BL: R² =", round(pairwise_adonis_dis[3, "R2"], 3),
      ", FDR =", pairwise_adonis_dis[3, "p.adjusted"]
    ),
    x = Inf, y = Inf, hjust = 1.05, vjust = 0.8)
p.dis.we3

## C. Barplots----
# phyloseq with tax_fix
pseq_pang2 <- pseq_pang %>% tax_fix() 

# Reorder disease names
dis_order <- c("Non-diseased", "Diseased", "Algal.bloom")
sample_data(pseq_pang2)$Reported_disease <- factor(sample_data(pseq_pang2)$Reported_disease, levels = dis_order)

# Rank = Family
# my palette
myPal.fam.16s <- tax_palette(
  data = pseq_pang2, rank = "Family", n = 20, pal = "brewerPlus",
  add = c(Other = "lightgrey"))
myPal.fam.16s["Chthoniobacteraceae"] <- "#00D9FF" # Override existing color if any color is not good
#myPal.fam.16s["Cyanobiaceae"] <- "#ABE496"
myPal.fam.16s["Clostridiaceae"] <- "#EA70FF"
#myPal.fam.16s["PeM15 Order"] <- "#CC00A7"
tax_palette_plot(myPal.fam.16s)

# set up for alphabetical sorting
topTaxa.pang.fam <- pseq_pang2 %>%
  tax_top(n = 20, rank = "Family") %>%
  sort() # this makes them alphabetical

## plot with alphabetical sorting
dis.fam <- pseq_pang2 %>%
  tax_sort(by = sum, at = "Family") %>% # this orders all genera by abundance
  ps_select(Crop_species, Reported_disease) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Reported_disease") %>%
  comp_barplot(tax_order = topTaxa.pang.fam, # this brings the named taxa to the front
               tax_level = "Family", n_taxa = 20, # RA > 2%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.fam.16s) +
  #coord_flip() +
  theme(plot.title = element_text(hjust = 0, size = 18, margin = margin(b = 15)), # Increase margin for plot title
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
        legend.position = "right",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18, face = "italic")) +
  labs(title = "Prokaryotes: Pangasius ponds",
       y = "Relative abundance") +
  guides(fill = guide_legend(ncol = 1)) + # Adjust the number of rows in the legend
  scale_x_discrete(labels = c("Non-diseased" = "ND", 
                              "Diseased" = "DS", 
                              "Algal.bloom" = "BL"))
dis.fam

### Compare taxa between groups----
#### Family proportion-----
fam.abd.pang <- taxa_proportions(pseq_pang2, 'Family', treatment = "Reported_disease")
# Save the data frame as an Excel file

### Compare taxa between groups----
### Make a dataframe with dominant phyla
pseq_pang3 <- pseq_pang2 %>%
  tax_transform(trans = "compositional", rank = "Family") %>%
  ps_get() %>%
  ps_otu2samdat(taxa = c("Bacteria Kingdom","Beijerinckiaceae","Burkholderiaceae",         
                         "Chthoniobacteraceae","Clostridiaceae", "Comamonadaceae",           
                         "Cyanobacteriia Class", "Cyanobiaceae",  "Gammaproteobacteria Class",
                         "Ilumatobacteraceae", "Isosphaeraceae","Kapabacteriales Order",    
                         "MWH-UniP1 aquatic group","Pedosphaeraceae", "PeM15 Order",              
                         "Peptostreptococcaceae", "Phycisphaeraceae", "Pirellulaceae",            
                         "Rubinisphaeraceae", "Sporichthyaceae")) %>% 
  samdat_tbl()

# run kw and dunn test for dominant taxa
kruskal.test(`Cyanobacteriia Class` ~ Reported_disease, data = pseq_pang3) # p-value = 0.002289
dunn_test(`Cyanobacteriia Class` ~ Reported_disease, p.adjust.method = "BH",data = pseq_pang3)

kruskal.test(Comamonadaceae ~ Reported_disease, data = pseq_pang3) # p-value = 0.02095
dunn_test(Comamonadaceae ~ Reported_disease, p.adjust.method = "BH",data = pseq_pang3)

kruskal.test(Pirellulaceae ~ Reported_disease, data = pseq_pang3) # p-value =  0.002896
dunn_test(Pirellulaceae ~ Reported_disease, p.adjust.method = "BH",data = pseq_pang3)

kruskal.test(Phycisphaeraceae ~ Reported_disease, data = pseq_pang3) # p-value = 0.0003602
dunn_test(Phycisphaeraceae ~ Reported_disease, p.adjust.method = "BH",data = pseq_pang3)

kruskal.test(Cyanobiaceae ~ Reported_disease, data = pseq_pang3) # p-value = 0.01165
dunn_test(Cyanobiaceae ~ Reported_disease, p.adjust.method = "BH",data = pseq_pang3)

kruskal.test(Phycisphaeraceae ~ Reported_disease, data = pseq_pang3) # p-value = 0.0003602
dunn_test(Phycisphaeraceae ~ Reported_disease, p.adjust.method = "BH",data = pseq_pang3)

kruskal.test(Sporichthyaceae ~ Reported_disease, data = pseq_pang3) #p-value =  0.003389
dunn_test(Sporichthyaceae ~ Reported_disease, p.adjust.method = "BH",data = pseq_pang3)
# All are significantly different
# except Sporichthyaceae, others were > 5%

# Microeukaryotes----
# Load the phyloseq object
ps.18s <- readRDS("phyloseq_18S_filtered_with_tree_pr2_90-150bp_20240416.rds")
ps.18s # 5390 taxa and 872 samples

# Pangaisu ponds
pseq_pang.18s <- ps.18s %>%
  ps_filter(
    Crop == "Pangasius",
    Crop_species != "Shing") #%>% 
#tax_fix() #%>% 
pseq_pang.18s # 3947 taxa and 421 samples

## A. Estimate Alpha diversity----
set.seed(1234)#The set. seed() function in R is used to create reproducible results when writing code that involves creating variables that take on random values. By using the set. seed() function, you guarantee that the same random values are produced each time you run the code
ps_rarefy_dis.18s <- rarefy_even_depth(pseq_pang.18s, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#285OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_dis.18s <- estimate_richness(ps_rarefy_dis.18s, measures = c("Chao1", "Shannon"))
alpha_estimates_dis.18s <- cbind(alpha_estimates_dis.18s, sample_data(pseq_pang.18s))

# Reorder variables
alpha_estimates_dis.18s$Reported_disease <- factor(alpha_estimates_dis.18s$Reported_disease, levels = c("Non-diseased", "Diseased", "Algal.bloom"))

### Shannon----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
p.dis.s.18s <- ggplot(data = alpha_estimates_dis.18s, aes(y = Shannon, x = Reported_disease)) +
  geom_boxplot(aes(color = Reported_disease), fill = NA) +  # Outline color, no fill for box
  geom_jitter(aes(color = Reported_disease), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  scale_color_manual(values = dis.colors, labels = c("ND", "DS", "BL")) +  # Ensure points and box outlines use the same colors
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
    panel.grid = element_blank()) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_x_discrete(labels = c("Non-diseased" = "ND", 
                              "Diseased" = "DS", 
                              "Algal.bloom" = "BL"))  + # Insert line breaks in category names
  labs(color = "Reported disease",  # Customize legend title
       y = "Shannon",  # Add y-axis title to the plot
       x = "Reported disease",
       title = "Alpha diversity") +
  guides(color = guide_legend(nrow = 1))  # Set number of columns in the legend
p.dis.s.18s

#### Perform Kw test----
kruskal.test(Shannon ~ Reported_disease, data = alpha_estimates_dis.18s)
#Kruskal-Wallis chi-squared = 9.2783, df = 2, p-value = 0.009666

# perform Dunn test (https://www.youtube.com/watch?v=pyLQmUfrel8)
stat.test.d.s.18s <- dunn_test(Shannon ~ Reported_disease, data = alpha_estimates_dis.18s, p.adjust.method = "BH")

# add the value on the plotBox plot
stat.test.d.s2.18s <- stat.test.d.s.18s %>% add_xy_position(x = "Reported_disease")
p.dis.s2.18s <- p.dis.s.18s + stat_pvalue_manual(stat.test.d.s2.18s, label = "p.adj.signif", 
                                                 step.increase = 0.05, tip.length = 0.005, 
                                                 size = 8)  # Adjust the size of the asterisk here
#hide.ns = TRUE)
p.dis.s2.18s

## B. Beta diversity (weighted unifrac)----
# wunifrac:weighted UniFrac took about #01:59 - 02:12 = 13 minutes to calculate the distance
#dis_wuni.18s <- pseq_pang.18s %>%
#  tax_transform(rank = "unique", trans = "compositional") %>%
#  dist_calc(dist = "wunifrac") 
dis_wuni.18s <- readRDS("weighted_unifrac_beta_diversity_pseq_pang_18s.rds")

# Reported_disease
p.dis.we.18s <- dis_wuni.18s %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Reported_disease", 
    #fill = "Reported_disease",
    shape = 16, 
    alpha = 1,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(colour = Reported_disease)
  ) +
  scale_color_manual(values = dis.colors, labels = c("ND", "DS", "BL")) +  # Ensure points and box outlines use the same colors
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
    #panel.border = element_blank(),
    #panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  scale_x_discrete(labels = c("Non-diseased" = "ND", 
                              "Diseased" = "DS", 
                              "Algal.bloom" = "BL"))  + # Insert line breaks in category names
  labs(color = "Reported disease",  # Customize legend title
       #y = "Chao1",  # Add y-axis title to the plot
       #x = "Reported disease",
       title = "PCoA") +
  guides(color = guide_legend(nrow = 1)) + # Set number of columns in the legend
  xlim(-2, 1.5) + ylim(-1.5, 1.5)
p.dis.we.18s

### PERMANOVA Test----
#Plot PERMANOVA with phyloseq
metadata_pang.dis.18s <- sample_data(pseq_pang.18s) %>%
  data.frame() %>%
  tibble()
wuni_dist.dis.18s = phyloseq::distance(pseq_pang.18s, method="wunifrac")
PERM_dis_wuni.18s <- adonis2(wuni_dist.dis.18s ~ Reported_disease, data = metadata_pang.dis.18s)

# to add test result on the plot
p.dis.we2.18s <- p.dis.we.18s + annotate(geom = "label",
                                         label = paste("PERMANOVA: R² = ", round(PERM_dis_wuni.18s["Reported_disease","R2"], 3), 
                                                       ", p = ", PERM_dis_wuni.18s["Reported_disease", "Pr(>F)"], sep = ""),
                                         x=Inf, y=Inf, hjust = 1, vjust = 1)
p.dis.we2.18s

#### Pairwise adonis----
pairwise_adonis_dis.18s <- pairwise.adonis(wuni_dist.dis.18s, metadata_pang.dis.18s$Reported_disease, p.adjust.m = "BH")
pairwise_adonis_dis.18s

# Add both PERMANOVA and pairwise adonis results on the plot
p.dis.we3.18s <- p.dis.we.18s + 
  annotate(
    geom = "label",
    label = paste(
      "PERMANOVA: R² = ", round(PERM_dis_wuni.18s["Reported_disease","R2"], 3), 
      ", p =", PERM_dis_wuni.18s["Reported_disease", "Pr(>F)"], "\n",
      "DS vs ND: R² =", round(pairwise_adonis_dis.18s[1, "R2"], 3),
      ", FDR =", pairwise_adonis_dis.18s[1, "p.adjusted"], "\n",
      "DS vs BL: R² =", round(pairwise_adonis_dis.18s[2, "R2"], 3),
      ", FDR =", pairwise_adonis_dis.18s[2, "p.adjusted"], "\n",
      "ND vs BL: R² =", round(pairwise_adonis_dis.18s[3, "R2"], 3),
      ", FDR =", pairwise_adonis_dis.18s[3, "p.adjusted"]
    ),
    x = Inf, y = Inf, hjust = 1.05, vjust = 0.8)
p.dis.we3.18s

## C. Barplot----
## subset only pangasius (18s)
pseq_pang.18s2 <- ps.18s %>% tax_fix()
pseq_pang.18s2 # 3947 taxa and 421 samples

# Reorder disease names
dis_order <- c("Non-diseased", "Diseased", "Algal.bloom")
sample_data(pseq_pang.18s2)$Reported_disease <- factor(sample_data(pseq_pang.18s2)$Reported_disease, levels = dis_order)

# Rank = Family
# my palette
myPal.fam.18s <- tax_palette(
  data = pseq_pang.18s2, rank = "Family", n = 35, pal = "brewerPlus",
  add = c(Other = "lightgrey"))
myPal.fam.18s["Alveolata Division"] <- "#CC00A7"
myPal.fam.18s["Chlamydomonadales_X"] <- "#EA70FF"
#myPal.fam.18s[" Oxytrichidae"] <- "#EA70FF"
tax_palette_plot(myPal.fam.18s)

# set up for alphabetical sorting
topTaxa.pang.fam.18s <- pseq_pang.18s2 %>%
  tax_top(n = 20, rank = "Family") %>%
  sort() # this makes them alphabetical

## plot with alphabetical sorting
dis.fam.18s <- pseq_pang.18s2 %>%
  tax_sort(by = sum, at = "Family") %>% # this orders all genera by abundance
  ps_select(Crop_species, Reported_disease) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Reported_disease") %>%
  comp_barplot(tax_order = topTaxa.pang.fam.18s, # this brings the named taxa to the front
               tax_level = "Family", n_taxa = 20, 
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.fam.18s) +
  #coord_flip() +
  theme(plot.title = element_text(hjust = 0, size = 18, margin = margin(b = 15)), # Increase margin for plot title
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
        legend.position = "right",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18, face = "italic")) +
  scale_x_discrete(labels = c("Non-diseased" = "ND", 
                              "Diseased" = "DS", 
                              "Algal.bloom" = "BL"))  + # Insert line breaks in category names
  labs(title = "Microeukaryotic composition at Family level",
       y = "Relative abundance") +
  guides(fill = guide_legend(ncol = 1))  # Adjust the number of rows in the legend
dis.fam.18s

#### Family proportion-----
fam.abd.pang.18s <- taxa_proportions(pseq_pang.18s2, 'Family', treatment = "Reported_disease")
# Save the data frame as an Excel file

### Compare taxa between groups----
### Make a dataframe with dominant phyla
pseq_pang.18s3 <- pseq_pang.18s2 %>%
  tax_transform(trans = "compositional", rank = "Family") %>%
  ps_get() %>%
  ps_otu2samdat(taxa = c("Alveolata Division","Aulacoseiraceae","Chlamydomonadales_X",   
                         "Cryptomonadales_X", "Cryptophyceae Class","Dinophyceae Class",     
                         "Euglenaceae", "Euglenozoa Subdivision", "Eukaryota Domain",      
                         "Fungi Subdivision","Gyrista Subdivision","Kinetoplastea Class",   
                         "Oxytrichidae", "Peridiniales Order", "Raphidophyceae_XX",     
                         "Spirotrichea Class", "Stephanodiscaceae","Thalassiosirales Order",
                         "Tintinnidiidae", "TSAR Supergroup")) %>% 
  samdat_tbl()

# run kw and dunn test for dominant taxa
kruskal.test(Euglenaceae ~ Reported_disease, data = pseq_pang.18s3) # p-value = 1.69e-07
dunn_test(Euglenaceae ~ Reported_disease, p.adjust.method = "BH",data = pseq_pang.18s3)

kruskal.test(`Chlamydomonadales_X` ~ Reported_disease, data = pseq_pang.18s3) # p-value = 0.002026
dunn_test(`Chlamydomonadales_X` ~ Reported_disease, p.adjust.method = "BH",data = pseq_pang.18s3)

kruskal.test(Stephanodiscaceae ~ Reported_disease, data = pseq_pang.18s3) # p-value = 0.001009
dunn_test(Stephanodiscaceae ~ Reported_disease, p.adjust.method = "BH",data = pseq_pang.18s3)

kruskal.test(`Nephroselmidales_X` ~ Reported_disease, data = pseq_pang.18s3) # p-value = 0.0003602
dunn_test(`Nephroselmidales_X` ~ Reported_disease, p.adjust.method = "BH",data = pseq_pang.18s3)
# this is not included in the top 20 or even top 40 families, because in diseased and non-diseased, it is in very low quanity while in bloom it's above 10%. For now, just ingone this

kruskal.test(Aulacoseiraceae ~ Reported_disease, data = pseq_pang.18s3) # p-value = 0.009936
dunn_test(Aulacoseiraceae ~ Reported_disease, p.adjust.method = "BH",data = pseq_pang.18s3)
# All are significantly different

# Combine plots----
# Combine barplots
p1 <- cowplot::plot_grid(dis.fam,
                         dis.fam.18s,
                         labels = c("A", "D"),
                         align = "v",
                         ncol = 1)
# Combine diversity plots
p2 <- cowplot::plot_grid(p.dis.s2 + theme(legend.position = "none", axis.title.x = element_blank()),
                         p.dis.we3 + theme(legend.position = "none"),
                         p.dis.s2.18s + theme(legend.position = "none", axis.title.x = element_blank()),
                         p.dis.we3.18s + theme(legend.position = "none"),
                         labels = c("B","C", "E", "F"),
                         align = "h",
                         ncol = 2)
# Combine both
p <- cowplot::plot_grid(p1, p2, ncol = 2, 
                        rel_widths = c(4,6),
                        align = "h")
# Get legend
legend.dis <- get_legend(p.dis.s2)

# Final plot 
fig_7 <- cowplot::plot_grid(p,
                              legend.dis,
                              rel_heights = c(2, 0.2),
                              ncol = 1)
fig_7
# Save as 2000*1200

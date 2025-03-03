# 17/09/2024 
# Sanjit Debnath
# This script is to plot alpha and beta diversity of prokaryotes and microeukaryotes of pangasius and tilapia pond samples separately, to see the effect of different geographical locations

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
library(dunn.test); packageVersion("dunn.test")
library(VennDiagram); packageVersion("VennDiagram") # to plot venn diagram
# to combine venn diagram with other ggplot
library(gridExtra); packageVersion("gridExtra")

## Setup working dictionary first
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/statistical_analysis")

# Set theme for ggplot
theme_set(theme_bw())
set.seed(1234)

# Load the phyloseq object
ps <- readRDS("phyloseq_metadata_6_v3_20240419.rds")
ps #10523 taxa and 891 samples

# Define colors for each upazila
upazila.colors <- c("Jamalpur_Sadar" = "#1B9E77FF", "Muktagacha" = "#E6AB02FF", "Tarakanda" = "#490092FF")
upazila.colors2 <- c("JS" = "#1B9E77FF", "MG" = "#E6AB02FF", "TK" = "#490092FF")

# 1. Prokaryotes----
# Pangasius only----
pseq_pang <- ps %>%
  ps_filter(
    #Reported_disease == "Non-diseased",
    Crop == "Pangasius" & Crop_species != "Shing")
#Crop_species != "P.no.fish") #%>% 
#tax_fix()# %>% 
pseq_pang # 8091 taxa and 421 samples

## A1. Alpha diversity----
set.seed(1234)#The set. seed() function in R is used to create reproducible results when writing code that involves creating variables that take on random values. By using the set. seed() function, you guarantee that the same random values are produced each time you run the code
ps_rarefy_2.1 <- rarefy_even_depth(pseq_pang, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#1346OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_2.1 <- estimate_richness(ps_rarefy_2.1, measures = c("Chao1", "Shannon"))
alpha_estimates_2.1 <- cbind(alpha_estimates_2.1, sample_data(pseq_pang))

## upazila
### Shannon----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
p.up.s <- ggplot(data = alpha_estimates_2.1, aes(y = Shannon, x = Upazila)) +
  geom_boxplot(aes(color = Upazila)) +
  geom_jitter(aes(color = Upazila), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = upazila.colors,
                     labels = c("Jamalpur_Sadar" = "Jamalpur Sadar")) +
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
    legend.text = element_text(size = 18),
    panel.grid = element_blank()) +
  guides(fill = guide_legend(nrow = 1)) + # Set number of columns in the legend
  labs (y = "Shannon", # Add y-axis title to the plot
        x = "Upazila") +
  #title = "a") +
  scale_x_discrete(labels = c("Jamalpur_Sadar" = "Jamalpur Sadar")) + # Insert line break in category names
  ylim(2.5, 6)
p.up.s

#### Perform Kw test----
kruskal.test(Shannon ~ Upazila, data = alpha_estimates_2.1)
#data:  Shannon by Upazila
#Kruskal-Wallis chi-squared = 23.914, df = 1, p-value = 1.008e-06

##### Perform Dunn test (https://www.youtube.com/watch?v=pyLQmUfrel8)----
stat.test.pvs <- dunn_test(Shannon ~ Upazila, data = alpha_estimates_2.1, p.adjust.method = "BH")
# actually for two variables i don't need dunn test but to automatically add the value on 
# plot, i did it. It doesn't change the value

# add the value on the plotBox plot
stat.test.pvs <- stat.test.pvs %>% add_xy_position(x = "Upazila")
p.up.s2 <- p.up.s + stat_pvalue_manual(stat.test.pvs, label = "p.adj.signif", 
                                       step.increase = 0.05, tip.length = 0.005, 
                                       size = 8,  # Adjust the size of the asterisk here
                                       hide.ns = TRUE)
p.up.s2

## B1. Beta diversity (weighted unifrac)----
# wunifrac:weighted UniFrac took about #03:10 - 03:21 = 11 minutes to calculate the distance
#up_wuni <- pseq_pang %>%
#  tax_transform(rank = "unique", trans = "compositional") %>%
#  dist_calc(dist = "wunifrac") 
#saveRDS(up_wuni, "weighted_unifrac_beta_diversity_for_upazila_pangasius.rds")
up_wuni <- readRDS("weighted_unifrac_beta_diversity_for_upazila_pangasius.rds")

# Upazila_pangasius pond
p.up.we <- up_wuni %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Upazila", 
    fill = "Upazila",
    shape = 16, 
    alpha = 1,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(colour = Upazila)
  ) +
  scale_color_manual(values = upazila.colors,
                     labels = c("Jamalpur_Sadar" = "Jamalpur Sadar")) +
  scale_fill_manual(values = upazila.colors,
                    labels = c("Jamalpur_Sadar" = "Jamalpur Sadar")) +
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
    legend.text = element_text(size = 18),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  guides(fill = guide_legend(nrow = 1)) #+ # Set number of columns in the legend
#labs( #x = "PCoA1 [25.7%]",
#y = "PCoA2 [16.2%]",
#  title = "b")
p.up.we

### PERMANOVA----
metadata_pang <- sample_data(pseq_pang) %>%
  data.frame() %>%
  tibble()
wuni_dist = phyloseq::distance(pseq_pang, method="wunifrac")
PERM_up_wuni <- adonis2(wuni_dist ~ Upazila, data = metadata_pang)

# to add test result on the plot
p.up.we2 <- p.up.we + annotate(geom = "label",
                               label = paste("PERMANOVA: R² = ", round(PERM_up_wuni["Upazila","R2"], 3), 
                                             ", p = ", PERM_up_wuni["Upazila", "Pr(>F)"], sep = ""),
                               x=Inf, y=Inf, hjust = 1, vjust = 1)
p.up.we2


# Combine pangasius plots----
pang.up <- cowplot::plot_grid(p.up.s2 + theme(legend.position = "none"), 
                              #axis.text.x = element_blank()),
                              p.up.we2,# + theme(legend.position = "none"),#, axis.text.x = element_blank()),
                              align = "h")
pang.up

# Tilapia ponds----
pseq_tila <- ps %>% 
  ps_filter(#Reported_disease != "Diseased",
    Crop == "Tilapia",
    Crop_species != "Gulsha.Carp",
    Crop_species != "Gulsha.Pabda")
#Crop_species != "T.no.fish") %>% 
#tax_fix() #%>% 
pseq_tila # 8257 taxa and 371 samples

## A2. Alpha diversity----
ps_rarefy_2.2 <- rarefy_even_depth(pseq_tila, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#1716OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_2.2 <- estimate_richness(ps_rarefy_2.2, measures = c("Chao1", "Shannon"))
alpha_estimates_2.2 <- cbind(alpha_estimates_2.2, sample_data(pseq_tila))

# Upazila_Tilapia ponds
### Shannon----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
t.up.s <- ggplot(data = alpha_estimates_2.2, aes(y = Shannon, x = Upazila)) +
  geom_boxplot(aes(color = Upazila)) +
  geom_jitter(aes(color = Upazila), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = upazila.colors,
                     labels = c("Jamalpur_Sadar" = "Jamalpur Sadar", "Tarakanda" = "Tarakanda")) +
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
    legend.text = element_text(size = 18),
    panel.grid = element_blank()) +
  guides(fill = guide_legend(nrow = 1)) + # Set number of columns in the legend
  labs (y = "Shannon", # Add y-axis title to the plot
        x = "Upazila")+
  #title = "d") + 
  scale_x_discrete(labels = c("Jamalpur_Sadar" = "Jamalpur Sadar")) +
  ylim(2.5, 6)
t.up.s

#### Perform Kw test----
kruskal.test(Shannon ~ Upazila, data = alpha_estimates_2.2)
#data:  Shannon by Upazila
#Kruskal-Wallis chi-squared = 5.0478, df = 1, p-value = 0.02466

# perform Dunn test (https://www.youtube.com/watch?v=pyLQmUfrel8)
stat.test.tus <- dunn_test(Shannon ~ Upazila, data = alpha_estimates_2.2, p.adjust.method = "BH")

# add the value on the plotBox plot
stat.test.tus <- stat.test.tus %>% add_xy_position(x = "Upazila")
t.up.s2 <- t.up.s + stat_pvalue_manual(stat.test.tus, label = "p.adj.signif", 
                                       step.increase = 0.05, tip.length = 0.005, 
                                       size = 8)  # Adjust the size of the asterisk here
#hide.ns = TRUE)
t.up.s2

## B2. Beta diversity (weighted unifrac)----
# wunifrac:weighted UniFrac took about #04:33 - 04:45 = 12 minutes to calculate the distance
#t_wuni <- pseq_tila %>%
#  tax_transform(rank = "unique", trans = "compositional") %>%
#  dist_calc(dist = "wunifrac") 
#saveRDS(t_wuni, "weighted_unifrac_beta_diversity_for_upazila_tilapia.rds")
t_wuni <- readRDS("weighted_unifrac_beta_diversity_for_upazila_tilapia.rds")

# Upazila
t.up.we <- t_wuni %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Upazila", fill = "Upazila",
    shape = 16, 
    alpha = 1.5,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(colour = Upazila)
  ) +
  scale_color_manual(values = upazila.colors,
                     labels = c("Jamalpur_Sadar" = "Jamalpur Sadar", "Tarakanda" = "Tarakanda")) +
  scale_fill_manual(values = upazila.colors,
                    labels = c("Jamalpur_Sadar" = "Jamalpur Sadar", "Tarakanda" = "Tarakanda")) +
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
    legend.text = element_text(size = 18),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  guides(fill = guide_legend(nrow = 1)) #+ # Set number of columns in the legend labs( #x = "PCoA1 [25.7%]",
#y = "PCoA2 [16.2%]",
#title = "e")
t.up.we

#### PERMANOVA ----
metadata_tila <- sample_data(pseq_tila) %>%
  data.frame() %>%
  tibble()
wuni_dist.t = phyloseq::distance(pseq_tila, method="wunifrac")
PERM_up_wuni.t <- adonis2(wuni_dist.t ~ Upazila, data = metadata_tila)

# to add test result on the plot
t.up.we2 <- t.up.we + annotate(geom = "label",
                               label = paste("PERMANOVA: R² = ", round(PERM_up_wuni.t["Upazila","R2"], 3), 
                                             ", p = ", PERM_up_wuni.t["Upazila", "Pr(>F)"], sep = ""),
                               x=Inf, y=Inf, hjust = 1, vjust = 1)
t.up.we2

# Combine tilapia----
tila.up <- cowplot::plot_grid(t.up.s2 + theme(legend.position = "none"),# axis.text.x = element_blank()),
                              t.up.we2,# + theme(legend.position = "none"), #axis.text.x = element_blank()),
                              align = "h")
tila.up

# 2. Microeukaryotes----
# Load the 18s data
ps.18s <- readRDS("phyloseq_18S_filtered_with_tree_pr2_90-150bp_20240416.rds")
ps.18s # 5390 taxa and 872 samples

# Pangasius only----
pseq_pang.18s <- ps.18s %>%
  ps_filter(
    #Reported_disease == "Non-diseased",
    Crop == "Pangasius" & Crop_species != "Shing")
#Crop_species != "P.no.fish") #%>% 
#tax_fix()# %>% 
pseq_pang.18s # 3947 taxa and 421 samples

## A3. Alpha diversity----
set.seed(1234)#The set. seed() function in R is used to create reproducible results when writing code that involves creating variables that take on random values. By using the set. seed() function, you guarantee that the same random values are produced each time you run the code
ps_rarefy_2.1.18s <- rarefy_even_depth(pseq_pang.18s, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#285OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_2.1.18s <- estimate_richness(ps_rarefy_2.1.18s, measures = c("Chao1", "Shannon"))
alpha_estimates_2.1.18s <- cbind(alpha_estimates_2.1.18s, sample_data(pseq_pang.18s))

## upazila
### Shannon----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
p.up.s.18s <- ggplot(data = alpha_estimates_2.1.18s, aes(y = Shannon, x = Upazila)) +
  geom_boxplot(aes(color = Upazila)) +
  geom_jitter(aes(color = Upazila), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = upazila.colors,
                     labels = c("Jamalpur_Sadar" = "Jamalpur Sadar")) +
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
    legend.text = element_text(size = 18),
    panel.grid = element_blank()) +
  guides(fill = guide_legend(nrow = 1)) + # Set number of columns in the legend
  labs (y = "Shannon", # Add y-axis title to the plot
        x = "Upazila")+
  #title = "a") +
  scale_x_discrete(labels = c("Jamalpur_Sadar" = "Jamalpur Sadar")) + # Insert line break in category names
  ylim(0, 6)
p.up.s.18s

#### Perform Kw test----
kruskal.test(Shannon ~ Upazila, data = alpha_estimates_2.1.18s)
#data:  Shannon by Upazila
#Kruskal-Wallis chi-squared = 3.1727, df = 1, p-value = 0.07488

##### Perform Dunn test (https://www.youtube.com/watch?v=pyLQmUfrel8)----
stat.test.pvs.18s <- dunn_test(Shannon ~ Upazila, data = alpha_estimates_2.1.18s, p.adjust.method = "BH")
# actually for two variables i don't need dunn test but to automatically add the value on 
# plot, i did it. It doesn't change the value

# add the value on the plotBox plot
stat.test.pvs.18s <- stat.test.pvs.18s %>% add_xy_position(x = "Upazila")
p.up.s2.18s <- p.up.s.18s + stat_pvalue_manual(stat.test.pvs.18s, label = "p.adj.signif", 
                                               step.increase = 0.05, tip.length = 0.005, 
                                               size = 8,  # Adjust the size of the asterisk here
                                               hide.ns = TRUE)
p.up.s2.18s

## B3. Beta diversity (weighted unifrac)----
# wunifrac:weighted UniFrac took about #03:10 - 03:21 = 11 minutes to calculate the distance
#up_wuni.18s <- pseq_pang.18s %>%
#  tax_transform(rank = "unique", trans = "compositional") %>%
#  dist_calc(dist = "wunifrac") 
#saveRDS(up_wuni.18s, "weighted_unifrac_beta_diversity_for_upazila_pangasius.18s.rds")
up_wuni.18s <- readRDS("weighted_unifrac_beta_diversity_for_upazila_pangasius.18s.rds")

# Upazila_pangasius pond
p.up.we.18s <- up_wuni.18s %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Upazila", 
    fill = "Upazila",
    shape = 16, 
    alpha = 1,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(colour = Upazila)
  ) +
  scale_color_manual(values = upazila.colors,
                     labels = c("Jamalpur_Sadar" = "Jamalpur Sadar")) +
  scale_fill_manual(values = upazila.colors,
                    labels = c("Jamalpur_Sadar" = "Jamalpur Sadar")) +
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
    legend.text = element_text(size = 18),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  guides(fill = guide_legend(nrow = 1)) #+ # Set number of columns in the legend
#labs( #x = "PCoA1 [25.7%]",
#y = "PCoA2 [16.2%]",
#  title = "b")
p.up.we.18s

### PERMANOVA----
metadata_pang.18s <- sample_data(pseq_pang.18s) %>%
  data.frame() %>%
  tibble()
wuni_dist.18s = phyloseq::distance(pseq_pang.18s, method="wunifrac")
PERM_up_wuni.18s <- adonis2(wuni_dist.18s ~ Upazila, data = metadata_pang.18s)

# to add test result on the plot
p.up.we2.18s <- p.up.we.18s + annotate(geom = "label",
                                       label = paste("PERMANOVA: R² = ", round(PERM_up_wuni.18s["Upazila","R2"], 3), 
                                                     ", p = ", PERM_up_wuni.18s["Upazila", "Pr(>F)"], sep = ""),
                                       x=Inf, y=Inf, hjust = 1, vjust = 1)
p.up.we2.18s

# Combine pangasius 18s----
pang.up.18s <- cowplot::plot_grid(p.up.s2.18s + theme(legend.position = "none"),# axis.text.x = element_blank()),
                                  p.up.we2.18s,# + theme(legend.position = "none"),# axis.text.x = element_blank()),
                                  align = "h")

# Tilapia ponds----
pseq_tila.18s <- ps.18s %>% 
  ps_filter(#Reported_disease != "Diseased",
    Crop == "Tilapia",
    Crop_species != "Gulsha.Carp",
    Crop_species != "Gulsha.Pabda")
#Crop_species != "T.no.fish") %>% 
#tax_fix() #%>% 
pseq_tila.18s # 4463 taxa and 376 samples

## A4. Alpha diversity----
ps_rarefy_2.2.18s <- rarefy_even_depth(pseq_tila.18s, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#287OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_2.2.18s <- estimate_richness(ps_rarefy_2.2.18s, measures = c("Chao1", "Shannon"))
alpha_estimates_2.2.18s <- cbind(alpha_estimates_2.2.18s, sample_data(pseq_tila.18s))

# Upazila_Tilapia ponds
### Shannon----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
t.up.s.18s <- ggplot(data = alpha_estimates_2.2.18s, aes(y = Shannon, x = Upazila)) +
  geom_boxplot(aes(color = Upazila)) +
  geom_jitter(aes(color = Upazila), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = upazila.colors,
                     labels = c("Jamalpur_Sadar" = "Jamalpur Sadar", "Tarakanda" = "Tarakanda")) +
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
    legend.text = element_text(size = 18),
    panel.grid = element_blank()) +
  guides(fill = guide_legend(nrow = 1)) + # Set number of columns in the legend
  labs (y = "Shannon", # Add y-axis title to the plot
        x = "Upazila")+
  #title = "d") + 
  scale_x_discrete(labels = c("Jamalpur_Sadar" = "Jamalpur Sadar")) +
  ylim(0, 6)
t.up.s.18s

#### Perform Kw test----
kruskal.test(Shannon ~ Upazila, data = alpha_estimates_2.2.18s)
#data:  Shannon by Upazila
#Kruskal-Wallis chi-squared = 28.2, df = 1, p-value = 1.094e-07

# perform Dunn test (https://www.youtube.com/watch?v=pyLQmUfrel8)
stat.test.tus.18s <- dunn_test(Shannon ~ Upazila, data = alpha_estimates_2.2.18s, p.adjust.method = "BH")

# add the value on the plotBox plot
stat.test.tus.18s <- stat.test.tus.18s %>% add_xy_position(x = "Upazila")
t.up.s2.18s <- t.up.s.18s + stat_pvalue_manual(stat.test.tus.18s, label = "p.adj.signif", 
                                               step.increase = 0.05, tip.length = 0.005, 
                                               size = 8)  # Adjust the size of the asterisk here
#hide.ns = TRUE)
t.up.s2.18s


## B4. Beta diversity (weighted unifrac)----
# wunifrac:weighted UniFrac took about #04:33 - 04:45 = 12 minutes to calculate the distance
#t_wuni.18s <- pseq_tila.18s %>%
#  tax_transform(rank = "unique", trans = "compositional") %>%
#  dist_calc(dist = "wunifrac") 
#saveRDS(t_wuni.18s, "weighted_unifrac_beta_diversity_for_upazila_tilapia.18s.rds")
t_wuni.18s <- readRDS("weighted_unifrac_beta_diversity_for_upazila_tilapia.18s.rds")

# Upazila
t.up.we.18s <- t_wuni.18s %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Upazila", fill = "Upazila",
    shape = 16, 
    alpha = 1.5,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(colour = Upazila)
  ) +
  scale_color_manual(values = upazila.colors,
                     labels = c("Jamalpur_Sadar" = "Jamalpur Sadar", "Tarakanda" = "Tarakanda")) +
  scale_fill_manual(values = upazila.colors,
                    labels = c("Jamalpur_Sadar" = "Jamalpur Sadar", "Tarakanda" = "Tarakanda")) +
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
    legend.text = element_text(size = 18),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  guides(fill = guide_legend(nrow = 1)) #+ # Set number of columns in the legend labs( #x = "PCoA1 [25.7%]",
#y = "PCoA2 [16.2%]",
#title = "e")
t.up.we.18s

#### PERMANOVA ----
metadata_tila.18s <- sample_data(pseq_tila.18s) %>%
  data.frame() %>%
  tibble()
wuni_dist.t.18s = phyloseq::distance(pseq_tila.18s, method="wunifrac")
PERM_up_wuni.t.18s <- adonis2(wuni_dist.t.18s ~ Upazila, data = metadata_tila.18s)

# to add test result on the plot
t.up.we2.18s <- t.up.we.18s + annotate(geom = "label",
                                       label = paste("PERMANOVA: R² = ", round(PERM_up_wuni.t.18s["Upazila","R2"], 3), 
                                                     ", p = ", PERM_up_wuni.t.18s["Upazila", "Pr(>F)"], sep = ""),
                                       x=Inf, y=Inf, hjust = 1, vjust = 1)
t.up.we2.18s

# Combine tilapia 18s----
tila.up.18s <- cowplot::plot_grid(t.up.s2.18s + theme(legend.position = "none"),# axis.text.x = element_blank()),
                                  t.up.we2.18s,# + theme(legend.position = "none"), #axis.text.x = element_blank()),
                                  align = "h")
tila.up.18s

# Same plot but without legend
pang.up2 <- cowplot::plot_grid(p.up.s2 + theme(legend.position = "none"), 
                               #axis.text.x = element_blank()),
                               p.up.we2 + theme(legend.position = "none"),#, axis.text.x = element_blank()),
                               align = "h")

pang.up.18s2 <- cowplot::plot_grid(p.up.s2.18s + theme(legend.position = "none"),# axis.text.x = element_blank()),
                                   p.up.we2.18s + theme(legend.position = "none"),# axis.text.x = element_blank()),
                                   align = "h")

tila.up2 <- cowplot::plot_grid(t.up.s2 + theme(legend.position = "none"),# axis.text.x = element_blank()),
                               t.up.we2 + theme(legend.position = "none"), #axis.text.x = element_blank()),
                               align = "h")

tila.up.18s2 <- cowplot::plot_grid(t.up.s2.18s + theme(legend.position = "none"),# axis.text.x = element_blank()),
                                   t.up.we2.18s + theme(legend.position = "none"), #axis.text.x = element_blank()),
                                   align = "h")

# Combine two upazilas----
comp.pt2 <- cowplot::plot_grid(pang.up2,
                               pang.up.18s2,
                               tila.up2,
                               tila.up.18s2,
                               labels = "AUTO",
                               ncol = 1)
comp.pt2
# Save as svg file so that can edit the image using software like inkscape
#ggsave("Effect of Upazila-pt.18s.svg", combined_plot.up.pt.18s2, width = 20, height = 12)

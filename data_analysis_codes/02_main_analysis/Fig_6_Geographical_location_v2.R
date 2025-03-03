# 02/10/2024
# Sanjit Debnath
# This script is to plot microbial composition (prokaryotes and microeukaryotes), alpha, beta diversity, and shared ASVs  to see the effect of different geographical locations

# Load Libraries
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(microViz); packageVersion("microViz")
library(vegan); packageVersion("vegan") # needed for PERMANOVA test
library(pairwiseAdonis); packageVersion("pairwiseAdonis") # to perform pairwise comparisom (https://github.com/pmartinezarbizu/pairwiseAdonis)
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
library(phylosmith); packageVersion("phylosmith")

## Setup working dictionary first
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/statistical_analysis")

# Set theme for ggplot
theme_set(theme_bw())
set.seed(1234)

# Define colors for each upazila
upazila.colors <- c("Jamalpur_Sadar" = "#1B9E77FF", "Muktagacha" = "#E6AB02FF", "Tarakanda" = "#490092FF")
upazila.colors2 <- c("JS" = "#1B9E77FF", "MG" = "#E6AB02FF", "TK" = "#490092FF")

# Load the phyloseq object
ps <- readRDS("phyloseq_metadata_6_v3_20240419.rds")
ps #10523 taxa and 891 samples

# Prokaryotics----
# Pangasius and tilapia combined
pseq_pt <- ps %>%
  ps_filter(
    Crop_species != "Shing",
    Crop_species != "Gulsha.Carp" & Crop_species != "Gulsha.Pabda"
  ) %>%
  ps_mutate(
    Upazila = recode(Upazila, # recode from dplyr package to replace the existing names
                     "Muktagacha" = "MG",
                     "Tarakanda" = "TK",
                     "Jamalpur_Sadar" = "JS"))
pseq_pt # 10378 taxa and 812 samples

## Estimate alpha diversity----
set.seed(1234)#The set. seed() function in R is used to create reproducible results when writing code that involves creating variables that take on random values. By using the set. seed() function, you guarantee that the same random values are produced each time you run the code
ps_rarefy_2 <- rarefy_even_depth(pseq_pt, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#1318OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_2 <- estimate_richness(ps_rarefy_2, measures = c("Chao1", "Shannon"))
alpha_estimates_2 <- cbind(alpha_estimates_2, sample_data(pseq_pt))


### Shannon----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
pt.up.s <- ggplot(data = alpha_estimates_2, aes(y = Shannon, x = Upazila)) +
  geom_boxplot(aes(color = Upazila)) +
  geom_jitter(aes(color = Upazila), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = upazila.colors2) +
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
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  guides(fill = guide_legend(nrow = 1)) + # Set number of columns in the legend
  labs (fill = "Upazila", # to add customized legend title
        y = "Shannon", # Add y-axis title to the plot
        x = "Upazila")+
  #title = "b") +
  #scale_x_discrete(labels = c("Jamalpur_Sadar" = "Jamalpur Sadar")) + # Insert line break in category names
  ylim(2.5, 7)
pt.up.s

#### Perform Kw test----
kruskal.test(Shannon ~ Upazila, data = alpha_estimates_2)
#data:  Shannon by Upazila
#Kruskal-Wallis chi-squared = 50.229, df = 2, p-value = 1.239e-11

stat.test.pts <- dunn_test(Shannon ~ Upazila, data = alpha_estimates_2, p.adjust.method = "BH")

# add the value on the plotBox plot
stat.test.pts <- stat.test.pts %>% add_xy_position(x = "Upazila")
pt.up.s2 <- pt.up.s + stat_pvalue_manual(stat.test.pts, label = "p.adj.signif", 
                                         step.increase = 0.05, tip.length = 0.005, 
                                         size = 8)  # Adjust the size of the asterisk here
#hide.ns = TRUE)
pt.up.s2

## Beta diversity (weighted unifrac)----
# wunifrac:weighted UniFrac took about #04:33 - 04:45 = 12 minutes to calculate the distance
#pt_wuni <- pseq_pt %>%
#  tax_transform(rank = "unique", trans = "compositional") %>%
#  dist_calc(dist = "wunifrac") 
#saveRDS(pt_wuni, "weighted_unifrac_beta_diversity_for_upazila_pangasius_tilapia.rds")
pt_wuni <- readRDS("weighted_unifrac_beta_diversity_for_upazila_pangasius_tilapia.rds")

# Upazila
pt.up.we <- pt_wuni %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Upazila", fill = "Upazila",
    shape = 16, 
    alpha = 1,
    size =3
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
    legend.text = element_text(size = 14),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  guides(fill = guide_legend(nrow = 1)) #+ # Set number of columns in the legendlabs( #x = "PCoA1 [25.7%]",
pt.up.we

### PERMANOVA----
metadata_pt <- sample_data(pseq_pt) %>%
  data.frame() %>%
  tibble()
wuni_dist.pt = phyloseq::distance(pseq_pt, method="wunifrac")
PERM_up_wun.pt <- adonis2(wuni_dist.pt ~ Upazila, data = metadata_pt)

# to add test result on the plot
pt.up.we2 <- pt.up.we + annotate(geom = "label",
                                 label = paste("PERMANOVA: R² = ", round(PERM_up_wun.pt["Upazila","R2"], 3), 
                                               ", p = ", PERM_up_wun.pt["Upazila", "Pr(>F)"], sep = ""),
                                 x=Inf, y=Inf, hjust = 1, vjust = 1)
pt.up.we2

##### Pairwise adonis----
pairwise_adonis_result <- pairwise.adonis(wuni_dist.pt, metadata_pt$Upazila,  p.adjust.m = "BH")
pairwise_adonis_result

# Add both PERMANOVA and pairwise adonis results on the plot
pt.up.we3 <- pt.up.we + 
  annotate(
    geom = "label",
    label = paste(
      "PERMANOVA, R² =", round(PERM_up_wun.pt["Upazila","R2"], 3), 
      ", p =", PERM_up_wun.pt["Upazila", "Pr(>F)"], "\n",
      "MG vs TK: R² =", round(pairwise_adonis_result[1, "R2"], 3),
      ", FDR =", pairwise_adonis_result[1, "p.adjusted"], "\n",
      "MG vs JS: R² =", round(pairwise_adonis_result[2, "R2"], 3),
      ", FDR =", pairwise_adonis_result[2, "p.adjusted"], "\n",
      "TK vs JS: R² =", round(pairwise_adonis_result[3, "R2"], 3),
      ", FDR =", pairwise_adonis_result[3, "p.adjusted"]
    ),
    x = Inf, y = Inf, hjust = 1, vjust = 1)
pt.up.we3

## Shared taxa (Venn diagram)----
# Subset Jamalpur sadar
js.pt <- pseq_pt %>% 
  ps_filter(
    Upazila == "JS") %>% 
tax_fix()
js.pt # 7078 taxa and 277 samples

# Subset Muktagacha
muk.pt <- pseq_pt %>% 
  ps_filter(
    Upazila == "MG") %>% 
tax_fix()
muk.pt # 7418 taxa and 365 samples

# Subset Tarakanda
tara.pt <- pseq_pt %>% 
  ps_filter(Upazila == "TK")%>% 
  tax_fix()
tara.pt #5859 taxa and 170 samples

js.pt_asvs <- js.pt %>% otu_table() %>%  rownames() 
muk.pt_asvs <- muk.pt %>% otu_table() %>%  rownames()
tara.pt_asvs <- tara.pt %>% otu_table() %>%  rownames()

venn_list_up.pt <- list("JS" = js.pt_asvs, "MG" = muk.pt_asvs, "TK" = tara.pt_asvs)

# Plot
venn_plot_up.pt <- venn.diagram(venn_list_up.pt, filename = NULL,
                                #circles
                                lwd = 2, lty = 'blank', fill = c("#1B9E77FF", "#E6AB02FF", "#490092FF"),
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
pt <- grid.arrange(venn_plot_up.pt)
pt <- grid.draw(venn_plot_up.pt)

## Compositional barplot Geographical location----
# pangasius and tilapia combined
pseq_pt2 <- pseq_pt %>% 
  tax_fix()
pseq_pt2 # 10378 taxa and 812 samples

#### Family proportion-----
#### Between upazila-----
tax.abd.pt.up.fam <- taxa_proportions(pseq_pt2, 'Family') # Overall to check who are dominant 
tax.abd.pt.up2.fam <- taxa_proportions(pseq_pt2, 'Family', treatment = "Upazila")

# Family
# my palette
myPal.fam.16s <- tax_palette(
  data = pseq_pt2, rank = "Family", n = 35, pal = "brewerPlus",
  add = c(Other = "lightgrey"))
#myPal.fam.16s["Sporichthyaceae"] <- "#00D9FF" # Override existing color if any color is not good
#myPal.fam.16s["Ilumatobacteraceae"] <- "#ABE496"
#myPal.fam.16s["Microbacteriaceae"] <- "#EA70FF" 
myPal.fam.16s["Sporichthyaceae"] <-  "#CC00A7"
myPal.fam.16s["Clostridiaceae"] <- "#B2DF8A"
tax_palette_plot(myPal.fam.16s)

# set up for alphabetical sorting
topTaxa.pt.fam <- pseq_pt2 %>%
  tax_top(n = 10, rank = "Family") %>%
  sort() # this makes them alphabetical

## plot with alphabetical sorting
up.famTop10 <- pseq_pt2 %>%
  tax_sort(by = sum, at = "Family") %>% # this orders all genera by abundance
  ps_select(Crop_species, Upazila) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Upazila") %>%
  comp_barplot(tax_order = topTaxa.pt.fam, # this brings the named taxa to the front
    tax_level = "Family", n_taxa = 10, 
    tax_sort = "name",
    sample_order = "asis", 
    merge_other = TRUE, other_name = "Other",
    bar_width = 0.9,
    bar_outline_colour = NA,
    palette = myPal.fam.16s
  ) +
  #coord_flip() +
  theme(plot.title = element_text(hjust = 0, size = 18, margin = margin(b = 15)), # Increase margin for plot title
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_blank(),  # to remove x axis title
        axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
        legend.position = "right",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18, face = "italic")) +
  labs(title = "Prokaryotes: Upazilas",
       y = "Relative abundance") +
  #scale_x_discrete(labels = function(x) gsub("_", "\n", x)) +
  guides(fill = guide_legend(ncol = 1))  # Adjust the number of rows in the legend
up.famTop10

### Compare taxa between groups----
### Make a dataframe with dominant phyla
pseq_pt3.fam <- pseq_pt2 %>%
  tax_transform(trans = "compositional", rank = "Family") %>%
  ps_get() %>%
  ps_otu2samdat(taxa = c("Pirellulaceae","Cyanobiaceae","Sporichthyaceae", # RA > 5.19%
                         "Phycisphaeraceae","Chthoniobacteraceae", "Clostridiaceae", # RA > 3.2%          
                         "Isosphaeraceae", "Kapabacteriales Order","PeM15 Order",
                         "Pedosphaeraceae"
                         )) %>% 
  samdat_tbl()

# run kw and dunn test for dominant taxa
kruskal.test(Pirellulaceae ~ Upazila, data = pseq_pt3.fam) # p-value < 2.2e-16
dunn_test(Pirellulaceae ~ Upazila, p.adjust.method = "BH",data = pseq_pt3.fam)

kruskal.test(Cyanobiaceae ~ Upazila, data = pseq_pt3.fam) # p-value  < 2.2e-16
dunn_test(Cyanobiaceae ~ Upazila, p.adjust.method = "BH",data = pseq_pt3.fam)

kruskal.test(Sporichthyaceae ~ Upazila, data = pseq_pt3.fam) # p-value = 1.328e-10
dunn_test(Sporichthyaceae ~ Upazila, p.adjust.method = "BH",data = pseq_pt3.fam)

kruskal.test(Phycisphaeraceae ~ Upazila, data = pseq_pt3.fam) # p-value < 2.2e-16
dunn_test(Phycisphaeraceae ~ Upazila, p.adjust.method = "BH",data = pseq_pt3.fam)

kruskal.test(Chthoniobacteraceae ~ Upazila, data = pseq_pt3.fam) # p-value = 6.933e-07
dunn_test(Chthoniobacteraceae ~ Upazila, p.adjust.method = "BH",data = pseq_pt3.fam)

kruskal.test(Clostridiaceae ~ Upazila, data = pseq_pt3.fam) # p-value < 2.2e-16
dunn_test(Clostridiaceae ~ Upazila, p.adjust.method = "BH",data = pseq_pt3.fam)

#### Combine plots----
pt.up <- cowplot::plot_grid(up.famTop10 + theme(legend.position = "right") + guides(fill = guide_legend(ncol = 1)) ,
                            pt.up.s2 + theme(legend.position = "none"),
                            pt.up.we3 + theme(legend.position = "none"),
                            align = "h",
                            rel_widths = c(3.6,3.2,3.2),
                            labels = c("A","B","C D"),
                            ncol = 3)
pt.up2 <- grid.arrange(pt.up,
                       venn_plot_up.pt,
                       widths = c(0.75, 0.25),# in grid.arrange, it's widths instead of rel_widths
                       ncol = 2)  # Set the number of columns as per your preference

# Microeukaryotics----
# Load phyloseq object
ps.18s <- readRDS("phyloseq_18S_filtered_with_tree_pr2_90-150bp_20240416.rds")
ps.18s # 5390 taxa and 872 samples

# Remove samples without pangasius and tilapia
pseq_pt.18s <- ps.18s %>% 
  ps_filter(
    Crop_species != "Shing",
    Crop_species != "Gulsha.Carp" & Crop_species != "Gulsha.Pabda"
  ) %>%
  ps_mutate(
    Upazila = recode(Upazila,
                     "Muktagacha" = "MG",
                     "Tarakanda" = "TK",
                     "Jamalpur_Sadar" = "JS"))
pseq_pt.18s # 5210 taxa and 797 samples

## Estimate alpha diversity----
set.seed(1234)#The set. seed() function in R is used to create reproducible results when writing code that involves creating variables that take on random values. By using the set. seed() function, you guarantee that the same random values are produced each time you run the code
ps_rarefy_2.18s <- rarefy_even_depth(pseq_pt.18s, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#210OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_2.18s <- estimate_richness(ps_rarefy_2.18s, measures = c("Chao1", "Shannon"))
alpha_estimates_2.18s <- cbind(alpha_estimates_2.18s, sample_data(pseq_pt.18s))

### Shannon----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
pt.up.s.18s <- ggplot(data = alpha_estimates_2.18s, aes(y = Shannon, x = Upazila)) +
  geom_boxplot(aes(color = Upazila)) +
  geom_jitter(aes(color = Upazila), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = upazila.colors2) +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, face = "bold", margin = margin(t = 5)),
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
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  guides(fill = guide_legend(ncol = 1)) + # Set number of columns in the legend
  labs (fill = "Upazila", # to add customized legend title
        y = "Shannon", # Add y-axis title to the plot
        x = "Upazila")
#title = "a") +
#scale_x_discrete(labels = c("Jamalpur_Sadar" = "Jamalpur Sadar")) #+ # Insert line break in category names
#ylim(2.5, 6)
pt.up.s.18s

#### Perform Kw test----
kruskal.test(Shannon ~ Upazila, data = alpha_estimates_2.18s)
#data:  Shannon by Upazila
#Kruskal-Wallis chi-squared = 40.087, df = 2, p-value = 1.973e-09

stat.test.pts.18s <- dunn_test(Shannon ~ Upazila, data = alpha_estimates_2.18s, p.adjust.method = "BH")

# add the value on the plotBox plot
stat.test.pts.18s <- stat.test.pts.18s %>% add_xy_position(x = "Upazila")
pt.up.s2.18s <- pt.up.s.18s + stat_pvalue_manual(stat.test.pts.18s, label = "p.adj.signif", 
                                                 step.increase = 0.05, tip.length = 0.005, 
                                                 size = 8)  # Adjust the size of the asterisk here
#hide.ns = TRUE)
pt.up.s2.18s


## Beta diversity (weighted unifrac)----
# wunifrac:weighted UniFrac took about #04:33 - 04:45 = 12 minutes to calculate the distance
#pt_wuni.18s <- pseq_pt.18s %>%
#  tax_transform(rank = "unique", trans = "compositional") %>%
#  dist_calc(dist = "wunifrac") 
#saveRDS(pt_wuni.18s, "weighted_unifrac_beta_diversity_for_upazila_pangasius_tilapia_18s.rds")
pt_wuni.18s <- readRDS("weighted_unifrac_beta_diversity_for_upazila_pangasius_tilapia_18s.rds")

# Upazila
pt.up.we.18s <- pt_wuni.18s %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Upazila", fill = "Upazila",
    shape = 21, 
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
    legend.text = element_text(size = 14),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  guides(fill = guide_legend(nrow = 1))#+ # Set number of columns in the legend
pt.up.we.18s

### PERMANOVA----
metadata_pt.18s <- sample_data(pseq_pt.18s) %>%
  data.frame() %>%
  tibble()
wuni_dist.pt.18s = phyloseq::distance(pseq_pt.18s, method="wunifrac")
PERM_up_wun.pt.18s <- adonis2(wuni_dist.pt.18s ~ Upazila, data = metadata_pt.18s)

# to add test result on the plot
pt.up.we2.18s <- pt.up.we.18s + annotate(geom = "label",
                                         label = paste("PERMANOVA: R² = ", round(PERM_up_wun.pt.18s["Upazila","R2"], 3), 
                                                       ", p = ", PERM_up_wun.pt.18s["Upazila", "Pr(>F)"], sep = ""),
                                         x=Inf, y=Inf, hjust = 1, vjust = 1)
pt.up.we2.18s

##### Pairwise adonis----
pairwise_adonis_result.18s <- pairwise.adonis(wuni_dist.pt.18s, metadata_pt.18s$Upazila,  p.adjust.m = "BH")
pairwise_adonis_result.18s

# Add both PERMANOVA and pairwise adonis results on the plot
pt.up.we3.18s <- pt.up.we.18s + 
  annotate(
    geom = "label",
    label = paste(
      "PERMANOVA, R² =", round(PERM_up_wun.pt.18s["Upazila","R2"], 3), 
      ", p =", PERM_up_wun.pt.18s["Upazila", "Pr(>F)"], "\n",
      "MG vs TK: R² =", round(pairwise_adonis_result.18s[1, "R2"], 3),
      ", FDR =", pairwise_adonis_result.18s[1, "p.adjusted"], "\n",
      "MG vs JS: R² =", round(pairwise_adonis_result.18s[2, "R2"], 3),
      ", FDR =", pairwise_adonis_result.18s[2, "p.adjusted"], "\n",
      "TK vs JS: R² =", round(pairwise_adonis_result.18s[3, "R2"], 3),
      ", FDR =", pairwise_adonis_result.18s[3, "p.adjusted"]
    ),
    x = Inf, y = Inf, hjust = 1, vjust = 1)
pt.up.we3.18s

## Shared taxa (Venn diagram)----
# Subset Jamalpur sadar
js.pt.18s <- pseq_pt.18s %>% 
  ps_filter(
    Upazila == "JS") %>% 
tax_fix()
js.pt.18s # 3933 taxa and 268 samples

# Subset Muktagacha
muk.pt.18s <- pseq_pt.18s %>% 
  ps_filter(
    Upazila == "MG") %>% 
tax_fix()
muk.pt.18s # 3633 taxa and 367 samples 

# Subset tarakanda
tara.pt.18s <- pseq_pt.18s %>% 
  ps_filter(Upazila == "TK")%>% 
  tax_fix()
tara.pt.18s # 3173 taxa and 162 samples

js.pt_asvs.18s <- js.pt.18s %>% otu_table() %>%  rownames() # already did before
muk.pt_asvs.18s <- muk.pt.18s %>% otu_table() %>%  rownames()
tara.pt_asvs.18s <- tara.pt.18s %>% otu_table() %>%  rownames()

venn_list_up.pt.18s <- list("JS" = js.pt_asvs.18s, "MG" = muk.pt_asvs.18s, "TK" = tara.pt_asvs.18s)

# Plot
venn_plot_up.pt.18s <- venn.diagram(venn_list_up.pt.18s, filename = NULL,
                                    #circles
                                    lwd = 2, lty = 'blank', fill = c("#1B9E77FF", "#E6AB02FF", "#490092FF"),
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
                                    auto_scale = FALSE)#, # it resize the circle
#main = "c") # main title set to empty string

# Create a new plot
plot.new()

# Draw the Venn diagram
pt.18s <- grid.arrange(venn_plot_up.pt.18s)
pt.18s <- grid.draw(venn_plot_up.pt.18s)

## Compositional barplot----
# pangasius and tilapia combined
pseq_pt.18s2 <- pseq_pt.18s %>% 
  tax_fix()
pseq_pt.18s2 # 5210 taxa and 797 samples

#### Family proportion-----
#### Between upazila-----
tax.abd.pt.18s.up.fam <- taxa_proportions(pseq_pt.18s2, 'Family') # Overall to check who are dominant 
tax.abd.pt.18s.up2.fam <- taxa_proportions(pseq_pt.18s2, 'Family', treatment = "Upazila")

# my palette
myPal.fam.18s <- tax_palette(
  data = pseq_pt.18s2, rank = "Family", n = 25, pal = "greenArmytage", #pal = "brewerPlus",
  add = c(Other = "lightgrey")) 
myPal.fam.18s["Alveolata Division"] <- "orange"
myPal.fam.18s["Euglenaceae"] <- "#FFE100"  
myPal.fam.18s["Raphidophyceae_XX"] <- "#5EF1F2"
myPal.fam.18s["Stephanodiscaceae"] <- "#00998F"
myPal.fam.18s["Gyrista Subdivision"] <- "#F0A3FF"
tax_palette_plot(myPal.fam.18s)

# set up for alphabetical sorting
topTaxa.pt.fam.18s <- pseq_pt.18s2 %>%
  tax_top(n = 10, rank = "Family") %>%
  sort() # this makes them alphabetical

## plot with alphabetical sorting
up.fam.18s <- pseq_pt.18s2 %>%
  tax_sort(by = sum, at = "Family") %>% # this orders all genera by abundance
  ps_select(Crop_species, Upazila) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "Upazila") %>%
  comp_barplot(tax_order = topTaxa.pt.fam.18s, # this brings the named taxa to the front
               tax_level = "Family", n_taxa = 10, 
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.fam.18s
  ) +
  #coord_flip() +
  theme(plot.title = element_text(hjust = 0, size = 18, margin = margin(b = 15)), # Increase margin for plot title
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_blank(),  # to remove x axis title
        axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
        legend.position = "right",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18, face = "italic")) +
  labs(title = "Microeukaryotes: Upazilas",
       y = "Relative abundance") +
  #scale_x_discrete(labels = function(x) gsub("_", "\n", x)) +
  guides(fill = guide_legend(ncol = 1))  # Adjust the number of rows in the legend
up.fam.18s

### Compare taxa between groups----
### Make a dataframe with dominant phyla
pseq_pt.18s3 <- pseq_pt.18s2 %>%
  tax_transform(trans = "compositional", rank = "Family") %>%
  ps_get() %>%
  ps_otu2samdat(taxa = c("Stephanodiscaceae","Eukaryota Domain","Gyrista Subdivision", # RA > 4.37%
                         "Fungi Subdivision","Alveolata Division","Dinophyceae Class", 
                         "Raphidophyceae_XX", "Euglenaceae", "Cryptophyceae Class",
                         "Cryptomonadales_X")) %>% 
  samdat_tbl()

# run kw and dunn test for dominant taxa
kruskal.test(Stephanodiscaceae ~ Upazila, data = pseq_pt.18s3) # p-value = 2.364e-09
dunn_test(Stephanodiscaceae ~ Upazila, p.adjust.method = "BH",data = pseq_pt.18s3)

kruskal.test(Euglenaceae ~ Upazila, data = pseq_pt.18s3) # p-value < 2.2e-16
dunn_test(Euglenaceae ~ Upazila, p.adjust.method = "BH",data = pseq_pt.18s3)

#### Combine plots----
pt.up.18s <- cowplot::plot_grid(up.fam.18s + theme(legend.position = "right") + guides(fill = guide_legend(ncol = 1)) ,
                                #up.phy, # + theme(legend.position = "none"),
                                pt.up.s2.18s + theme(legend.position = "none", axis.title.x = element_blank()),
                                pt.up.we3.18s + theme(legend.position = "none"),
                                align = "h",
                                rel_widths = c(3.6,3.2,3.2),
                                labels = c("E","F","G H"),
                                ncol = 3)
pt.up.18s2 <- grid.arrange(pt.up.18s,
                           venn_plot_up.pt.18s,
                           widths = c(0.75, 0.25),# in grid.arrange, it's widths instead of rel_widths
                           ncol = 2)  # Set the number of columns as per your preference

# Combine prokaryotic and microeukaryotic----
comb.pt.up <- cowplot::plot_grid(pt.up2,
                                 pt.up.18s2,
                                 ncol = 1)
comb.pt.up
# Save as 2000*1200 (and give final touch in inkscape)


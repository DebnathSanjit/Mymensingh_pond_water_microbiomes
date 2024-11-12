# 22/09/2024
# Sanjit Debnath

# This script is for microbial (prokaryotic and microeukaryotic) seasonal dynamics including shannon diversity, glm and harmonic regression at family rank with a relative abundance > 1%
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
library(rstatix); packageVersion("rstatix")
#library(dunn.test); packageVersion("dunn.test")
library(openxlsx) # to write workbook in excel
library(emmeans); packageVersion("emmeans") # for gamma glm
library(reshape); packageVersion("reshape")

## Setup working dictionary first
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/statistical_analysis")

# set seed
set.seed(1234)

# Set theme for ggplot
theme_set(theme_bw())

# Prokaryotes----
# Alpha diversity----
# Load the phyloseq object
ps <- readRDS("phyloseq_metadata_6_v3_20240419.rds")
ps # 10523 taxa and 891 samples

# Set season color
season.colors <- c("Monsoon" = "#EC823C","Winter" ="#1B9E77","Pre-monsoon" = "#100AFF")

# Manually define the correct order of months
month_order <- c("Jun-16","Oct-16","Nov-16","Dec-16","Jan-17","Feb-17","Mar-17",
                 "Apr-17","May-17","Jun-17","Jul-17","Aug-17","Sep-17","Oct-17", 
                 "Nov-17","Dec-17","Jan-18","Feb-18","Mar-18","Apr-18","May-18")

# Access the sample data from the phyloseq object and modify the Sampling_months column
sample_data(ps)$Sampling_months <- factor(sample_data(ps)$Sampling_months, levels = month_order)

# Pangasius ponds only----
pseq_pang <- ps %>% 
  ps_filter(
            Crop == "Pangasius" & Crop_species != "Shing") %>% 
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
pseq_pang # 8091 taxa and 421 samples

# Reorder season
season_order2 <- c("Monsoon", "Winter", "Pre-monsoon")
season_order3 <- c("Monsoon1", "Winter1", "Pre_monsoon1","Monsoon2", "Winter2", "Pre_monsoon2")

# Access the sample data from the phyloseq object and modify the Sampling_months column
sample_data(pseq_pang)$Season2 <- factor(sample_data(pseq_pang)$Season2, levels = season_order2)
sample_data(pseq_pang)$Season3 <- factor(sample_data(pseq_pang)$Season3, levels = season_order3)

## A1.Alpha diversity----
#set.seed(1234)#The set. seed() function in R is used to create reproducible results when writing code that involves creating variables that take on random values. By using the set. seed() function, you guarantee that the same random values are produced each time you run the code
ps_rarefy_2.1 <- rarefy_even_depth(pseq_pang, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#1346OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_2.1 <- estimate_richness(ps_rarefy_2.1, measures = c("Chao1", "Shannon"))
alpha_estimates_2.1 <- cbind(alpha_estimates_2.1, sample_data(pseq_pang))

### Shannon----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
p.sm.s <- ggplot(data = alpha_estimates_2.1, aes(y = Shannon, x = Sampling_months)) +
  geom_boxplot(aes(color = Season2)) +
  geom_jitter(aes(color = Season2), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = season.colors) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 18, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 18, margin = margin(r = 5)),
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
  guides(fill = guide_legend(ncol = 1)) + # Set number of columns in the legend
  labs (color = "Season",  # Customize legend title
    y = "Shannon", # Add y-axis title to the plot
    x = "Sampling months",
    title = "Prokaryotes: Pangasius ponds") +
  ylim(2.5, 6)
p.sm.s

### Statistical test----
# Since month is not an independept variable, the microbial community from one month can not be totally different from another month. Therefore, doing t-test (KW) is not enough. we need to do linear modelling (GLM, harmonic)

#### Check Normality of Shannon Diversity----
# Histogram
hist(alpha_estimates_2.1$Shannon, main="Histogram of Shannon Diversity", xlab="Shannon Index", col="lightblue", breaks=20)

# Q-Q plot
qqnorm(alpha_estimates_2.1$Shannon)
qqline(alpha_estimates_2.1$Shannon, col = "red")

# Shapiro-Wilk test
shapiro.test(alpha_estimates_2.1$Shannon)
# p-value < 2.2e-16 # not normally distributed

# Skewness and Kurtosis
library(moments); packageVersion("moments")
skewness(alpha_estimates_2.1$Shannon)
# -3.543162 # The skewness value of -3.54 suggests that the distribution is negatively skewed (a long tail on the left). which we can see from the histogram

kurtosis(alpha_estimates_2.1$Shannon)
# 36.14627 # The kurtosis value of 36.15 is much higher than 3, indicating that the distribution has heavy tails and is more peaked compared to a normal distribution.
# The Shannon index is not normally distributed. Since the normality assumption is violated, a standard Gaussian GLM may not be the best fit.

#### Gamma Distribution GLM (with log link)----
#Since the Shannon diversity is positive and skewed, a Gamma distribution with a log link might be more appropriate. This assumes that the Shannon index follows a gamma distribution, which handles positive, skewed data.
glm_gamma_pang <- glm(Shannon ~ Season3, family = Gamma(link = "log"), data = alpha_estimates_2.1)
glm_pang_summary <- summary(glm_gamma_pang)

# Extract coefficients
glm_coeff_pang <- as.data.frame(glm_pang_summary$coefficients)

# Add the row names (e.g., '(Intercept)', 'Season02.Winter', etc.) as a new column called 'Variable'
glm_coeff_pang$Variable <- rownames(glm_coeff_pang)

# Move 'Variable' to the first position
glm_coeff_pang <- glm_coeff_pang[, c("Variable", colnames(glm_coeff_pang)[1:4])]

# Format the p-values to avoid rounding off in Excel
#glm_coeff_pang$`Pr(>|t|)` <- format(glm_coeff_pang$`Pr(>|t|)`, scientific = TRUE)

#### Pairwise comparison----
pairwise_glm_pang <- emmeans(glm_gamma_pang, pairwise ~ Season3, adjust = "BH")
summary(pairwise_glm_pang)

# Extract pairwise comparison results
emmeans_table_pang <- as.data.frame(pairwise_glm_pang$emmeans)
contrasts_table_pang <- as.data.frame(pairwise_glm_pang$contrasts)

# Optional: You can also check diagnostic plots for the GLM
par(mfrow = c(2, 2))  # Plot 4 diagnostic plots
plot(glm_gamma_pang)

# Create a new Excel workbook
wb_pang <- createWorkbook()

# Add sheets for GLM summary, emmeans, and contrasts
addWorksheet(wb_pang, "GLM Summary")
addWorksheet(wb_pang, "EMMeans")
addWorksheet(wb_pang, "Contrasts")

# Write data to the respective sheets
writeData(wb_pang, sheet = "GLM Summary", glm_coeff_pang)
writeData(wb_pang, sheet = "EMMeans", emmeans_table_pang)
writeData(wb_pang, sheet = "Contrasts", contrasts_table_pang)
# Save the workbook to a file
#saveWorkbook(wb_pang, file = "glm_and_pairwise_results_pseq_pang.xlsx", overwrite = TRUE)

# Tilapia ponds only----
pseq_tila <- ps %>% 
  ps_filter(Crop == "Tilapia",
            Crop_species != "Gulsha.Carp" & Crop_species != "Gulsha.Pabda") %>% 
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
pseq_tila #  8380 taxa and 391 samples

# Reorder season
season_order2 <- c("Monsoon", "Winter", "Pre-monsoon")
season_order3 <- c("Monsoon1", "Winter1", "Pre_monsoon1","Monsoon2", "Winter2", "Pre_monsoon2")

# Access the sample data from the phyloseq object and modify the Sampling_months column
sample_data(pseq_tila)$Season2 <- factor(sample_data(pseq_tila)$Season2, levels = season_order2)
sample_data(pseq_tila)$Season3 <- factor(sample_data(pseq_tila)$Season3, levels = season_order3)

## A2. Alpha diversity----
ps_rarefy_2.2 <- rarefy_even_depth(pseq_tila, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#1716OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_2.2 <- estimate_richness(ps_rarefy_2.2, measures = c("Chao1", "Shannon"))
alpha_estimates_2.2 <- cbind(alpha_estimates_2.2, sample_data(pseq_tila))

### Shannon----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
t.sm.s <- ggplot(data = alpha_estimates_2.2, aes(y = Shannon, x = Sampling_months)) +
  geom_boxplot(aes(color = Season2)) +
  geom_jitter(aes(color = Season2), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = season.colors) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 18, margin = margin(t = 5)),
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
  guides(fill = guide_legend(nrow = 1)) + # Set number of columns in the legend
  labs (color = "Season",  # Customize legend title
    y = "Shannon", # Add y-axis title to the plot
    x = "Sampling months",
    title = "Prokaryotes: Tilapia ponds") +
  ylim(2.5, 6)
t.sm.s

### Statistical test----
#### Check Normality of Shannon Diversity----
# Histogram
hist(alpha_estimates_2.2$Shannon, main="Histogram of Shannon Diversity", xlab="Shannon Index", col="lightblue", breaks=20)

# Q-Q plot
qqnorm(alpha_estimates_2.2$Shannon)
qqline(alpha_estimates_2.2$Shannon, col = "red")

# Shapiro-Wilk test
shapiro.test(alpha_estimates_2.2$Shannon)
# p-value = 6.553e-07 # not normally distributed

# Skewness and Kurtosis
skewness(alpha_estimates_2.2$Shannon)
# -0.7369513 # The skewness value of -0.7369513 suggests that the distribution is negatively moderately left-skewed. which we can see from the histogram

kurtosis(alpha_estimates_2.2$Shannon)
#  4.341973 # The kurtosis value of  4.341973 is much higher than 3 (kurtosis value), indicating that the distribution has a sharper peak and heavier tails than a normal distribution.
# The Shannon index is not normally distributed. Since the normality assumption is violated, a standard Gaussian GLM may not be the best fit.

#### Gamma Distribution GLM (with log link)----
#Since the Shannon diversity is positive and skewed, a Gamma distribution with a log link might be more appropriate. This assumes that the Shannon index follows a gamma distribution, which handles positive, skewed data.
glm_gamma_tila <- glm(Shannon ~ Season3, family = Gamma(link = "log"), data = alpha_estimates_2.2)
glm_tila_summary <- summary(glm_gamma_tila)

# Extract coefficients
glm_coeff_tila <- as.data.frame(glm_tila_summary$coefficients)

# Add the row names (e.g., '(Intercept)', 'Season02.Winter', etc.) as a new column called 'Variable'
glm_coeff_tila$Variable <- rownames(glm_coeff_tila)

# Move 'Variable' to the first position
glm_coeff_tila <- glm_coeff_tila[, c("Variable", colnames(glm_coeff_tila)[1:4])]

# Format the p-values to avoid rounding off in Excel
glm_coeff_tila$`Pr(>|t|)` <- format(glm_coeff_tila$`Pr(>|t|)`, scientific = TRUE)

#### Pairwise comparison----
pairwise_glm_tila <- emmeans(glm_gamma_tila, pairwise ~ Season3, adjust = "BH")
summary(pairwise_glm_tila)

# Extract pairwise comparison results
emmeans_table_tila <- as.data.frame(pairwise_glm_tila$emmeans)
contrasts_table_tila <- as.data.frame(pairwise_glm_tila$contrasts)

# Optional: You can also check diagnostic plots for the GLM
par(mfrow = c(2, 2))  # Plot 4 diagnostic plots
plot(glm_gamma_tila)

# Create a new Excel workbook
wb_tila <- createWorkbook()

# Add sheets for GLM summary, emmeans, and contrasts
addWorksheet(wb_tila, "GLM Summary")
addWorksheet(wb_tila, "EMMeans")
addWorksheet(wb_tila, "Contrasts")

# Write data to the respective sheets
writeData(wb_tila, sheet = "GLM Summary", glm_coeff_tila)
writeData(wb_tila, sheet = "EMMeans", emmeans_table_tila)
writeData(wb_tila, sheet = "Contrasts", contrasts_table_tila)
# Save the workbook to a file
#saveWorkbook(wb_tila, file = "glm_and_pairwise_results_pseq_tila.xlsx", overwrite = TRUE)

# Harmonic regression----
# I have adapted this script from Luis to do a harmonic regression to see seasonal trend of taxa.
# As I have so many taxa, if I plot each ASV, the plot becomes very noisy and not very helpful. So, I aggregate at family level, adjusted p value using BH method 
# follow: https://github.com/lbolanos32/SAR11_BATS_WEC_2022/blob/main/Rscripts/Fig5_LinearModels.R

# Set family colors
family_colors <- c("Beijerinckiaceae" = "blue","Burkholderiaceae" = "#E6AB02",
                   "Chthoniobacteraceae" = "#CC00A7","Cyanobiaceae" = "#da6130",
                   "MWH-UniP1 aquatic group" = "#A6CEE3","Pirellulaceae" = "#1ff8ff",
                   "Pedosphaeraceae" = "#6A3D9A","Rubinisphaeraceae" = "#B15928",      
                   "Ilumatobacteraceae" = "#FFFF99", "Comamonadaceae" = "#66A61E" ,           
                   "Chitinophagaceae" = "#5f7b35",
                   "Planctomycetales Order" = "#EA70FF","Verrucomicrobiae Class" = "#5c47b8",
                   "Bacteria Kingdom" = "#7570B3", "Clostridiaceae" = "#33A02C", 
                   "Gammaproteobacteria Class" = "#E7298A","Mycobacteriaceae" = "#1F58B4",
                   "PeM15 Order" = "#9DCC00","Peptostreptococcaceae" = "#1B9E77",
                   "Rhizobiales Incertae Sedis" = "#000080", "Sporichthyaceae" = "#15A7E6")

### Pangasius ponds----
pseq_pang <- ps %>%
  ps_filter(Sampling_point > 1, # for seasonal trend, remove
            Crop == "Pangasius" & Crop_species != "Shing") %>% 
  tax_fix()
pseq_pang # 8044 taxa and 409 samples

##### convert reads to relative abundance
pang.rel <- transform_sample_counts(pseq_pang, function(x){x / sum(x)})

# Aggregate at Family rank
#pang.glom <- tax_glom(pang.rel, taxrank = "Family") # Since it takes long, save it for later use
#saveRDS(pang.glom, "pseq_pang_aggregated_at_family.rds")
pang.glom <- readRDS("pseq_pang_aggregated_at_family.rds")
pang.glom # 575 taxa and 409 samples

# Extract the taxonomic information
pang.tax <- as.data.frame(tax_table(pang.glom)[,5])
pang.tax$ASV <- rownames(pang.tax)

# Extract the ASV table
pang_ASV <- as.data.frame(otu_table(pang.glom))
# Transpose the ASV table
pang_ASVt <- t(pang_ASV) # transpose

#####Create a data.frame cols=ASV, rows= samples -> add date to this and estimate the sin/cos thing 
# Extract Family Taxonomy:
pang.tax <- as.data.frame(tax_table(pang.glom)[,5]) ###Keep this one to add to the melted object
# Convert ASV Table to Dataframe:
pang_ASV <- as.data.frame(otu_table(pang.glom))
# Transpose the ASV Table:
ASV_pang_trans <- as.data.frame(x = t(pang_ASV), stringsAsFactors = FALSE)
# Extract Metadata (Sampling Date):
md_to_add.pang <- as.data.frame(sample_data(pang.glom))[,c(12)] #column number for Sampling_date in metadata
# Combine ASV Data and Metadata:
final_2a.pang <- cbind(ASV_pang_trans,md_to_add.pang)
# Convert Sampling Dates to Date Format:
final_2a.pang$Sampling_date <- as.Date(final_2a.pang$Sampling_date,"%d/%m/%Y") # should i change as %Y/%m/%d
# Sort Data by Sampling Date:
final_2a.pang_sort <- final_2a.pang[order(final_2a.pang$Sampling_date),] ##sort from least recent to most recent
# Convert the Date Format Again:
final_2a.pang_sort$Sampling_date <- as.Date(final_2a.pang_sort$Sampling_date,"%Y-%m-%d")

##Add the Julian day # The Julian day starts from 1 for January 1st and goes up to 365 (or 366 in leap years) for December 31st.
final_2a.pang_sort$Jul <- format(final_2a.pang_sort$Sampling_date, "%j")

###Serial Day
#Starting from
pang_start <- as.Date('2016-01-01') # should be first day of the sampling year
final_2a.pang_sort$Serial <- as.numeric(-difftime(pang_start,final_2a.pang_sort$Sampling_date)+1)

##Days since winter solstice (21 December): This marks the shortest day and longest night of the year in Bangladesh, as well as the official beginning of winter in astronomical terms.
##Add the winter solstice 
final_2a.pang_sort$wintersols[final_2a.pang_sort$Sampling_date>as.Date("2016-01-01") & final_2a.pang_sort$Sampling_date<as.Date("2016-12-21")] <- "2015-12-21"
final_2a.pang_sort$wintersols[final_2a.pang_sort$Sampling_date>as.Date("2017-01-01") & final_2a.pang_sort$Sampling_date<as.Date("2017-12-21")] <- "2016-12-21"
final_2a.pang_sort$wintersols[final_2a.pang_sort$Sampling_date>as.Date("2018-01-01") & final_2a.pang_sort$Sampling_date<as.Date("2018-12-21")] <- "2017-12-21"

final_2a.pang_sort$DaysSincewintersols <- as.numeric(-difftime(final_2a.pang_sort$wintersols,final_2a.pang_sort$Sampling_date) + 1) #Add the numeric column

# Remove Rows with Missing Values:
final_2a.pang_sort <- final_2a.pang_sort[!is.na(final_2a.pang_sort$DaysSincewintersols), ]

## Convert Julian date to numeric:
final_2a.pang_sort$Jul <- as.numeric(final_2a.pang_sort$Jul)

# Generate sine and cosine terms for harmonic regression:
final_2a.pang_sort$cosDX1 <- as.numeric(cos(2*pi*((final_2a.pang_sort$DaysSincewintersols)/365))) #Peaking in midwinter
final_2a.pang_sort$sinDX1 <- as.numeric(sin(2*pi*((final_2a.pang_sort$Jul)/365))) #Peaking other time
# check final_2a.pang_sort for number of variables which will be needed below. Remove last 7 columns
#[576] "Sampling_date" [577] "Jul"[578]"Serial" [579]"wintersols" [580]"DaysSincewintersols" [581] "cosDX1" [582] "senDX1"

# Fit harmonic regression models:
storage.pang <- list()
for(i in names(final_2a.pang_sort)[1:575]){ #number of taxa present in final_2a.pang_sort (582) - 7
  storage.pang[[i]] <- lm(get(i) ~ sinDX1+cosDX1, final_2a.pang_sort)
}

# Extract p-values for significance testing:
pvalues_lm.pang <- sapply(storage.pang, function(x){
  ff <- summary(x)$fstatistic
  pf(ff[1], df1 = ff[2], df2 = ff[3], lower.tail = FALSE)
})

# Adjust the p-values using the BH method
pvalues_lm.pang_adjusted <- p.adjust(pvalues_lm.pang, method = "BH")

# Filter significant taxa based on adjusted p-values
pvalues_lm.pang_adjusted05 <- pvalues_lm.pang_adjusted[pvalues_lm.pang_adjusted < 0.05] # 143

# Extract significant taxa names after p-value adjustment
sign_names.pang <- str_remove(names(pvalues_lm.pang_adjusted05), ".value")

#Number of the column that matches a list of names, get the position of the names <0.05
signll.pang <- match(sign_names.pang,names(final_2a.pang_sort))

#subset list of lists 
storage05.pang <- storage.pang[signll.pang]

#Get the constant, sin and cos of the models:
# check if storage05.pang contains any NAs and if so, Filter out NA values from storage05.pang
filtered_storage.pang <- storage05.pang[!sapply(storage05.pang, is.null)]

# Apply sapply to filtered storage
coeffs05.pang <- data.frame(coeffs = sapply(filtered_storage.pang, function(item) item$coefficients))

#coeffs05.pang<-data.frame(coeffs=sapply(storage05.pang, FUN=function(item){item$coefficients}))
colnames(coeffs05.pang) <- str_remove(colnames(coeffs05.pang), "coeffs.")

#Subset original day, sin, cos
newdf_model.pang <- final_2a.pang_sort[, c("Serial","cosDX1","sinDX1")]

# Automatize this :P 
# use the function below to automatically do it, and add the value in newdf_model.pang table
automate_harmonic_regression.pang <- function(final_2a.pang_sort, coeffs05.pang, sin_col, cos_col, prefix = "ASV_") {
  # Subset original sin and cos columns along with Serial column
  newdf_model.pang <- final_2a.pang_sort[, c("Serial", sin_col, cos_col)]
  
  # Automate creation of new columns for each ASV
  for (col_name in colnames(coeffs05.pang)) {
    if (startsWith(col_name, "ASV_")) {
      new_col_name <- paste0(prefix, str_extract(col_name, "\\d+"))
      new_col <- coeffs05.pang[[col_name]][1] + 
        (newdf_model.pang$sinDX1 * coeffs05.pang[[col_name]][2]) + 
        (newdf_model.pang$cosDX1 * coeffs05.pang[[col_name]][3])
      newdf_model.pang[[new_col_name]] <- new_col
    }
  }
  
  return(newdf_model.pang)
}

newdf_model.pang <- automate_harmonic_regression.pang(final_2a.pang_sort, coeffs05.pang, "sinDX1", "cosDX1")
# check newdf_model.pang for number of variables, 142

# subset of newdf_model.pang with only the first column and columns 4th to last
Tomelt.pang <- (newdf_model.pang[c(1,4:142)]) # for p < 0.05, 139 family have seasonal trend

#Add taxa to this melted 
newdf_model.pang_melted <- melt(Tomelt.pang,id.vars="Serial") # 52264 observation

#Add taxa extracting sign_names.pang variable from the Phyloseq object, thereafter extract the tax table only with family
##### Subset from phyloseq 
fam_sign05.pang <- prune_taxa(sign_names.pang, pang.glom)
fam_sign05.pang # 139 taxa and 409 samples remains after p-value adjusted

# Extract Taxonomic Information:
fam_pang.sign <- data.frame(tax_table(fam_sign05.pang))[,5,drop=FALSE] # 5th column (Family)
# Add ASV Information:
fam_pang.sign$ASV <- rownames(fam_pang.sign) #fam_pang.sign contains the necessary to add taxa to melted

#change column name from "variable" to "ASV" in newdf_model.pang_melted 
names(newdf_model.pang_melted)[2] <- "ASV"

#merge by ASV and add the column Genus 
newdf_model.pang_melted_taxa <- merge(newdf_model.pang_melted, fam_pang.sign[, c("Family", "ASV")], by="ASV")

#Subset 1 to 365 (one year)
newdf_model.pang_melted_taxa365 <- newdf_model.pang_melted_taxa[newdf_model.pang_melted_taxa$Serial > 365 & newdf_model.pang_melted_taxa$Serial < 730,]

# consider all the years
newdf_model.pang_melted_taxa365$Serial365 <- (newdf_model.pang_melted_taxa365$Serial-364)

# The code from Luis, was considering sum of values for each ASVs across all samples but was not considering relative abundance. So I used the below code to consider family with mean relative abundance greater than 1%
# Calculate the mean relative abundance for each ASV
mean_abundance_pang <- colMeans(pang_ASVt)

# Filter ASVs with mean relative abundance > 1% (0.01)
h.pang <- colnames(pang_ASVt)[mean_abundance_pang > 0.01] # RA > 1%

# Remove any NAs
h.pang <- h.pang[!is.na(h.pang)]

# Filter the dataset to only include these ASVs
pangfamHigh2 <- newdf_model.pang_melted_taxa[newdf_model.pang_melted_taxa$ASV %in% h.pang, ] # 4136 observation

# check unique family names in the dataframe
unique(pangfamHigh2$Family) # 11 classified, unclassified and candidate family


#### Plot month wise----
# Plotting by day wise makes it difficult to explain the seasonal trend is a bit confusing. So in x axis, plot according to sampling month
# Define the breaks for the first day of each month
month_breaks <- c(275, 306, 336, 367, 
                  398, 426, 457, 487, 518, 548, 579, 610, 640, 671, 701, 732, 
                  763, 791, 822, 852)

# Define the corresponding labels for the breaks
month_labels <- c("Oct-16", "Nov-16", "Dec-16", "Jan-17", "Feb-17", "Mar-17", "Apr-17", 
                  "May-17", "Jun-17", "Jul-17", "Aug-17", "Sep-17", "Oct-17", "Nov-17", "Dec-17", 
                  "Jan-18", "Feb-18", "Mar-18", "Apr-18", "May-18")

# Add a new column "Sampling_month" based on the Serial column
pangfamHigh2 <- pangfamHigh2 %>%
  mutate(Sampling_month = case_when(
    Serial >= 275 & Serial < 306  ~ "Oct-16",
    Serial >= 306 & Serial < 336  ~ "Nov-16",
    Serial >= 336 & Serial < 367  ~ "Dec-16",
    Serial >= 367 & Serial < 398  ~ "Jan-17",
    Serial >= 398 & Serial < 426  ~ "Feb-17",
    Serial >= 426 & Serial < 457  ~ "Mar-17",
    Serial >= 457 & Serial < 487  ~ "Apr-17",
    Serial >= 487 & Serial < 518  ~ "May-17",
    Serial >= 518 & Serial < 548  ~ "Jun-17",
    Serial >= 548 & Serial < 579  ~ "Jul-17",
    Serial >= 579 & Serial < 610  ~ "Aug-17",
    Serial >= 610 & Serial < 640  ~ "Sep-17",
    Serial >= 640 & Serial < 671  ~ "Oct-17",
    Serial >= 671 & Serial < 701  ~ "Nov-17",
    Serial >= 701 & Serial < 732  ~ "Dec-17",
    Serial >= 732 & Serial < 763  ~ "Jan-18",
    Serial >= 763 & Serial < 791  ~ "Feb-18",
    Serial >= 791 & Serial < 822  ~ "Mar-18",
    Serial >= 822 & Serial < 852  ~ "Apr-18",
    Serial >= 852 & Serial < 883  ~ "May-18",
    TRUE ~ NA_character_  # Handle cases outside the defined ranges
  ))
# View the updated dataframe
head(pangfamHigh2)

# Plot using custom x-axis breaks and labels
pang.fam.high2 <- ggplot(pangfamHigh2, aes(x = Serial, y = value, colour = Family, 
                                           group = Family, label = Family)) + 
  geom_line(linewidth = 1) +  
  #directlabels::geom_dl(aes(label = Family),  method = list(cex = 0.9, "smart.grid")) + 
  scale_color_manual(values = family_colors) + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, margin = margin(b = 15)), 
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 18, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 18, margin = margin(r = 5)),
    axis.title.x = element_blank(),#text(size = 18, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18, face = "italic"),
    legend.spacing.y = unit(1.3, "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black")) + 
  labs(x = "Sampling months", y = "Relative abundance", 
       title = "Prokaryotes: Pangasius ponds") + 
  scale_x_continuous(breaks = month_breaks, labels = month_labels)
pang.fam.high2

##### Interactive plot to see which line is which family----
library(plotly); library("plotly")
ggplotly(pang.fam.high2)

### Tilapia----
# Subset tilapia
pseq_tila <- ps %>%
  ps_filter(Sampling_point > 1,
            Crop == "Tilapia",
            Crop_species != "Gulsha.Carp",
            Crop_species != "Gulsha.Pabda") %>% 
  tax_fix()
pseq_tila #  8338 taxa and 385 samples

##### convert reads to %
tila.rel <- transform_sample_counts(pseq_tila, function(x){x / sum(x)})

# Aggregate at Family rank
#tila.glom <- tax_glom(tila.rel, taxrank = "Family") # Save it and use later as it takes long time to process
#saveRDS(tila.glom, "pseq_tila_aggregated_at_family.rds") 
tila.glom <- readRDS("pseq_tila_aggregated_at_family.rds")
tila.glom # 566 taxa and 385 samples

# Extract the taxonomic information
tila.tax<-as.data.frame(tax_table(tila.glom)[,5]) # column 5 (Family)
tila.tax$ASV <- rownames(tila.tax)

# Extract the ASV table
tila_ASV <- as.data.frame(otu_table(tila.glom))
# Transpose the ASV table
tila_ASVt <- t(tila_ASV) # transpose

#####Create a data.frame cols=ASV, rows= samples -> add date to this and estimate the sin/cos thing 
# Extract Family Taxonomy:
tila.tax <- as.data.frame(tax_table(tila.glom)[,5]) ###Keep this one to add to the melted object
# Convert ASV Table to Dataframe:
tila_ASV <- as.data.frame(otu_table(tila.glom))
# Transpose the ASV Table:
ASV_tila_trans <- as.data.frame(x = t(tila_ASV), stringsAsFactors = FALSE)
# Extract Metadata (Sampling Date):
md_to_add_tila <- as.data.frame(sample_data(tila.glom))[,c(12)] #culumn number of Sampling_date in metadata
# Combine ASV Data and Metadata:
final_2a_tila <- cbind(ASV_tila_trans,md_to_add_tila)
# Convert Sampling Dates to Date Format:
final_2a_tila$Sampling_date <- as.Date(final_2a_tila$Sampling_date,"%d/%m/%Y") # should i change as %Y/%m/%d
# Sort Data by Sampling Date:
final_2a_tila_sort <- final_2a_tila[order(final_2a_tila$Sampling_date),] ##sort from least recent to most recent
# Convert the Date Format Again:
final_2a_tila_sort$Sampling_date <- as.Date(final_2a_tila_sort$Sampling_date,"%Y-%m-%d")

##Add the Julian day 
final_2a_tila_sort$Jul <- format(final_2a_tila_sort$Sampling_date, "%j")

###Serial Day
#Starting from
tila_start <- as.Date('2016-01-01') # should be first day of the sampling year

final_2a_tila_sort$Serial <- as.numeric(-difftime(tila_start,final_2a_tila_sort$Sampling_date)+1)

##Days since winter solstice (21 December)
##Add the winter solstice 
final_2a_tila_sort$wintersols[final_2a_tila_sort$Sampling_date>as.Date("2016-01-01") & final_2a_tila_sort$Sampling_date<as.Date("2016-12-21")] <- "2015-12-21"
final_2a_tila_sort$wintersols[final_2a_tila_sort$Sampling_date>as.Date("2017-01-01") & final_2a_tila_sort$Sampling_date<as.Date("2017-12-21")] <- "2016-12-21"
final_2a_tila_sort$wintersols[final_2a_tila_sort$Sampling_date>as.Date("2018-01-01") & final_2a_tila_sort$Sampling_date<as.Date("2018-12-21")] <- "2017-12-21"

final_2a_tila_sort$DaysSincewintersols <- as.numeric(-difftime(final_2a_tila_sort$wintersols,final_2a_tila_sort$Sampling_date) + 1) #Add the numeric column

# Remove Rows with Missing Values:
final_2a_tila_sort <- final_2a_tila_sort[!is.na(final_2a_tila_sort$DaysSincewintersols), ]

## Convert Julian date to numeric:
final_2a_tila_sort$Jul <- as.numeric(final_2a_tila_sort$Jul)

# Generate sine and cosine terms for harmonic regression:
final_2a_tila_sort$cosDX1 <- as.numeric(cos(2*pi*((final_2a_tila_sort$DaysSincewintersols)/365))) #Peaking in midwinter
final_2a_tila_sort$sinDX1 <- as.numeric(sin(2*pi*((final_2a_tila_sort$Jul)/365))) #Peaking other time
# check final_2a_tila_sort for number of variables which will be needed below. remove last 7 column

# Fit harmonic regression models:
storage.tila <- list()
for(i in names(final_2a_tila_sort)[1:566]){ #number of taxa present in final_2a_tila_sort (573) - 7
  storage.tila[[i]] <- lm(get(i) ~ sinDX1+cosDX1, final_2a_tila_sort)
}

# Extract p-values for significance testing:
pvalues_lm.tila <- sapply(storage.tila, function(x){
  ff <- summary(x)$fstatistic
  pf(ff[1], df1 = ff[2], df2 = ff[3], lower.tail = FALSE)
})

# Adjust the p-values using the BH method
pvalues_lm.tila_adjusted <- p.adjust(pvalues_lm.tila, method = "BH")

# Filter significant taxa based on adjusted p-values
pvalues_lm.tila_adjusted05 <- pvalues_lm.tila_adjusted[pvalues_lm.tila_adjusted < 0.05] # 156

# Extract significant taxa names after p-value adjustment
sign_names.tila <- str_remove(names(pvalues_lm.tila_adjusted05), ".value")

#Number of the column that matches a list of names, get the position of the names <0.05
signll.tila <- match(sign_names.tila,names(final_2a_tila_sort))

#subset list of lists 
storage.tila01 <- storage.tila[signll.tila]

#Get the constant, sin and cos of the models:
# check if storage.tila05 contains any NAs and if so, Filter out NA values from storage.tila05
filtered_storage.tila <- storage.tila01[!sapply(storage.tila01, is.null)]

# Apply sapply to filtered storage.tila
coeffs05.tila <- data.frame(coeffs = sapply(filtered_storage.tila, function(item) item$coefficients))
colnames(coeffs05.tila) <- str_remove(colnames(coeffs05.tila), "coeffs.")

#Subset original day, sin, cos
newdf_model.tila <- final_2a_tila_sort[, c("Serial","cosDX1","sinDX1")]

# Automatize this :P 
# use the function below to automatically do it, and add the value in newdf_model.tila table
automate_harmonic_regression.tila <- function(final_2a_tila_sort, coeffs05.tila, sin_col, cos_col, prefix = "ASV_") {
  # Subset original sin and cos columns along with Serial column
  newdf_model.tila <- final_2a_tila_sort[, c("Serial", sin_col, cos_col)]
  
  # Automate creation of new columns for each ASV
  for (col_name in colnames(coeffs05.tila)) {
    if (startsWith(col_name, "ASV_")) {
      new_col_name <- paste0(prefix, str_extract(col_name, "\\d+"))
      new_col <- coeffs05.tila[[col_name]][1] + 
        (newdf_model.tila$sinDX1 * coeffs05.tila[[col_name]][2]) + 
        (newdf_model.tila$cosDX1 * coeffs05.tila[[col_name]][3])
      newdf_model.tila[[new_col_name]] <- new_col
    }
  }
  
  return(newdf_model.tila)
}

newdf_model.tila <- automate_harmonic_regression.tila(final_2a_tila_sort, coeffs05.tila, "sinDX1", "cosDX1")
# check newdf_model.tila for number of variables 159 with p < 0.05

# subset of newdf_model.tila with only the first column and columns 4th to last
Tomelt.tila <- (newdf_model.tila[c(1,4:159)]) # newdf_model.tila - 3 = 156 families were significant for p < 0.05

#Add taxa to this melted 
newdf_model.tila_melted <- melt(Tomelt.tila,id.vars="Serial")

#Add taxa extracting sign_names.tila variable from the Phyloseq object, thereafter extract the tax table only with family
##### Subset from phyloseq 
fam_sign05.tila <- prune_taxa(sign_names.tila, tila.glom)
fam_sign05.tila # 156 taxa and 385 samples with significant trend

# Extract Taxonomic Information:
tila_fam.sign <- data.frame(tax_table(fam_sign05.tila))[,5,drop=FALSE] # 5th column (Family)
# Add ASV Information:
tila_fam.sign$ASV <- rownames(tila_fam.sign) #tila_fam.sign contains the neccesary to add taxa to melted

#change column name from "variable" to "ASV" in newdf_model.tila_melted 
names(newdf_model.tila_melted)[2] <- "ASV"

#merge by ASV and add the column phylum
newdf_model.tila_melted_taxa <- merge(newdf_model.tila_melted, tila_fam.sign[, c("Family", "ASV")], by="ASV")
# 56316 observation

#Subset 1 to 365 
newdf_model.tila_melted_taxa365 <- newdf_model.tila_melted_taxa[newdf_model.tila_melted_taxa$Serial > 365 & newdf_model.tila_melted_taxa$Serial < 730,]

# Condiser all the years
newdf_model.tila_melted_taxa365$Serial365 <- (newdf_model.tila_melted_taxa365$Serial-364) # 39780 observation with p < 0.05

# Consider only those with a mean relative abundance >1% across all samples
# Calculate the mean relative abundance for each ASV
mean_abundance_tila <- colMeans(tila_ASVt)

# Filter ASVs with mean relative abundance > 1% (0.01)
h.tila <- colnames(tila_ASVt)[mean_abundance_tila > 0.01]

# Remove any NAs
h.tila <- h.tila[!is.na(h.tila)]

# Filter the dataset to only include these ASVs
tilafamHigh2 <- newdf_model.tila_melted_taxa[newdf_model.tila_melted_taxa$ASV %in% h.tila, ] # 5776 observation

# check unique family names in the dataframe
unique(tilafamHigh2$Family) # 16 families


#### Plot month wise----
# Add a new column "Sampling_month" based on the Serial column
tilafamHigh2 <- tilafamHigh2 %>%
  mutate(Sampling_month = case_when(
    Serial >= 275 & Serial < 306  ~ "Oct-16",
    Serial >= 306 & Serial < 336  ~ "Nov-16",
    Serial >= 336 & Serial < 367  ~ "Dec-16",
    Serial >= 367 & Serial < 398  ~ "Jan-17",
    Serial >= 398 & Serial < 426  ~ "Feb-17",
    Serial >= 426 & Serial < 457  ~ "Mar-17",
    Serial >= 457 & Serial < 487  ~ "Apr-17",
    Serial >= 487 & Serial < 518  ~ "May-17",
    Serial >= 518 & Serial < 548  ~ "Jun-17",
    Serial >= 548 & Serial < 579  ~ "Jul-17",
    Serial >= 579 & Serial < 610  ~ "Aug-17",
    Serial >= 610 & Serial < 640  ~ "Sep-17",
    Serial >= 640 & Serial < 671  ~ "Oct-17",
    Serial >= 671 & Serial < 701  ~ "Nov-17",
    Serial >= 701 & Serial < 732  ~ "Dec-17",
    Serial >= 732 & Serial < 763  ~ "Jan-18",
    Serial >= 763 & Serial < 791  ~ "Feb-18",
    Serial >= 791 & Serial < 822  ~ "Mar-18",
    Serial >= 822 & Serial < 852  ~ "Apr-18",
    Serial >= 852 & Serial < 883  ~ "May-18",
    TRUE ~ NA_character_  # Handle cases outside the defined ranges
  ))
# View the updated dataframe
head(tilafamHigh2)

# Plot using custom x-axis breaks and labels
tila.fam.high2 <- ggplot(tilafamHigh2, aes(x = Serial, y = value, colour = Family, group = Family, label = Family)) + 
  geom_line(linewidth = 1) +
  #directlabels::geom_dl(aes(label = Family),  method = list(cex = 0.9, "smart.grid")) + 
  scale_color_manual(values = family_colors) + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, margin = margin(b = 15)), 
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 18, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 18, margin = margin(r = 5)),
    axis.title.x = element_blank(),#text(size = 18, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18, face = "italic"),
    legend.spacing.y = unit(1.3, "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black")) + 
  labs(x = "day", y = "Relative abundance", 
       title = "Prokaryotes: Tilapia ponds") + 
  scale_x_continuous(breaks = month_breaks, labels = month_labels)
tila.fam.high2

##### Interactive plot----
ggplotly(tila.fam.high2)


# Microeukaryotes----
# This script is for bacterial seasonal dynamics including shannon diversity, glm and harmonic regression at family rank
# Alpha diversity----
# Load phyloseq object
ps.18s <- readRDS("phyloseq_18S_filtered_with_tree_pr2_90-150bp_20240416.rds")
ps.18s # 5390 taxa and 872 samples

# Manually define the correct order of months
month_order <- c("Jun-16","Oct-16","Nov-16","Dec-16","Jan-17","Feb-17","Mar-17",
                 "Apr-17","May-17","Jun-17","Jul-17","Aug-17","Sep-17","Oct-17", 
                 "Nov-17","Dec-17","Jan-18","Feb-18","Mar-18","Apr-18","May-18")

# Access the sample data from the phyloseq object and modify the Sampling_months column
sample_data(ps.18s)$Sampling_months <- factor(sample_data(ps.18s)$Sampling_months, levels = month_order)

# Pangasius ponds only----
pseq_pang.18s <- ps.18s %>% 
  ps_filter(
    Crop == "Pangasius" & Crop_species != "Shing") %>% 
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
pseq_pang.18s # 3947 taxa and 421 samples

# Access the sample data from the phyloseq object and modify the Sampling_months column
sample_data(pseq_pang.18s)$Season2 <- factor(sample_data(pseq_pang.18s)$Season2, levels = season_order2)
sample_data(pseq_pang.18s)$Season3 <- factor(sample_data(pseq_pang.18s)$Season3, levels = season_order3)

## A3. Alpha diversity----
#set.seed(1234)#The set. seed() function in R is used to create reproducible results when writing code that involves creating variables that take on random values. By using the set. seed() function, you guarantee that the same random values are produced each time you run the code
ps_rarefy_2.1.18s <- rarefy_even_depth(pseq_pang.18s, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#285OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_2.1.18s <- estimate_richness(ps_rarefy_2.1.18s, measures = c("Chao1", "Shannon"))
alpha_estimates_2.1.18s <- cbind(alpha_estimates_2.1.18s, sample_data(pseq_pang.18s))

### Shannon----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
p.sm.s.18s <- ggplot(data = alpha_estimates_2.1.18s, aes(y = Shannon, x = Sampling_months)) +
  geom_boxplot(aes(color = Season2)) +
  geom_jitter(aes(color = Season2), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = season.colors) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 18, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 18, margin = margin(r = 5)),
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
  guides(fill = guide_legend(nrow = 1)) + # Set number of columns in the legend
  labs (color = "Season",  # Customize legend title
        y = "Shannon", # Add y-axis title to the plot
        x = "Sampling months",
        title = "Microeukaryotes: Pangasius ponds") +
  ylim(0, 6)
p.sm.s.18s

### statistical test----
#### Check Normality of Shannon Diversity----
# Histogram
hist(alpha_estimates_2.1.18s$Shannon, main="Histogram of Shannon Diversity", xlab="Shannon Index", col="lightblue", breaks=20)

# Q-Q plot
qqnorm(alpha_estimates_2.1.18s$Shannon)
qqline(alpha_estimates_2.1.18s$Shannon, col = "red")

# Shapiro-Wilk test
shapiro.test(alpha_estimates_2.1.18s$Shannon)
# p-value < 2.131e-09 # not normally distributed

# Skewness and Kurtosis
skewness(alpha_estimates_2.1.18s$Shannon)
# -0.7737533 # The skewness value of -0.7737533 suggests that the distribution is negatively skewed (a tail on the left). which we can see from the histogram

kurtosis(alpha_estimates_2.1.18s$Shannon)
# 3.46102 # The kurtosis value of 3.46102 is slightly higher than 3 (kurtosis value), indicating that the distribution has heavy tails and is more peaked compared to a normal distribution.
# The Shannon index is not normally distributed. Since the normality assumption is violated, a standard Gaussian GLM may not be the best fit.

#### Gamma Distribution GLM (with log link)----
#Since the Shannon diversity is positive and skewed, a Gamma distribution with a log link might be more appropriate. This assumes that the Shannon index follows a gamma distribution, which handles positive, skewed data.
glm_gamma_pang.18s <- glm(Shannon ~ Season3, family = Gamma(link = "log"), data = alpha_estimates_2.1.18s)
glm_pang_summary.18s <- summary(glm_gamma_pang.18s)

# Extract coefficients
glm_coeff_pang.18s <- as.data.frame(glm_pang_summary.18s$coefficients)

# Add the row names (e.g., '(Intercept)', 'Season02.Winter', etc.) as a new column called 'Variable'
glm_coeff_pang.18s$Variable <- rownames(glm_coeff_pang.18s)

# Move 'Variable' to the first position
glm_coeff_pang.18s <- glm_coeff_pang.18s[, c("Variable", colnames(glm_coeff_pang.18s)[1:4])]

# Format the p-values to avoid rounding off in Excel
#glm_coeff_pang$`Pr(>|t|)` <- format(glm_coeff_pang$`Pr(>|t|)`, scientific = TRUE)

#### Pairwise comparison----
pairwise_glm_pang.18s <- emmeans(glm_gamma_pang.18s, pairwise ~ Season3, adjust = "BH")
summary(pairwise_glm_pang.18s)

# Extract pairwise comparison results
emmeans_table_pang.18s <- as.data.frame(pairwise_glm_pang.18s$emmeans)
contrasts_table_pang.18s <- as.data.frame(pairwise_glm_pang.18s$contrasts)

# Optional: You can also check diagnostic plots for the GLM
par(mfrow = c(2, 2))  # Plot 4 diagnostic plots
plot(glm_gamma_pang.18s)

# Create a new Excel workbook
wb_pang.18s <- createWorkbook()

# Add sheets for GLM summary, emmeans, and contrasts
addWorksheet(wb_pang.18s, "GLM Summary")
addWorksheet(wb_pang.18s, "EMMeans")
addWorksheet(wb_pang.18s, "Contrasts")

# Write data to the respective sheets
writeData(wb_pang.18s, sheet = "GLM Summary", glm_coeff_pang.18s)
writeData(wb_pang.18s, sheet = "EMMeans", emmeans_table_pang.18s)
writeData(wb_pang.18s, sheet = "Contrasts", contrasts_table_pang.18s)
# Save the workbook to a file
#saveWorkbook(wb_pang.18s, file = "glm_and_pairwise_results_pseq_pang_18s.xlsx", overwrite = TRUE)

# Tilapia ponds only----
pseq_tila.18s <- ps.18s %>% 
  ps_filter(Crop == "Tilapia",
            Crop_species != "Gulsha.Carp" & Crop_species != "Gulsha.Pabda") %>% 
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
pseq_tila.18s #  8380 taxa and 391 samples

# Access the sample data from the phyloseq object and modify the Sampling_months column
sample_data(pseq_tila.18s)$Season2 <- factor(sample_data(pseq_tila.18s)$Season2, levels = season_order2)
sample_data(pseq_tila.18s)$Season3 <- factor(sample_data(pseq_tila.18s)$Season3, levels = season_order3)

## A4. Alpha diversity----
ps_rarefy_2.2.18s <- rarefy_even_depth(pseq_tila.18s, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#287OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_2.2.18s <- estimate_richness(ps_rarefy_2.2.18s, measures = c("Chao1", "Shannon"))
alpha_estimates_2.2.18s <- cbind(alpha_estimates_2.2.18s, sample_data(pseq_tila.18s))

### Shannon----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
t.sm.s.18s <- ggplot(data = alpha_estimates_2.2.18s, aes(y = Shannon, x = Sampling_months)) +
  geom_boxplot(aes(color = Season2)) +
  geom_jitter(aes(color = Season2), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  #geom_point() +
  scale_color_manual(values = season.colors) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 18, margin = margin(t = 5)),
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
  guides(fill = guide_legend(nrow = 1)) + # Set number of columns in the legend
  labs (color = "Season",  # Customize legend title
        y = "Shannon", # Add y-axis title to the plot
        x = "Sampling months",
        title = "Microeukaryotes: Tilapia ponds") +
  ylim(0, 6)
t.sm.s.18s

#### Check Normality of Shannon Diversity----
# Histogram
hist(alpha_estimates_2.2.18s$Shannon, main="Histogram of Shannon Diversity", xlab="Shannon Index", col="lightblue", breaks=20)

# Q-Q plot
qqnorm(alpha_estimates_2.2.18s$Shannon)
qqline(alpha_estimates_2.2.18s$Shannon, col = "red")

# Shapiro-Wilk test
shapiro.test(alpha_estimates_2.2.18s$Shannon)
# p-value =  7.039e-12 # not normally distributed

# Skewness and Kurtosis
skewness(alpha_estimates_2.2.18s$Shannon)
# -1.073949 # The skewness value of --1.073949 suggests that the distribution is negatively highly left-skewed. which we can see from the histogram

kurtosis(alpha_estimates_2.2.18s$Shannon)
#  4.421637 # The kurtosis value of  4.421637 is higher than 3 (kurtosis value), indicating that the distribution has a sharper peak and heavier tails than a normal distribution.
# The Shannon index is not normally distributed. Since the normality assumption is violated, a standard Gaussian GLM may not be the best fit.

#### Gamma Distribution GLM (with log link)----
#Since the Shannon diversity is positive and skewed, a Gamma distribution with a log link might be more appropriate. This assumes that the Shannon index follows a gamma distribution, which handles positive, skewed data.
glm_gamma_tila.18s <- glm(Shannon ~ Season3, family = Gamma(link = "log"), data = alpha_estimates_2.2.18s)
glm_tila_summary.18s <- summary(glm_gamma_tila.18s)

# Extract coefficients
glm_coeff_tila.18s <- as.data.frame(glm_tila_summary.18s$coefficients)

# Add the row names (e.g., '(Intercept)', 'Season02.Winter', etc.) as a new column called 'Variable'
glm_coeff_tila.18s$Variable <- rownames(glm_coeff_tila.18s)

# Move 'Variable' to the first position
glm_coeff_tila.18s <- glm_coeff_tila.18s[, c("Variable", colnames(glm_coeff_tila.18s)[1:4])]

# Format the p-values to avoid rounding off in Excel
#glm_coeff_tila.18s$`Pr(>|t|)` <- format(glm_coeff_tila.18s$`Pr(>|t|)`, scientific = TRUE)

#### Pairwise comparison----
pairwise_glm_tila.18s <- emmeans(glm_gamma_tila.18s, pairwise ~ Season3, adjust = "BH")
summary(pairwise_glm_tila.18s)

# Extract pairwise comparison results
emmeans_table_tila.18s <- as.data.frame(pairwise_glm_tila.18s$emmeans)
contrasts_table_tila.18s <- as.data.frame(pairwise_glm_tila.18s$contrasts)

# Optional: You can also check diagnostic plots for the GLM
par(mfrow = c(2, 2))  # Plot 4 diagnostic plots
plot(glm_gamma_tila.18s)

# Create a new Excel workbook
wb_tila.18s <- createWorkbook()

# Add sheets for GLM summary, emmeans, and contrasts
addWorksheet(wb_tila.18s, "GLM Summary")
addWorksheet(wb_tila.18s, "EMMeans")
addWorksheet(wb_tila.18s, "Contrasts")

# Write data to the respective sheets
writeData(wb_tila.18s, sheet = "GLM Summary", glm_coeff_tila.18s)
writeData(wb_tila.18s, sheet = "EMMeans", emmeans_table_tila.18s)
writeData(wb_tila.18s, sheet = "Contrasts", contrasts_table_tila.18s)

# Save the workbook to a file
#saveWorkbook(wb_tila.18s, file = "glm_and_pairwise_results_pseq_tila.18s.xlsx", overwrite = TRUE)

# Harmonic regression----
# Define colours for family (these are mostly save as the barplot but only for a few families, i changed the color as they were very faded in this plot)
family_colors.18s <- c("Alveolata Division" = "blue","Aulacoseiraceae" = "#cfd251",
                       "Cryptomonadales_X" = "#1ff8ff", "Gyrista Subdivision" = "#da6130",
                       "Hemiselmidaceae" = "#A6761D","Kathablepharidida_X" = "#69c86c",
                       "Oxytrichidae" = "#6A3D9A","Peronosporales" = "#532a5a",
                       "Raphidophyceae_XX" = "#E6AB02", "Spirotrichea Class" = "#E7298A",
                       "Stephanodiscaceae" =  "#66A61E", "Thalassiosirales Order" = "#FB9A99", 
                       "Cryptophyceae Class" = "#c29545","Dinophyceae Class" = "#FF7F00",
                       "Saccharomycetales" = "#83d5af",
                       "Tintinnidiidae" = "#A6CEE3","TSAR Supergroup" = "#D95F02")

### Pangasius ponds----
# Subset healthy pangasius (remove algal bloom) & also remove ponds with no fish
pseq_pang.18s <- ps.18s %>%
  ps_filter(Sampling_point > 1,
            Crop == "Pangasius" & Crop_species != "Shing") %>% 
  tax_fix()
pseq_pang.18s # 3936 taxa and 409 samples

##### convert reads to relative abundance
pang.rel.18s <- transform_sample_counts(pseq_pang.18s, function(x){x / sum(x)})

# Aggregate at Family rank
#pang.glom.18s <- tax_glom(pang.rel.18s, taxrank = "Family")
#saveRDS(pang.glom.18s, "pseq_pang_18s_aggregated_at_family.rds")
pang.glom.18s <- readRDS("pseq_pang_18s_aggregated_at_family.rds")
pang.glom.18s # 292 taxa and 409 samples

# Extract the taxonomic information
pang.tax.18s <-as.data.frame(tax_table(pang.glom.18s)[,7]) # 7th column (Family)
pang.tax.18s$ASV <- rownames(pang.tax.18s)

# Extract the ASV table
pang_ASV.18s <- as.data.frame(otu_table(pang.glom.18s))
# Transpose the ASV table
pang_ASVt.18s <- t(pang_ASV.18s) # transpose

#####Create a data.frame cols=ASV, rows= samples -> add date to this and estimate the sin/cos thing 
pang.tax.18s <- as.data.frame(tax_table(pang.glom.18s)[,7]) ###Keep this one to add to the melted object

pang_ASV.18s <- as.data.frame(otu_table(pang.glom.18s))

ASV_pang_trans.18s <- as.data.frame(x = t(pang_ASV.18s), stringsAsFactors = FALSE)
md_to_add.pang.18s <- as.data.frame(sample_data(pang.glom.18s))[,c(12)] #column number for Sampling_date in metadata

final_2a.pang.18s <- cbind(ASV_pang_trans.18s, md_to_add.pang.18s)
final_2a.pang.18s$Sampling_date <- as.Date(final_2a.pang.18s$Sampling_date,"%d/%m/%Y") # should i change as %Y/%m/%d

final_2a.pang_sort.18s <- final_2a.pang.18s[order(final_2a.pang.18s$Sampling_date),] ##sort from least recent to most recent
final_2a.pang_sort.18s$Sampling_date <- as.Date(final_2a.pang_sort.18s$Sampling_date,"%Y-%m-%d")

##Add the Julian day # what is this
final_2a.pang_sort.18s$Jul <- format(final_2a.pang_sort.18s$Sampling_date, "%j")

###Serial Day
#Starting from
pang_start.18s <- as.Date('2016-01-01') # should be first day of the sampling year

final_2a.pang_sort.18s$Serial <- as.numeric(-difftime(pang_start.18s,final_2a.pang_sort.18s$Sampling_date)+1)

##Days since winter solstice (21 December)
##Add the winter solstice 
final_2a.pang_sort.18s$wintersols[final_2a.pang_sort.18s$Sampling_date>as.Date("2016-01-01") & final_2a.pang_sort.18s$Sampling_date<as.Date("2016-12-21")] <- "2015-12-21"
final_2a.pang_sort.18s$wintersols[final_2a.pang_sort.18s$Sampling_date>as.Date("2017-01-01") & final_2a.pang_sort.18s$Sampling_date<as.Date("2017-12-21")] <- "2016-12-21"
final_2a.pang_sort.18s$wintersols[final_2a.pang_sort.18s$Sampling_date>as.Date("2018-01-01") & final_2a.pang_sort.18s$Sampling_date<as.Date("2018-12-21")] <- "2017-12-21"

final_2a.pang_sort.18s$DaysSincewintersols <- as.numeric(-difftime(final_2a.pang_sort.18s$wintersols,final_2a.pang_sort.18s$Sampling_date) + 1) #Add the numeric column

# Remove Rows with Missing Values:
final_2a.pang_sort.18s <- final_2a.pang_sort.18s[!is.na(final_2a.pang_sort.18s$DaysSincewintersols), ]

## Convert Julian date to numeric:
final_2a.pang_sort.18s$Jul <- as.numeric(final_2a.pang_sort.18s$Jul)

# Generate sine and cosine terms for harmonic regression:
final_2a.pang_sort.18s$cosDX1 <- as.numeric(cos(2*pi*((final_2a.pang_sort.18s$DaysSincewintersols)/365))) #Peaking in midwinter
final_2a.pang_sort.18s$sinDX1 <- as.numeric(sin(2*pi*((final_2a.pang_sort.18s$Jul)/365))) #Peaking other time
# check final_2a.pang_sort.18s for number of variables which will be needed below. Remove last 7 columns

# Fit harmonic regression models:
storage.pang.18s <- list()
for(i in names(final_2a.pang_sort.18s)[1:292]){ #number of taxa present in final_2a.pang_sort.18s (299) -7
  storage.pang.18s[[i]] <- lm(get(i) ~ sinDX1+cosDX1, final_2a.pang_sort.18s)
}

# Extract p-values for significance testing:
pvalues_lm.pang.18s <- sapply(storage.pang.18s, function(x){
  ff <- summary(x)$fstatistic
  pf(ff[1], df1 = ff[2], df2 = ff[3], lower.tail = FALSE)
})

# Adjust the p-values using the BH method
pvalues_lm.pang_adjusted.18s <- p.adjust(pvalues_lm.pang.18s, method = "BH")

# Filter significant taxa based on adjusted p-values
pvalues_lm.pang_adjusted05.18s <- pvalues_lm.pang_adjusted.18s[pvalues_lm.pang_adjusted.18s < 0.05] # 78

# Extract significant taxa names after p-value adjustment
sign_names.pang.18s <- str_remove(names(pvalues_lm.pang_adjusted05.18s), ".value")

#Number of the column that matches a list of names, get the position of the names <0.05
signll.pang.18s <- match(sign_names.pang.18s, names(final_2a.pang_sort.18s))

#subset list of lists 
storage05.pang.18s <- storage.pang.18s[signll.pang.18s]

#Get the constant, sin and cos of the models:
# check if storage05.pang contains any NAs and if so, Filter out NA values from storage05.pang
filtered_storage.pang.18s <- storage05.pang.18s[!sapply(storage05.pang.18s, is.null)]

# Apply sapply to filtered storage
coeffs05.pang.18s <- data.frame(coeffs = sapply(filtered_storage.pang.18s, function(item) item$coefficients))
colnames(coeffs05.pang.18s) <- str_remove(colnames(coeffs05.pang.18s), "coeffs.")

#Subset original day, sin, cos
newdf_model.pang.18s <- final_2a.pang_sort.18s[, c("Serial","cosDX1","sinDX1")]

# Automatize this :P 
# use the function below to automatically do it, and add the value in newdf_model.pang table
automate_harmonic_regression.pang.18s <- function(final_2a.pang_sort.18s, 
                                                  coeffs05.pang.18s, sin_col, cos_col, prefix = "ASV_") {
  # Subset original sin and cos columns along with Serial column
  newdf_model.pang.18s <- final_2a.pang_sort.18s[, c("Serial", sin_col, cos_col)]
  
  # Automate creation of new columns for each ASV
  for (col_name in colnames(coeffs05.pang.18s)) {
    if (startsWith(col_name, "ASV_")) {
      new_col_name <- paste0(prefix, str_extract(col_name, "\\d+"))
      new_col <- coeffs05.pang.18s[[col_name]][1] + 
        (newdf_model.pang.18s$sinDX1 * coeffs05.pang.18s[[col_name]][2]) + 
        (newdf_model.pang.18s$cosDX1 * coeffs05.pang.18s[[col_name]][3])
      newdf_model.pang.18s[[new_col_name]] <- new_col
    }
  }
  
  return(newdf_model.pang.18s)
}

newdf_model.pang.18s <- automate_harmonic_regression.pang.18s(final_2a.pang_sort.18s, coeffs05.pang.18s, "sinDX1", "cosDX1")
# check newdf_model.pang.18s for number of variables
Tomelt.pang.18s <- (newdf_model.pang.18s[c(1,4:78)]) # for p < 0.05, 75 families

#Add taxa to this melted 
newdf_model.pang_melted.18s <- melt(Tomelt.pang.18s, id.vars="Serial") # 28350 observation

#Add taxa extracting sign_names.pang variable from the Phyloseq object, thereafter extract the tax table only with genus 
##### Subset from physeq 
fam_sign05.pang.18s <- prune_taxa(sign_names.pang.18s, pang.glom.18s)
fam_sign05.pang.18s # 75 taxa and 409 samples with significant seasonal trend

fam_pang.sign.18s <- data.frame(tax_table(fam_sign05.pang.18s))[,7,drop=FALSE] # 7th column (Family)

fam_pang.sign.18s$ASV <- rownames(fam_pang.sign.18s) #fam_pang.sign.18s contains the neccesary to add taxa to melted

#change column name from "variable" to "ASV" in newdf_model.pang_melted 
names(newdf_model.pang_melted.18s)[2] <- "ASV"

#merge by ASV and add the column Genus 
newdf_model.pang_melted_taxa.18s <- merge(newdf_model.pang_melted.18s, 
                                          fam_pang.sign.18s[, c("Family", "ASV")], by="ASV")

#Subset 1 to 365 (one year)
newdf_model.pang_melted_taxa365.18s <- newdf_model.pang_melted_taxa.18s[newdf_model.pang_melted_taxa.18s$Serial > 365 & newdf_model.pang_melted_taxa.18s$Serial < 730,]

# I will plot all the years i have
newdf_model.pang_melted_taxa365.18s$Serial365 <- (newdf_model.pang_melted_taxa365.18s$Serial-364)

# Keep only those with a relative abundance > 1%
# Calculate the mean relative abundance for each ASV
mean_abundance_pang_18s <- colMeans(pang_ASVt.18s)

# Filter ASVs with mean relative abundance > 1% (0.01)
h.pang.18s <- colnames(pang_ASVt.18s)[mean_abundance_pang_18s > 0.01]

# Remove any NAs
h.pang.18s <- h.pang.18s[!is.na(h.pang.18s)]

### i want to plot all years instead of just one year
# instead of just one year, if i want to plot all the sampling months
pangfamHigh2.18s <- newdf_model.pang_melted_taxa.18s[newdf_model.pang_melted_taxa.18s$ASV %in% h.pang.18s, ] # 3780 observation

# check unique family names in the dataframe
unique(pangfamHigh2.18s$Family) # 12 families


#### Plot month wise----
# Add a new column "Sampling_month" based on the Serial column
pangfamHigh2.18s <- pangfamHigh2.18s %>%
  mutate(Sampling_month = case_when(
    Serial >= 275 & Serial < 306  ~ "Oct-16",
    Serial >= 306 & Serial < 336  ~ "Nov-16",
    Serial >= 336 & Serial < 367  ~ "Dec-16",
    Serial >= 367 & Serial < 398  ~ "Jan-17",
    Serial >= 398 & Serial < 426  ~ "Feb-17",
    Serial >= 426 & Serial < 457  ~ "Mar-17",
    Serial >= 457 & Serial < 487  ~ "Apr-17",
    Serial >= 487 & Serial < 518  ~ "May-17",
    Serial >= 518 & Serial < 548  ~ "Jun-17",
    Serial >= 548 & Serial < 579  ~ "Jul-17",
    Serial >= 579 & Serial < 610  ~ "Aug-17",
    Serial >= 610 & Serial < 640  ~ "Sep-17",
    Serial >= 640 & Serial < 671  ~ "Oct-17",
    Serial >= 671 & Serial < 701  ~ "Nov-17",
    Serial >= 701 & Serial < 732  ~ "Dec-17",
    Serial >= 732 & Serial < 763  ~ "Jan-18",
    Serial >= 763 & Serial < 791  ~ "Feb-18",
    Serial >= 791 & Serial < 822  ~ "Mar-18",
    Serial >= 822 & Serial < 852  ~ "Apr-18",
    Serial >= 852 & Serial < 883  ~ "May-18",
    TRUE ~ NA_character_  # Handle cases outside the defined ranges
  ))
# View the updated dataframe
head(pangfamHigh2.18s)

# Plot using custom x-axis breaks and labels
pang.fam.high2.18s <- ggplot(pangfamHigh2.18s, aes(x = Serial, y = value, colour = Family, 
                                                   group = Family, label = Family)) + 
  geom_line(linewidth = 1) +  
  #directlabels::geom_dl(aes(label = Family),  method = list(cex = 0.9, "smart.grid")) + 
  scale_color_manual(values = family_colors.18s) + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, margin = margin(b = 15)), 
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 18, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 18, margin = margin(r = 5)),
    axis.title.x = element_blank(),#text(size = 18, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18, face = "italic"),
    legend.spacing.y = unit(1.3, "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black")) + 
  labs(x = "day", y = "Relative abundance", 
       title = "Microeukaryotes: Pangasius ponds") + 
  scale_x_continuous(breaks = month_breaks, labels = month_labels)
pang.fam.high2.18s

##### Interactive plot----
ggplotly(pang.fam.high2.18s)

### Tilapia----
# Subset healthy pangasius (remove algal bloom) & also remove ponds with no fish
pseq_tila.18s <- ps.18s %>%
  ps_filter(Sampling_point > 1,
            Crop == "Tilapia",
            Crop_species != "Gulsha.Carp",
            Crop_species != "Gulsha.Pabda") %>% 
  tax_fix()
pseq_tila.18s #  4449 taxa and 370 samples

##### convert reads to relative abundance
tila.rel.18s <- transform_sample_counts(pseq_tila.18s, function(x){x / sum(x)})

# Aggregate at Family rank
#tila.glom.18s <- tax_glom(tila.rel.18s, taxrank = "Family")
#saveRDS(tila.glom.18s, "pseq_tila_18s_aggregated_at_family.rds")
tila.glom.18s <- readRDS("pseq_tila_18s_aggregated_at_family.rds")
tila.glom.18s # 305 taxa and 370 samples 

# Extract the taxonomic information
tila.tax.18s <-as.data.frame(tax_table(tila.glom.18s)[,7]) # 7th column (Family)
tila.tax.18s$ASV <- rownames(tila.tax.18s)

# Extract the ASV table
tila_ASV.18s <- as.data.frame(otu_table(tila.glom.18s))
# Transpose the ASV table
tila_ASVt.18s <- t(tila_ASV.18s) # transpose

#####Create a data.frame cols=ASV, rows= samples -> add date to this and estimate the sin/cos thing 
tila.tax.18s <- as.data.frame(tax_table(tila.glom.18s)[,7]) ###Keep this one to add to the melted object

tila_ASV.18s <- as.data.frame(otu_table(tila.glom.18s))

# Transpose
ASV_tila_18s_trans <- as.data.frame(x = t(tila_ASV.18s), stringsAsFactors = FALSE)
md_to_add.tila.18s <- as.data.frame(sample_data(tila.glom.18s))[,c(12)] #column number for Sampling_date in metadata

final_2a.tila.18s <- cbind(ASV_tila_18s_trans, md_to_add.tila.18s)
final_2a.tila.18s$Sampling_date <- as.Date(final_2a.tila.18s$Sampling_date,"%d/%m/%Y") # should i change as %Y/%m/%d

final_2a.tila_sort.18s <- final_2a.tila.18s[order(final_2a.tila.18s$Sampling_date),] ##sort from least recent to most recent
final_2a.tila_sort.18s$Sampling_date <- as.Date(final_2a.tila_sort.18s$Sampling_date,"%Y-%m-%d")

##Add the Julian day # what is this
final_2a.tila_sort.18s$Jul <- format(final_2a.tila_sort.18s$Sampling_date, "%j")

###Serial Day
#Starting from
tila_start.18s <- as.Date('2016-01-01') # should be first day of the sampling year

final_2a.tila_sort.18s$Serial <- as.numeric(-difftime(tila_start.18s,final_2a.tila_sort.18s$Sampling_date)+1)

##Days since winter solstice (21 Dic)
##Add the winter solstice 
final_2a.tila_sort.18s$wintersols[final_2a.tila_sort.18s$Sampling_date>as.Date("2016-01-01") & final_2a.tila_sort.18s$Sampling_date<as.Date("2016-12-21")] <- "2015-12-21"
final_2a.tila_sort.18s$wintersols[final_2a.tila_sort.18s$Sampling_date>as.Date("2017-01-01") & final_2a.tila_sort.18s$Sampling_date<as.Date("2017-12-21")] <- "2016-12-21"
final_2a.tila_sort.18s$wintersols[final_2a.tila_sort.18s$Sampling_date>as.Date("2018-01-01") & final_2a.tila_sort.18s$Sampling_date<as.Date("2018-12-21")] <- "2017-12-21"

final_2a.tila_sort.18s$DaysSincewintersols <- as.numeric(-difftime(final_2a.tila_sort.18s$wintersols,final_2a.tila_sort.18s$Sampling_date) + 1) #Add the numeric column
# This has missing values after 345, 

# Remove Rows with Missing Values:
final_2a.tila_sort.18s <- final_2a.tila_sort.18s[!is.na(final_2a.tila_sort.18s$DaysSincewintersols), ]

## Convert Julian date to numeric:
final_2a.tila_sort.18s$Jul <- as.numeric(final_2a.tila_sort.18s$Jul)

# Generate sine and cosine terms for harmonic regression:
final_2a.tila_sort.18s$cosDX1 <- as.numeric(cos(2*pi*((final_2a.tila_sort.18s$DaysSincewintersols)/365))) #Peaking in midwinter
final_2a.tila_sort.18s$sinDX1 <- as.numeric(sin(2*pi*((final_2a.tila_sort.18s$Jul)/365))) #Peaking other time
# check final_2a.tila_sort.18s for number of variables which will be needed below. Remove last 7 columns
#[306] "Sampling_date" [307] "Jul"[308]"Serial" [309]"wintersols" [310]"DaysSincewintersols" [311] "cosDX1" [312] "senDX1"

# Fit harmonic regression models:
storage.tila.18s <- list()
for(i in names(final_2a.tila_sort.18s)[1:305]){ #number of taxa present in final_2a.tila_sort.18s (312) -7
  storage.tila.18s[[i]] <- lm(get(i) ~ sinDX1+cosDX1, final_2a.tila_sort.18s)
}

# Extract p-values for significance testing:
pvalues_lm.tila.18s <- sapply(storage.tila.18s, function(x){
  ff <- summary(x)$fstatistic
  pf(ff[1], df1 = ff[2], df2 = ff[3], lower.tail = FALSE)
})

# Adjust the p-values using the BH method
pvalues_lm.tila_adjusted.18s <- p.adjust(pvalues_lm.tila.18s, method = "BH")

# Filter significant taxa based on adjusted p-values
pvalues_lm.tila_adjusted05.18s <- pvalues_lm.tila_adjusted.18s[pvalues_lm.tila_adjusted.18s < 0.05] # 71

# Extract significant taxa names after p-value adjustment
sign_names.tila.18s <- str_remove(names(pvalues_lm.tila_adjusted05.18s), ".value")

#Number of the column that matches a list of names, get the position of the names <0.01
signll.tila.18s <- match(sign_names.tila.18s, names(final_2a.tila_sort.18s))

#subset list of lists 
storage05.tila.18s <- storage.tila.18s[signll.tila.18s]

#I need the constant, sin and cos of the models
#Get the constant, sin and cos of the models:
# check if storage05.tila contains any NAs and if so, Filter out NA values from storage05.tila
filtered_storage.tila.18s <- storage05.tila.18s[!sapply(storage05.tila.18s, is.null)]

# Apply sapply to filtered storage
coeffs05.tila.18s <- data.frame(coeffs = sapply(filtered_storage.tila.18s, function(item) item$coefficients))

#coeffs05.tila<-data.frame(coeffs=sapply(storage05.tila, FUN=function(item){item$coefficients}))
colnames(coeffs05.tila.18s) <- str_remove(colnames(coeffs05.tila.18s), "coeffs.")

#Subset original day, sin, cos
newdf_model.tila.18s <- final_2a.tila_sort.18s[, c("Serial","cosDX1","sinDX1")]

# Automatize this :P 
# use the function below to automatically do it, and add the value in newdf_model.tila.18s table
automate_harmonic_regression.tila.18s <- function(final_2a.tila_sort.18s, 
                                                  coeffs05.tila.18s, sin_col, cos_col, prefix = "ASV_") {
  # Subset original sin and cos columns along with Serial column
  newdf_model.tila.18s <- final_2a.tila_sort.18s[, c("Serial", sin_col, cos_col)]
  
  # Automate creation of new columns for each ASV
  for (col_name in colnames(coeffs05.tila.18s)) {
    if (startsWith(col_name, "ASV_")) {
      new_col_name <- paste0(prefix, str_extract(col_name, "\\d+"))
      new_col <- coeffs05.tila.18s[[col_name]][1] + 
        (newdf_model.tila.18s$sinDX1 * coeffs05.tila.18s[[col_name]][2]) + 
        (newdf_model.tila.18s$cosDX1 * coeffs05.tila.18s[[col_name]][3])
      newdf_model.tila.18s[[new_col_name]] <- new_col
    }
  }
  
  return(newdf_model.tila.18s)
}

newdf_model.tila.18s <- automate_harmonic_regression.tila.18s(final_2a.tila_sort.18s, coeffs05.tila.18s, "sinDX1", "cosDX1")
# check newdf_model.tila.18s for number of variables
Tomelt.tila.18s <- (newdf_model.tila.18s[c(1,4:72)]) # for p < 0.05, 69 families showed seasonal trend

#Add taxa to this melted 
newdf_model.tila_melted.18s <- melt(Tomelt.tila.18s, id.vars="Serial") # 23943 observation

#Add taxa extracting sign_names.tila variable from the Phyloseq object, thereafter extract the tax table only with genus 
##### Subset from physeq 
fam_sign05.tila.18s <- prune_taxa(sign_names.tila.18s, tila.glom.18s)
# 37 taxa and 370 samples # showed seasonal trend

fam_tila.sign.18s <- data.frame(tax_table(fam_sign05.tila.18s))[,7,drop=FALSE] # 7th column (Family)

fam_tila.sign.18s$ASV <- rownames(fam_tila.sign.18s) #fam_tila.sign.18s contains the neccesary to add taxa to melted

#change column name from "variable" to "ASV" in newdf_model.tila_melted 
names(newdf_model.tila_melted.18s)[2] <- "ASV"

#merge by ASV and add the column Genus 
newdf_model.tila_melted_taxa.18s <- merge(newdf_model.tila_melted.18s, 
                                          fam_tila.sign.18s[, c("Family", "ASV")], by="ASV")
# for family, 12839 observation

#Subset 1 to 365 (one year)
newdf_model.tila_melted_taxa365.18s <- newdf_model.tila_melted_taxa.18s[newdf_model.tila_melted_taxa.18s$Serial > 365 & newdf_model.tila_melted_taxa.18s$Serial < 730,]

# consider all the years
newdf_model.tila_melted_taxa365.18s$Serial365 <- (newdf_model.tila_melted_taxa365.18s$Serial-364)

# Take those with a mean relative abundance >1%
# Calculate the mean relative abundance for each ASV
mean_abundance_tila.18s <- colMeans(tila_ASVt.18s)

# Filter ASVs with mean relative abundance > 2% (0.02)
h.tila.18s <- colnames(tila_ASVt.18s)[mean_abundance_tila.18s > 0.01]

# Remove any NAs
h.tila.18s <- h.tila.18s[!is.na(h.tila.18s)]

### I want to plot all years instead of just one year
# instead of just one year, if i want to plot all the sampling months
tilafamHigh2.18s <- newdf_model.tila_melted_taxa.18s[newdf_model.tila_melted_taxa.18s$ASV %in% h.tila.18s, ] # 2082 observation

# check unique family names in the dataframe
unique(tilafamHigh2.18s$Family) # 9 families


#### Plot month wise----
# Add a new column "Sampling_month" based on the Serial column
tilafamHigh2.18s <- tilafamHigh2.18s %>%
  mutate(Sampling_month = case_when(
    Serial >= 275 & Serial < 306  ~ "Oct-16",
    Serial >= 306 & Serial < 336  ~ "Nov-16",
    Serial >= 336 & Serial < 367  ~ "Dec-16",
    Serial >= 367 & Serial < 398  ~ "Jan-17",
    Serial >= 398 & Serial < 426  ~ "Feb-17",
    Serial >= 426 & Serial < 457  ~ "Mar-17",
    Serial >= 457 & Serial < 487  ~ "Apr-17",
    Serial >= 487 & Serial < 518  ~ "May-17",
    Serial >= 518 & Serial < 548  ~ "Jun-17",
    Serial >= 548 & Serial < 579  ~ "Jul-17",
    Serial >= 579 & Serial < 610  ~ "Aug-17",
    Serial >= 610 & Serial < 640  ~ "Sep-17",
    Serial >= 640 & Serial < 671  ~ "Oct-17",
    Serial >= 671 & Serial < 701  ~ "Nov-17",
    Serial >= 701 & Serial < 732  ~ "Dec-17",
    Serial >= 732 & Serial < 763  ~ "Jan-18",
    Serial >= 763 & Serial < 791  ~ "Feb-18",
    Serial >= 791 & Serial < 822  ~ "Mar-18",
    Serial >= 822 & Serial < 852  ~ "Apr-18",
    Serial >= 852 & Serial < 883  ~ "May-18",
    TRUE ~ NA_character_  # Handle cases outside the defined ranges
  ))

# View the updated dataframe
head(tilafamHigh2.18s)

# Plot using custom x-axis breaks and labels
tila.fam.high2.18s <- ggplot(tilafamHigh2.18s, aes(x = Serial, y = value, colour = Family, group = Family, label = Family)) + 
  geom_line(linewidth = 1) +
  #directlabels::geom_dl(aes(label = Family),  method = list(cex = 0.9, "smart.grid")) + 
  scale_color_manual(values = family_colors.18s) + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, margin = margin(b = 15)), 
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 18, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 18, margin = margin(r = 5)),
    axis.title.x = element_blank(),#text(size = 18, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18, face = "italic"),
    legend.spacing.y = unit(1.3, "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black")) + 
  labs(x = "day", y = "Relative abundance", 
       title = "Microeukaryotes: Tilapia ponds") + 
  scale_x_continuous(breaks = month_breaks, labels = month_labels)
tila.fam.high2.18s

##### Interactive plot----
ggplotly(tila.fam.high2.18s)

# Finalising plots----
# Combine all alpha diversity plots
season.alpha <- cowplot::plot_grid(p.sm.s +theme(legend.position = "none", axis.text.x = element_blank()),
                                   t.sm.s +theme(legend.position = "none", axis.text.x = element_blank()),
                                   p.sm.s.18s +theme(legend.position = "none", axis.text.x = element_blank()),
                                   t.sm.s.18s +theme(legend.position = "none"),# axis.text.x = element_blank()),
                                   labels = c("A", "B", "C", "D"),
                                   rel_heights = c(0.87,0.87,0.87, 1),
                                   ncol = 1)
# Combine harmonic regression
season.hr <- cowplot::plot_grid(pang.fam.high2 +theme(legend.position = "none", axis.text.x = element_blank()),#, plot.title = element_blank()),
                                tila.fam.high2 +theme(legend.position = "none",axis.text.x = element_blank()), 
                                pang.fam.high2.18s +theme(legend.position = "none", axis.text.x = element_blank()),#, plot.title = element_blank()),
                                tila.fam.high2.18s +theme(legend.position = "none"),#axis.text.x = element_blank()),
                                labels = c("E", "F", "G", "H"),
                                rel_heights = c(0.87,0.87,0.87, 1),
                                align = "v",
                                ncol = 1)
# Combine alpha and harmonic regression
comb.season <- cowplot::plot_grid(season.alpha,
                                  season.hr,
                                  ncol = 2)
# Get 16s plot legends
legend.season <- get_legend(p.sm.s)
legend.pang.fam <- get_legend(pang.fam.high2)
legend.tila.fam <- get_legend(tila.fam.high2)

# Combine 16s legends
legend.16s <- cowplot::plot_grid(legend.season,
                             legend.pang.fam,
                             legend.tila.fam, rel_heights = c(2,3,5),
                             ncol = 1)
# Get 18s plot legends
legend.season.18s <- get_legend(p.sm.s.18s)
legend.pang.fam.18s <- get_legend(pang.fam.high2.18s)
legend.tila.fam.18s <- get_legend(tila.fam.high2.18s)

# Combine legends
legend.18s <- cowplot::plot_grid(legend.season.18s,
                                 legend.pang.fam.18s,
                                 legend.tila.fam.18s, rel_heights = c(2,4,4),
                                 ncol = 1)
# Combine 16s and 18s legend
comb.legend <- cowplot::plot_grid(legend.16s,
                                  legend.18s, ncol = 1)

# Final figure
comb.seasonal.trend <- cowplot::plot_grid(comb.season,
                                          comb.legend,
                                          ncol = 2,
                                          rel_widths = c(2, 0.4))
comb.seasonal.trend
# Save as 2200*2000 and give a final touch for legends in inkscape











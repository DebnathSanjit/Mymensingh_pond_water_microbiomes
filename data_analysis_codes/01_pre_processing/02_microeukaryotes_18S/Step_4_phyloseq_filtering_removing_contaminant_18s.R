# 04/04/2024
# Sanjit Debnath

# Load Libraries
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")
library(microViz); packageVersion("microViz")
library(tidyverse); packageVersion("tidyverse")

## Setup working dictionary first----
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/pre_processing_18S/primer_removed")

#small function to tidy ps object after subsetting to remove any zero taxa or samples, 
#optionally transforms counts to relative abundance (JMcM)
tidyPS <- function(ps, RA = FALSE){
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  if(RA == TRUE){
    ps <- transform_sample_counts(ps, function(x){x / sum(x)})
  }
  return(ps)}

# Set theme for ggplot
theme_set(theme_bw())

# Load the phyloseq object
ps.18s <- readRDS("phyloseq_18S_90_150bp_pr2_v5.0.0_20240403.rds")
ps.18s # 9345 taxa and 955 samples
ps <- ps.18s

# Remove likely contaminates (in the extraction and PCR stage) from the amplicon data. 
# Then remove all Eukaryotes, Chloroplasts, Mitochondria and NA at Kingdom level, then 
# and also do some prevalence filtering. for prevalence filtering, can use 0.01-0.05 
# (1-5% of samples) threshold. Then add the refseq slot into the phyloseq object. 
# Then we can make the phylogenetic tree.

# Prevalence filtering: https://f1000research.com/articles/5-1492
# remove contaminant: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

#Inspect Library Sizes (number of reads) in each sample to check whether that sample was a true positive sample or a negative control
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

# Plot library size
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_neg)) + 
        geom_point()

# Check the contaminants for Water samples + mocks + negs
# prevalence filter---- (https://f1000research.com/articles/5-1492)
#This function will remove taxa (ASVs) with low prevalence, where prevalence is the fraction of total samples in which an ASV is observed.
# Compute prevalence of each feature, store as data.frame
ps_water_prevdf <- apply(X = otu_table(ps), MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
ps_water_prevdf <- data.frame(Prevalence = ps_water_prevdf, TotalAbundance = taxa_sums(ps), tax_table(ps))

# Define prevalence threshold as present ASV in 2 samples or more
prevalenceThreshold <- 2
keepTaxa <- rownames(ps_water_prevdf)[(ps_water_prevdf$Prevalence >= prevalenceThreshold)]
ps_water2 <- prune_taxa(keepTaxa, ps)
ps_water2    #5503 taxa and 955 samples

# Identify Contaminants - Prevalence
# for checking and removing contaminants, as prevalence seems works better, 
# I will only use prevalence filter (https://benjjneb.github.io/decontam/vignettes/decontam_intro.html)
sample_data(ps_water2)$is.neg <- sample_data(ps_water2)$Sample_or_neg == "neg"
ps_water_contamdf.prev <- isContaminant(ps_water2, method="prevalence", neg="is.neg") #default prevalence 0.1
table(ps_water_contamdf.prev$contaminant)
#FALSE  TRUE 
#5495     8 

# check the contaminats
which(ps_water_contamdf.prev$contaminant)
# 232  562 2829 2839 3382 4457 5149 5462

# Increase the prevalence to 0.5
ps_water_contamdf.prev05 <- isContaminant(ps_water2, method="prevalence", neg="is.neg", threshold=0.5)
table(ps_water_contamdf.prev05$contaminant)
#FALSE  TRUE 
#5483    20

which(ps_water_contamdf.prev05$contaminant)
# 80   81  101  109  156  175  232  562  760  969 1231 1489 1839 2339 2829 2839 3382
# 4457 5149 5462

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps_water2, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_neg == "neg", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_neg == "sample", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=ps_water_contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# Make data.frame of prevalence in positive and negative samples (threshold=0.5)
df.pa05 <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=ps_water_contamdf.prev05$contaminant)
ggplot(data=df.pa05, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")


# Remove likely contaminants from both phyloseq object
ps_water2  # 5503 taxa and 955 samples

#get vector of contaminants
wat.contam <- ps_water2@otu_table[ps_water_contamdf.prev05$contaminant,] %>% row.names() 
wat.contam
#[1] "ASV_80"   "ASV_81"   "ASV_101"  "ASV_109"  "ASV_156"  "ASV_175"  "ASV_232" 
#[8] "ASV_563"  "ASV_761"  "ASV_970"  "ASV_1233" "ASV_1491" "ASV_1844" "ASV_2347"
#[15] "ASV_2848" "ASV_2859" "ASV_3423" "ASV_4627" "ASV_5781" "ASV_7226"

view(ps_water2@tax_table[wat.contam,]) # to view, need to load tidyverse package


# We can manually check these asv if they are actual contaminat or not
### But we can also just remove these asvs without checking, as they don't make much difference 
# Remove contaminants----
water.nocontam <- prune_taxa(!(taxa_names(ps_water2) %in% wat.contam), ps_water2)
water.nocontam # 5483 taxa and 955 samples

# Save this as RDS file
saveRDS(water.nocontam, "phyloseq_18S_no_contam_pr2_90-150_20240403.rds")

# load latest phyloseq obect----
ps1 <- readRDS("phyloseq_18S_no_contam_pr2_90-150_20240403.rds")
ps1 # 5483 taxa and 955 samples

# subset only mock to check the accuracy later----
ps_mock <- ps1 %>% 
  ps_filter(Sample_type %in% c("Mock.DNA", "Mock.extr"))
ps_mock #5 taxa and 12 samples
saveRDS(ps_mock, "phyloseq_18S_mock_pr2_90-150_20240403.rds")

# as previously we have removed contaminants using negs, and also I have seperated 
#mock for checking accuracy later, we do not need them any more, so remove them
ps.no.control <- ps1 %>% 
  ps_filter(Sample_type %in% c ("Algal.bloom", "Pond_water"))
ps.no.control #  5468 taxa and 920 samples
saveRDS(ps.no.control, "phyloseq_18S_no_control_pr2_90-150_20240403.rds")

# load the rds file
ps.no.control <- readRDS("phyloseq_18S_no_control_pr2_90-150_20240403.rds")
ps.no.control

## Remove non-target sequences----
#Jamie suggests removing anything Class == craniata which removes fish and humans. 
# Ash suggest removing any family == Teleostei

# #removing Craniata----
ps_no_Craniata <- ps.no.control %>% subset_taxa(Class !="Craniata"  | is.na(Class))
ps_no_Craniata # 5456 taxa and 920 samples

#ps_no_Craniata %>%
#  sample_sums() %>%
#  quantile(probs = c(0, 0.05, 0.5, 0.95, 1))

#removing Teleostei
ps_no_Teleostei <- ps_no_Craniata %>% subset_taxa(Family !="Teleostei"  | is.na(Family))
ps_no_Teleostei # 5456 taxa and 920 samples

#removing NA
ps_no_NA <- ps_no_Teleostei %>% subset_taxa(Domain !='NA') 
ps_no_NA # 5404 taxa and 920 samples

# removing bacteria
ps_no_bac <- ps_no_NA %>% subset_taxa(Domain != "Bacteria")
ps_no_bac #5390 taxa and 920 samples

# check the count number per sample and add in the metadata
# Calculate the count number per sample
count_per_sample <- sample_sums(ps_no_bac)

# Add the count number to the sample metadata
sample_data(ps_no_bac)$count_per_sample <- count_per_sample

# Check the updated sample metadata
sample_data(ps_no_bac)

#ps_no_bac %>%
#  sample_sums() %>%
#  quantile(probs = c(0, 0.05, 0.5, 0.95, 1))

#Remove low read count samples
ps_high_counts <- prune_samples(sample_sums(ps_no_bac)>=2000, ps_no_bac)  #remove low read count (less than 2000 reads) samples that have failed
ps_high_counts #  5390 taxa and 872 samples 

# Rename this and saved as all_filtered as everything has been filtered from this 
ps_all_filtered <- ps_high_counts
ps_all_filtered #  5390 taxa and 872 samples

#adding refseq slot back again to phyloseq object
dna <- Biostrings::readDNAStringSet("asvs_fasta_18S_90-150bp_20240404.fasta")

ps <- merge_phyloseq(ps_all_filtered, dna)
ps

# Save this final phyloseq object and construct phologenetic tree for further downstream analysis
saveRDS(ps, "phyloseq_18S_all_filtered_pr2_90-150bp_20240404.rds")

# send this to Ash to make the phylogenetic tree


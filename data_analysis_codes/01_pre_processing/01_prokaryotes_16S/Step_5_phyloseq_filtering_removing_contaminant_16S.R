# Remove likely contaminates (in the extraction and PCR stage) from the amplicon data. 
# Then remove all Eukaryotes, Chloroplasts, Mitochondria and NA at Kingdom level, then 
# and also do some prevalence filtering. for prevalence filtering, can use 0.01-0.05 
# (1-5% of samples) threshold. Then add the refseq slot into the phyloseq object. 
# Then we can make the phylogenetic tree.

# Prevalence filtering: https://f1000research.com/articles/5-1492
# remove contaminant: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

# Load Libraries----
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")
library(microViz); packageVersion("microViz")
library(tidyverse); packageVersion("tidyverse")

## Setup working dictionary first----
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/Pre_processing_16s")

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
ps <- readRDS("phyloseq_all_samples_20240105.rds")
ps # 38250 taxa and 1066 samples

#Inspect Library Sizes (number of reads) in each sample to check whether that sample was a true positive sample or a negative control
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

# Plot library size
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_neg)) + 
        geom_point()

# Check the contaminants seperately for batch 4 (B4) and all other batches, as batch 4 has different type samples
# 1. Water samples + mocks + negs
ps_water <- ps %>% 
  ps_filter(Amplicon_batch != 4)
ps_water   #32706 taxa and 963 samples

# 2. Incubations + mock +negs (only batch 4)
ps_incu <- ps %>% 
  ps_filter(Amplicon_batch == 4)
ps_incu       #11014 taxa and 103 samples

# prevalence filter---- (https://f1000research.com/articles/5-1492)
#This function will remove taxa (ASVs) with low prevalence, where prevalence is the fraction of total samples in which an ASV is observed.
# Compute prevalence of each feature, store as data.frame
ps_water_prevdf <- apply(X = otu_table(ps_water), MARGIN = ifelse(taxa_are_rows(ps_water), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
ps_water_prevdf <- data.frame(Prevalence = ps_water_prevdf, TotalAbundance = taxa_sums(ps_water), tax_table(ps_water))

# Define prevalence threshold as present ASV in 2 samples or more
prevalenceThreshold <- 2
keepTaxa <- rownames(ps_water_prevdf)[(ps_water_prevdf$Prevalence >= prevalenceThreshold)]
ps_water2 <- prune_taxa(keepTaxa, ps_water)
ps_water2    #11687 taxa and 963 samples

# Identify Contaminants - Prevalence
# for checking and removing contaminants, as prevalence seems works better, 
# I will only use prevalence filter (https://benjjneb.github.io/decontam/vignettes/decontam_intro.html)
sample_data(ps_water2)$is.neg <- sample_data(ps_water2)$Sample_or_neg == "neg"
ps_water_contamdf.prev <- isContaminant(ps_water2, method="prevalence", neg="is.neg") #default prevalence 0.1
table(ps_water_contamdf.prev$contaminant)
#FALSE  TRUE 
#11667    20 

# check the contaminats
head(which(ps_water_contamdf.prev$contaminant))
# 68   85  144  213  506 1427

# Increase the prevalence to 0.5
ps_water_contamdf.prev05 <- isContaminant(ps_water2, method="prevalence", neg="is.neg", threshold=0.5)
table(ps_water_contamdf.prev05$contaminant)
#FALSE  TRUE 
#11662    25 

head(which(ps_water_contamdf.prev05$contaminant))
# 68  85  92 144 213 506

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


### Repeat the steps for batch 4----
# Compute prevalence of each feature, store as data.frame
ps_incu_prevdf <- apply(X = otu_table(ps_incu), MARGIN = ifelse(taxa_are_rows(ps_incu), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
ps_incu_prevdf <- data.frame(Prevalence = ps_incu_prevdf, TotalAbundance = taxa_sums(ps_incu), tax_table(ps_incu))

# Define prevalence threshold as present ASV in 2 samples or more
prevalenceThreshold <- 2
keepTaxa_b4 <- rownames(ps_incu_prevdf)[(ps_incu_prevdf$Prevalence >= prevalenceThreshold)]
ps_incu2 <- prune_taxa(keepTaxa_b4, ps_incu)
ps_incu2    #4286 taxa and 103 samples

# Identify Contaminants - Prevalence
# for checking and removing contaminants, as prevalence seems works better, 
# I will only use prevalence filter (https://benjjneb.github.io/decontam/vignettes/decontam_intro.html)
sample_data(ps_incu2)$is.neg <- sample_data(ps_incu2)$Sample_or_neg == "neg"
ps_incu_contamdf.prev <- isContaminant(ps_incu2, method="prevalence", neg="is.neg") #default prevalence 0.1
table(ps_incu_contamdf.prev$contaminant)
#FALSE  TRUE 
#4274    12  

# check the contaminats
head(which(ps_incu_contamdf.prev$contaminant))
# 212  477 1171 1389 1406 1896

# Increase the prevalence to 0.5
ps_incu_contamdf.prev05 <- isContaminant(ps_incu2, method="prevalence", neg="is.neg", threshold=0.5)
table(ps_incu_contamdf.prev05$contaminant)
#FALSE  TRUE 
#4274    12 

head(which(ps_incu_contamdf.prev05$contaminant))
# 212  477 1171 1389 1406 1896

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_b4 <- transform_sample_counts(ps_incu2, function(abund) 1*(abund>0))
ps.pa.neg_b4 <- prune_samples(sample_data(ps.pa_b4)$Sample_or_neg == "neg", ps.pa_b4)
ps.pa.pos_b4 <- prune_samples(sample_data(ps.pa_b4)$Sample_or_neg == "sample", ps.pa_b4)

# Make data.frame of prevalence in positive and negative samples
df.pa_b4 <- data.frame(pa.pos=taxa_sums(ps.pa.pos_b4), pa.neg=taxa_sums(ps.pa.neg_b4),
                    contaminant=ps_incu_contamdf.prev$contaminant)
ggplot(data=df.pa_b4, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# Make data.frame of prevalence in positive and negative samples (threshold=0.5)
df.pa05_b4 <- data.frame(pa.pos=taxa_sums(ps.pa.pos_b4), pa.neg=taxa_sums(ps.pa.neg_b4),
                      contaminant=ps_incu_contamdf.prev05$contaminant)
ggplot(data=df.pa05_b4, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")


# Remove likely contaminants from both phyloseq object
ps_water2  #11687 taxa and 963 samples
ps_incu2   #4286 taxa and 103 samples

#get vector of contaminants
wat.contam <- ps_water2@otu_table[ps_water_contamdf.prev05$contaminant,] %>% row.names() 
wat.contam
#[1] "ASV_68"    "ASV_85"    "ASV_92"    "ASV_144"   "ASV_215"   "ASV_510"  
#[7] "ASV_739"   "ASV_1466"  "ASV_2072"  "ASV_2102"  "ASV_2443"  "ASV_2697" 
#[13] "ASV_2953"  "ASV_4083"  "ASV_4314"  "ASV_4461"  "ASV_4742"  "ASV_5548" 
#[19] "ASV_5557"  "ASV_5858"  "ASV_5964"  "ASV_7250"  "ASV_8780"  "ASV_12225"
#[25] "ASV_15952"
view(ps_water2@tax_table[wat.contam,]) # to view, need to load tidyverse package

incu.contam <- ps_incu2@otu_table[ps_incu_contamdf.prev05$contaminant,] %>% row.names() 
incu.contam
#[1] "ASV_215"  "ASV_510"  "ASV_1466" "ASV_1828" "ASV_1851" "ASV_2824"
#[7] "ASV_2975" "ASV_3570" "ASV_5699" "ASV_6594" "ASV_7250" "ASV_8108"
view(ps_incu2@tax_table[incu.contam,]) # to view, need to load tidyverse package

# We can manually check these asv if they are actual contaminat or not

# If want to check individual likely contaminats----
# Check possible Contaminants in controls (mock, negs) as then often appear in control generally with higher counts compared to samples.
# check consistency of possible Contaminants in controls and samples as contaminates may show varying counts in different samples, while genuine biological signals are often more consistent
# check the presence of possible contaminants in negs, as they may be introduced during the DNA extraction or PCR steps.

# ASV_68_Escherichia-Shigella
ASV_68 <- prune_taxa("ASV_68", ps_water2) %>% tidyPS() %>% otu_table() 
ASV_68
# present in 10 mocks with high reads, in two negs. possible contam

# ASV_85_Paucibacter
ASV_85 <- prune_taxa("ASV_85", ps_water2) %>% tidyPS() %>% otu_table()
ASV_85
# present in 1 extract (5 reads), in 4 negs(3, 4, 5, 17,), in 167 samples relatively higher reads than neg. Not contam?

# ASV_92-Listeria
ASV_92 <- prune_taxa("ASV_92", ps_water2) %>% tidyPS() %>% otu_table()
ASV_92
# present in 10 mocks with very high read (2251-6686), in one extr neg with low reads (4), possible contam?

# ASV_144-Lactiplantibacillus
ASV_144 <- prune_taxa("ASV_144", ps_water2) %>% tidyPS() %>% otu_table()
ASV_144
# found in 4 pcr negs (5, 12, 546, 873,), 1 extract neg (122) and in 177 samples. possible contam?

#ASV_215-Ralstonia
ASV_215 <- prune_taxa("ASV_215", ps_water2) %>% tidyPS() %>% otu_table()
ASV_215 

# ASV_510-Sphingomonas
ASV_510 <- prune_taxa("ASV_510", ps_water2) %>% tidyPS() %>% otu_table()

# ASV_739-Pseudomonas
ASV_739 <- prune_taxa("ASV_739", ps_water2) %>% tidyPS() %>% otu_table()

#ASV_1466-Bradyrhizobium
ASV_1466 <- prune_taxa("ASV_1466", ps_water2) %>% tidyPS() %>% otu_table()
ASV_1466

#ASV_2072-Rhizobiaceae 
ASV_2072 <- prune_taxa("ASV_2072", ps_water2) %>% tidyPS() %>% otu_table()
ASV_2072

#ASV_2102-Alteromonas   
ASV_2102 <- prune_taxa("ASV_2102", ps_water2) %>% tidyPS() %>% otu_table()
ASV_2102

#ASV_2443   NA
ASV_2443 <- prune_taxa("ASV_2443", ps_water2) %>% tidyPS() %>% otu_table()
ASV_2443

#ASV_2697   "Methylomonas" 
ASV_2697 <- prune_taxa("ASV_2697", ps_water2) %>% tidyPS() %>% otu_table()
ASV_2697

#ASV_2953   "Pseudomonas"  
ASV_2953 <- prune_taxa("ASV_2953", ps_water2) %>% tidyPS() %>% otu_table()
ASV_2953

#ASV_4083   "NS9 marine group" 
ASV_4083 <- prune_taxa("ASV_4083", ps_water2) %>% tidyPS() %>% otu_table()
ASV_4083

#ASV_4314   "Alishewanella"  
ASV_4314 <- prune_taxa("ASV_4314", ps_water2) %>% tidyPS() %>% otu_table()
ASV_4314

#ASV_4461   "Algoriphagus"  
ASV_4461 <- prune_taxa("ASV_4461", ps_water2) %>% tidyPS() %>% otu_table()
ASV_4461

#ASV_4742   "Sphingomonas"  
ASV_4742 <- prune_taxa("ASV_4742", ps_water2) %>% tidyPS() %>% otu_table()
ASV_4742

#ASV_5548   "Thalassospira" 
ASV_5548 <- prune_taxa("ASV_5548", ps_water2) %>% tidyPS() %>% otu_table()
ASV_5548

#ASV_5557   "Shimia"    
ASV_5557 <- prune_taxa("ASV_5557", ps_water2) %>% tidyPS() %>% otu_table()
ASV_5557

#ASV_5858   "Rheinheimera"  
ASV_5858 <- prune_taxa("ASV_5858", ps_water2) %>% tidyPS() %>% otu_table()
ASV_5858 

#ASV_5964   "Rhodobacteraceae"
ASV_5964 <- prune_taxa("ASV_5964", ps_water2) %>% tidyPS() %>% otu_table()
ASV_5964

#ASV_7250    "Sediminibacterium"  
ASV_7250 <- prune_taxa("ASV_7250", ps_water2) %>% tidyPS() %>% otu_table()
ASV_7250

#ASV_8780    "Pseudomonas"    
ASV_8780 <- prune_taxa("ASV_8780", ps_water2) %>% tidyPS() %>% otu_table()
ASV_8780

#ASV_12225   "Methylophaga" 
ASV_12225 <- prune_taxa("ASV_12225", ps_water2) %>% tidyPS() %>% otu_table()
ASV_12225

#ASV_15952   "Balneola"
ASV_15952 <- prune_taxa("ASV_15952", ps_water2) %>% tidyPS() %>% otu_table()
ASV_15952

# similarly need to manually check for other phyloseq object (incubation, batch 4)

### But we can also just remove these asvs without checking, as they don't make much difference 
# Remove contaminants----
water.nocontam <- prune_taxa(!(taxa_names(ps_water2) %in% wat.contam), ps_water2)
water.nocontam # 11662 taxa and 963 samples

incu.nocontam <- prune_taxa(!(taxa_names(ps_incu2) %in% incu.contam), ps_incu2)
incu.nocontam # 4274 taxa and 103 samples

## Merge both phyloseq objects together
ps.nocontam <- merge_phyloseq(water.nocontam, incu.nocontam)
ps.nocontam  #13228 taxa and 1066 samples

# Save this as RDS file
#saveRDS(ps.nocontam, "phyloseq_no_contam_20240105.rds")

# load latest phyloseq obect----
ps1 <- readRDS("phyloseq_no_contam_20240105.rds")
ps1 # 13228 taxa and 1066 samples

# subset only mock to check the accuracy later----
ps_mock <- ps1 %>% 
  ps_filter(Sample_type == "Mock")
ps_mock
#saveRDS(ps_mock, "phyloseq_mock_20240106.rds")

# as previously we have removed contaminants using negs, and also I have seperated 
#mock for checking accuracy later, we do not need them any more, so remove them
ps.no.control <- ps1 %>% 
  ps_filter(Sample_type != "neg" & Sample_type != "Mock")
ps.no.control # 13226 taxa and 1032 samples

#removing eukaryotes:
ps_no_euk <- ps.no.control %>% subset_taxa(Kingdom !="Eukaryota"  | is.na(Kingdom))
ps_no_euk # 13226 taxa and 1032 samples
#ps.no.neg - ps_no_euk = no change of taxa as no Eukaryota was present

#removing chloroplast
ps_no_chl <- ps_no_euk %>% subset_taxa(Order !="Chloroplast"  | is.na(Order))
ps_no_chl
#ps_no_euk - ps_no_chl = 13226 - 638 = 12588 taxa and 1032 samples

#removing mitochondria
ps_no_mit <- ps_no_chl %>% subset_taxa(Family !='Mitochondria'  | is.na(Family)) 
ps_no_mit
# ps_no_chl - ps_no_mit = 12588 - 277 = 12311 taxa and 1032 samples

# removing NAs from phyloseq object
ps_no_NA <- ps_no_mit %>% subset_taxa(Kingdom !='NA') 
ps_no_NA
#ps_no_mit - ps_no_NA = 12311 - 5 = 12306 taxa and 1032 samples

#Remove low read count samples
ps_high_counts <- prune_samples(sample_sums(ps_no_NA)>=2000, ps_no_NA) %>% tidyPS() #remove low read count (less than 2000 reads) samples that have failed
ps_high_counts # 12304 taxa and 1000 samples 

# Save as rds file
#saveRDS(ps_high_counts, "phyloseq_high_counts_20240106.rds")
ps_high_counts <- readRDS("phyloseq_high_counts_20240106.rds")

# Check again for ASVs that don't occur in any sample (happens when samples are
# removed, especially mock community positive controls, which have different taxa)
#ps_high_counts2 <- filter_taxa(ps_high_counts, function(x) sum(x > 0) > 0, TRUE)
#ps_high_counts2

# Rename this and saved as all_filtered as everything has been filtered from this 
ps_all_filtered <- ps_high_counts
ps_all_filtered #12304 taxa and 1000 samples

# Add refseq slot back again to phyloseq object
dna <- Biostrings::readDNAStringSet("ASVs_fasta_all_batches_bayesian_20240101.fasta")

ps <- merge_phyloseq(ps_all_filtered, dna)
ps # 12304 taxa and 1000 samples

# Save this final phyloseq object and construct phologenetic tree for further downstream analysis
saveRDS(ps, "phyloseq_filtered_final_20240106.rds")

# sent this to Ashley to make the phylogenetic tree


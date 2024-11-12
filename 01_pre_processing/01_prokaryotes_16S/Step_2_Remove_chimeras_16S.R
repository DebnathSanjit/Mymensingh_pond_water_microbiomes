# 01/01/2024
# Sanjit Debnath

#Processing with DADA2 in R
# Workflow for amplicon data----
# First follow - https://benjjneb.github.io/dada2/tutorial.html
# Secondly follow - https://astrobiomike.github.io/amplicon/dada2_workflow_ex
# can have a look this video for explanation: https://www.youtube.com/watch?v=Qc3GF2RYUWw

## Load dada2 package
library(dada2); packageVersion("dada2")

## Setup working dictionary first----
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/Pre_processing_16S")

# Load sequence table----
seqtab_B1 <- readRDS("seqtab_B1_20231231.rds")
seqtab_B2 <- readRDS("seqtab_B2_20231231.rds")
seqtab_B4 <- readRDS("seqtab_B4_20231231.rds") # not part of this study
seqtab_B6 <- readRDS("seqtab_B6_20231231.rds")
seqtab_B7 <- readRDS("seqtab_B7_20240101.rds")
seqtab_B8 <- readRDS("seqtab_B8_20240101.rds")
seqtab_B9 <- readRDS("seqtab_B9_20240101.rds")
length(seqtab_B1)
length(seqtab_B2)
length(seqtab_B4) # not part of this study
length(seqtab_B6)
length(seqtab_B7)
length(seqtab_B8)
length(seqtab_B9)

#merge all sequence table from all sequence batch together
seqtabAll <- mergeSequenceTables(seqtab_B1, seqtab_B2, seqtab_B4, seqtab_B6, seqtab_B7, seqtab_B8, seqtab_B9)
#seqtabAll

#saveRDS(seqtabAll, "seqtab_All.rds")
seqtabAll <- readRDS("seqtab_All.rds")

# The approximate size for these amplicons should be 254 bp. Amplicons of 
# sizes other than this most likely to be non-specific, so need to remove them
# Remove non-target-length sequences from the sequence table and keep only those ranges 250-256 bp
seqtab2 <- seqtabAll[,nchar(colnames(seqtabAll)) %in% 250:256]
#seqtab2

#saveRDS(seqtab2, "seqtab_All_size_250-256.rds")
seqtabAll2 <- readRDS("seqtab_All_size_250-256.rds")

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtabAll)))
table(nchar(getSequences(seqtabAll2))) # We can see the difference

## Chimera identification/Remove chimeras----
seqtab.nochim <- removeBimeraDenovo(seqtabAll, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) # 1066 39302

sum(seqtab.nochim)/sum(seqtabAll)
#  0.9935064 # in this case we lost around 1% in terms of abundance

## Chimera identification/Remove chimeras----
seqtab.nochim2 <- removeBimeraDenovo(seqtabAll2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim2) # 1066 38250

sum(seqtab.nochim2)/sum(seqtabAll2)
#  0.9935028 # in this case we lost around 1% in terms of abundance

#sample.names <- row.names(seqtabAll2)

# Save as rds file
#saveRDS(seqtab.nochim2, "seqtab_nochimera_20240101.rds")
seqtab.nochim <- readRDS("seqtab_nochimera_20240101.rds")

# making sequence table from the separate batches already removed chimeras
# Load sequence tables with removed chimeras----
seqtab_B1_2 <- readRDS("seqtab_nochimera_B1_20240101.rds")
seqtab_B2_2 <- readRDS("seqtab_nochimera_B2_20240102.rds")
seqtab_B4_2 <- readRDS("seqtab_nochimera_B4_20240102.rds") # Not part of this study
seqtab_B6_2 <- readRDS("seqtab_nochimera_B6_20240102.rds")
seqtab_B7_2 <- readRDS("seqtab_nochimera_B7_20240102.rds")
seqtab_B8_2 <- readRDS("seqtab_nochimera_B8_20240102.rds")
seqtab_B9_2 <- readRDS("seqtab_nochimera_B9_20240104.rds")

#merge all sequence table from all sequence batch together
seqtab_nochimera <- mergeSequenceTables(seqtab_B1_2, seqtab_B2_2, seqtab_B4_2, seqtab_B6_2, seqtab_B7_2, seqtab_B8_2, seqtab_B9_2)
dim(seqtab_nochimera) # 1066 39441

#
#saveRDS(seqtab_nochimera, "seqtab_nochimera_removed_from_each_batch_20240105.rds")
seqtab_nochimera <- readRDS("seqtab_nochimera_removed_from_each_batch_20240105.rds")

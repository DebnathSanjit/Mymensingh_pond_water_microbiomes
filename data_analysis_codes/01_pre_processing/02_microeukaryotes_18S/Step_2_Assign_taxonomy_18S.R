# 31/03/2024
# Sanjit Debnath
# This script is to assign taxonomy

# Load library
library(dada2); packageVersion("dada2")

## Setup working dictionary first----
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/pre_processing_18S/primer_removed")

# Load sequence table----
seqtab_nochimera_B1 <- readRDS("seqtab_nochimera_B1_18S_20240330.rds")
seqtab_nochimera_B2 <- readRDS("seqtab_nochimera_B2_18S_20240401.rds")
seqtab_nochimera_B3 <- readRDS("seqtab_nochimera_B3_18S_20240330.rds")
seqtab_nochimera_B6 <- readRDS("seqtab_nochimera_B6_18S_20240330.rds")
length(seqtab_nochimera_B1) # 1354945
length(seqtab_nochimera_B2) #678990
length(seqtab_nochimera_B3) #590499
length(seqtab_nochimera_B6) #2041301

#merge all sequence table from all sequence batch together
seqtab.nochimera_All <- mergeSequenceTables(seqtab_nochimera_B1, seqtab_nochimera_B2, seqtab_nochimera_B3, seqtab_nochimera_B6)
seqtab.nochimera_All
#Error in mergeSequenceTables(seqtab_nochimera_B1, seqtab_nochimera_B2,  : 
# Duplicated sample names detected in the sequence table row names: 
# 10_TG2, 10_TG5, 10_TG8, 10_TH2, 10_TH5, 10_TH8, 11_TA2, 11_TA5, 11_TA8, 11_TB2, 11_TB5, 11_TB8, 
# 11_TC2, 11_TC5, 11_TC8, 11_TE2, 11_TE5, 11_TE8, 11_TF2, 11_TF5, 11_TF8, 11_TG2, 11_TG5, 11_TG8, 
# 11_TH2, 11_TH5, 11_TH8
# these were repeated sequenced in batch 6, so I removed these from batch 2
# on 20240331, I made the correct sequence table with the correct sequnces.
# so the sequence table is the correct one, don't need to save the above one, just did it to keep the record

# Save as rds file                            
#saveRDS(seqtab.nochimera_All, "seqtab_nochimera_all_18S_20240331.rds")
seqtab.nochimera_All <- readRDS("seqtab_nochimera_all_18S_20240331.rds")

# Assign taxonomy: use both Silva (for identifying bact & arch seqs) and PR2

## With Silva####
# Assign taxonomy using naive Bayesian classifier against the silva latest reference database, with a minimum bootstrap
# confidence of 60% to identify bact & arch sequences
# download the silva_nr99_v138.1_train_set.fa.gz file (https://zenodo.org/record/4587955#.ZFzADXbMLcs), and place it in the directory with the fq files.
taxa.silva <- assignTaxonomy(seqtab.nochimera_All, "silva_nr99_v138.1_train_set.fa.gz", 
                       minBoot = 60, multithread=TRUE)

# 9345 taxa identified
# Save taxa as RDS object
saveRDS(taxa.silva, "taxa_18S_90-150bp_silva_20240331.rds")

# Extract standard tables in universal formats rather than R-specific
# Give seq headers more manageable names (ASV_1, ASV_2...) instead of full seq
asv_seqs <- colnames(seqtab.nochimera_All)
asv_headers <- vector(dim(seqtab.nochimera_All)[2], mode="character")

for (i in 1:dim(seqtab.nochimera_All)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "asvs_fasta_18S_90-150bp_silva_20240331.fasta")
write(asv_fasta, "asvs_fasta_18S_90-150bp_20240404.fasta") # we need to add this after constructing the phyloseq object

# asv count table:
asv_tab <- t(seqtab.nochimera_All)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "asv_tab_18S_90-150bp_silva_20240331.tsv", sep="\t", 
            quote=F, col.names=NA)

# Save as RDS file as the asv table we will need to make the phyloseq object
saveRDS(asv_tab, "asv_tab_18S_90-150bp_silva_20240331.rds")

# Bayesian tax table:
asv_tax <- taxa.silva
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "tax_table_18S_90-150bp_silva_20240331.tsv", sep="\t", 
            quote=F, col.names=NA)

# Save as RDS file as the asv tax table we will need to make the phyloseq object
saveRDS(asv_tax, "asv_tax_tab_18S_90-150bp_silva_20240331.rds")

## With PR2####
# for 18S, we need to classify using PR2
# PR2 4.12 database, formatted for DADA2
# download the latest version: https://github.com/pr2database/pr2database/releases (download:pr2_version_5.0.0_SSU_dada2.fasta.gz)
# PR2 has different taxlevels so need to specify with taxLevels = xx
# Upgrade taxonomy from 8 levels (kingdom to species) to 9 levels (domain to species) with a new level subdivision 
# https://pr2-database.org/documentation/pr2-taxonomy-9-levels/
# minBoot = 60 because this amplicon is so short
taxa.pr2 <- assignTaxonomy(seqtab.nochimera_All, "pr2_version_5.0.0_SSU_dada2.fasta.gz", 
                                  taxLevels = c("Domain","Supergroup","Division","Subdivision","Class","Order","Family","Genus","Species"),
                                  minBoot = 60, multithread = TRUE)

saveRDS(taxa.pr2, "taxonomy_18S_90_150bp_pr2_v5.0.0_20240331.rds")

# Extract standard tables in universal formats rather than R-specific
# Give seq headers more manageable names (ASV_1, ASV_2...) instead of full seq
asv_seqs <- colnames(seqtab.nochimera_All)
asv_headers <- vector(dim(seqtab.nochimera_All)[2], mode="character")

for (i in 1:dim(seqtab.nochimera_All)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# PR2 tax table:
asv_tax.pr2 <- taxa.pr2
row.names(asv_tax.pr2) <- sub(">", "", asv_headers)
write.table(asv_tax.pr2, "asv_tax_table_18S_90_150bp_pr2_v5.0.0_20240331.tsv", sep="\t", 
            quote=F, col.names=NA)

# Save as RDS file as the asv tax table we will need to make the phyloseq object
saveRDS(asv_tax.pr2, "asv_tax_table_18S_90_150bp_pr2_v5.0.0_20240331.rds")




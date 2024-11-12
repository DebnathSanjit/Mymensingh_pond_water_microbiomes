# assign taxonomy_naive Bayesian classifier method
# Load library----
library(dada2); packageVersion("dada2")

## Setup working dictionary first----
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/Pre_processing_16S")

# Load sequence table
seqtab.nochim <- readRDS("seqtab_nochimera_20240101.rds") 

# Assign taxonomy using naive Bayesian classifier against the silva latest reference database, with a minimum bootstrap
# confidence of 80%
# download the silva_nr99_v138.1_train_set.fa.gz file (https://zenodo.org/record/4587955#.ZFzADXbMLcs), and place it in the directory with the fq files.
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", 
                       minBoot = 80, multithread=TRUE)

# Save taxa as RDS object
saveRDS(taxa, "taxa_all_batches_bayesian_20240101.rds")

# Add species assignment to ASVs with 100% match to references:
# download the silva_species_assignment_v138.1.fa.gz file (https://zenodo.org/record/4587955#.ZFzADXbMLcs), and place it in the directory with the fq files
taxa <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz")

# Save taxa as RDS object
saveRDS(taxa, "taxa_species_all_batches_bayesian_20240101.rds")

# How many have a species assignment?
length(taxa[, 7]) - sum(is.na(taxa[, 7]))  #he code calculates the number of non-missing values in the seventh column of the taxa object.

# Extract standard tables in universal formats rather than R-specific
# Give seq headers more manageable names (ASV_1, ASV_2...) instead of full seq
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs_fasta_all_batches_bayesian_20240101.fasta")

# ASV count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASV_tab_all_batches_bayesian_20240101.tsv", sep="\t", quote=F, col.names=NA)

# Save as RDS file as the asv table we will need to make the phyloseq object
saveRDS(asv_tab, "ASV_tab_all_batches_bayesian_20240101.rds")

# Bayesian tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "tax_table_all_batches_bayesian_20240101.tsv", sep="\t", quote=F, col.names=NA)

# Save as RDS file as the asv tax table we will need to make the phyloseq object
saveRDS(asv_tax, "ASV_tax_tab_all_batches_bayesian_20240101.rds")

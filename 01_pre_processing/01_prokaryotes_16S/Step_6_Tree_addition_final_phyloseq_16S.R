#Loading libraries----
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

#set working dictionary
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/Pre_processing_16S")

#small function to tidy ps object after subsetting to remove any zero taxa or samples, 
#optionally transforms counts to relative abundance (JMcM)
tidyPS <- function(ps, RA = FALSE){
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  if(RA == TRUE){
    ps <- transform_sample_counts(ps, function(x){x / sum(x)})
  }
  return(ps)}

# Load tree
tree <- read_tree ("ASVs.msa.fromAsh_20240112.treefile")

# Load phyloseq object
ps <- readRDS("phyloseq_filtered_final_20240106.rds")
ps #  12304 taxa and 1000 samples

#function to root the tree 
pick_new_outgroup <- function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape") # ape::Ntip
  # tablify parts of tree that we need.
  treeDT <-
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>%
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup) }

out_group <- pick_new_outgroup(tree)
rooted_tree <- ape::root(tree, outgroup=out_group, resolve.root=TRUE)

phy_tree(ps) <- rooted_tree
#saveRDS(ps, "phyloseq_filtered_with_tree_20240112.rds")
ps <- readRDS("phyloseq_filtered_with_tree_20240112.rds")
ps # 12304 taxa and 1000 samples


## Removing Archaea, I will not not use archaea, so remove them
ps_no_arc <- subset_taxa(ps, Kingdom !="Archaea")
ps_no_arc # 12029 taxa and 1000 samples

# Except pangasius and tilapia, I will not use other samples, so remove them
ps.pt <- ps_no_arc %>% 
  ps_filter(
    Sampling_point != 0,
    Amplicon_batch != 4,
    Pond_name != "PX"
  )
ps.pt # 10523 taxa and 892 samples

# I needed to add few columns in the metadata and with coding, I found it difficult. So, have updated the original metadata with the required columns 
# Load the new metadata
metadata.new <- read.csv("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/Metadata/Fish_pond_metadata_6_v3_20240419.csv")
# assign row names to match the sample names
rownames(metadata.new) <- metadata.new$Sample
# add the new metadata into the phyloseq object
sample_data(ps.pt) <- metadata.new
ps.pt # 10523 taxa and 891 samples

# Save the new phyloseq object and this is the final phyloseq object used for all downstream analysis
saveRDS(ps.pt, "phyloseq_metadata_6_v3_20240419.rds")



# 06/04/2024
# Sanjit Debnath
# After constructing tree, add to the phyloseq object

#Loading libraries----
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

#set working dictionary
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/Pre_processing_18S/primer_removed")

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
tree <- read_tree ("18S_tree_ASVs.msa_fromAsh_20240412.treefile")

# Load phyloseq object
ps <- readRDS("phyloseq_18S_all_filtered_pr2_90-150bp_20240404.rds")
ps #  5390 taxa and 872 samples

#function to root the tree (code from Jamie)
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

# Save as rds file
saveRDS(ps, "phyloseq_18S_filtered_with_tree_pr2_90-150bp_20240416.rds")
ps <- readRDS("phyloseq_18S_filtered_with_tree_pr2_90-150bp_20240416.rds")
ps # 5390 taxa and 872 samples 

# This is the final phyloseq object used for downstream analysis

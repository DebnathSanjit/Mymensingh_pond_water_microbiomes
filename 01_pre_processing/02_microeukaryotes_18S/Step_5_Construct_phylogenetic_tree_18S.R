# 06/04/2024
# Sanjit Debnath
# This script is to construct the phylogenetic tree. This can't be performed in normal computer, it need high tec server, Ash has access to that. He helped me to do this
#Loading libraries
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
#BiocManager::install("msa")
library(msa); packageVersion("msa")
library(phangorn); packageVersion("phangorn") # to construct a phylogenetic tree

#set working dictionary
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/pre_processing_18S/primer_removed")
getwd ()

seqtab <- readRDS("seqtab_nochimera_all_18S_20240331.rds")
#seqtab

# Construct the phylogenetic tree----
## perform a multiple-alignment of the inferred sequences----
seqs <- getSequences(seqtab) # requires dada2 package
names(seqs) <- seqs # This propagates to the tip labels of the tree
mult <- msa(seqs, method="ClustalW", type="dna", order="input") # took several hours to complete

# save this for later 
write(mult, "msa_18S_90-150bp_20240405.fasta")
saveRDS(mult, "msa_18S_90-150bp_20240405.rds")
mult <- readRDS("msa_18S_90-150bp_20240405.rds")

## first construct a neighbor-joining tree----
phang.align <- as.phyDat(mult, type="DNA", names=getSequence(seqtab))
dm <- dist.ml(phang.align) #07:27 - 07:45
treeNJ <- NJ(dm) # Note, tip order != sequence order # 07:45 - 08:10
fit = pml(treeNJ, data=phang.align)

## fit a GTR+G+I maximum likelihood tree using the neighbor-joining tree----
## negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0)) # after running for 26 hours, i stopped it,
detach("package:phangorn", unload=TRUE)





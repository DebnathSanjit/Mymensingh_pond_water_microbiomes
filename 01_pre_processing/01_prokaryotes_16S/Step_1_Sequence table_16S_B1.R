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
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/Sequences/16S")
getwd ()

## Show the path where the sequencing files are----
path <- "B1"

# Check if it has all the files you need
list.files (path)

# First need to set up a few variables that need to be used later on
# one holding the file names of all the forward reads
fnFs <- sort(list.files(path, pattern="_r1.fq.gz", full.names = TRUE))
# and one with the reverse
fnRs <- sort(list.files(path, pattern="_r2.fq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fq
# First, remove all characters after (and including) "_r"
sample.names <- sapply(strsplit(basename(fnFs), "_r"), `[`, 1)
# Then remove the sequencing run number "2689_" and anything preceding it
sample.names <- sapply(strsplit(basename(sample.names), "2689_"), `[`, 2)

# Write th sample names as text file (optional)
#write(sample.names, file="C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/Pre_processing/B1_sample_names.txt")

# view sample names
#fnFs <- fnFs
#fnRs <- fnRs
#sample.names <- sample.names
#sample.names

# Inspect read quality profiles (don't forget to run the dada2 package before starting ploting)
## Inspect forward read quality ----
#we can generate plots for all samples or for a subset of them
plotQualityProfile(fnFs[c(1, 35, 75, 90)])

## Inspect reverse read quality ----
plotQualityProfile(fnRs[c(1, 35, 75, 90)])

# a quality score of 40 and a quality score of 20 is an expected error rate of 1 in 10,000 vs 1 in 100.
# Truncating at 30 is ideal. We shouldn't truncate less than 30.
# For the forward sequences, we can truncate somewhat near 30 (~220 bp) and for reverse sequences
# we can truncate somewhat near 30 (~110 bp) where the quality distribution crashes.
# Although it seems we are removing most of our sequences, but it will cover the complete V4 with higher quality bases

## Quality Filter and trim----
# Place filtered files in filtered/ subdirectory
# Assign the filenames for the filtered fq.gz files.
filtFs <- file.path(path, "filtered", paste0(sample.names, "_f_filt.fq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_r_filt.fq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
#sample.names
#filtFs
#filtRs
#
#any(duplicated(c(fnFs, fnRs)))
#any(duplicated(c(filtFs, filtRs)))
#
#length(fnFs)
#length(fnRs)

# Quality filtering is done by 'filterAndTrim()' function
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,110), #truncates reads at this defined base
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,    # maxN = maximum number of ambiguous nucleotides (reads with 1 or more ambiguous necleotides will be removed)
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE, on Mac set multithread=TRUE
# maxEE = maximum number of estimated errors allowed in our reads, c(2,2) means F and R reads with more than 2 errors will be discarded
# truncQ = truncates the read at the first nucleotide with a specific quality score. truncQ = 2 means that the probability of the base being incorrect is 63%
head(out)

#Now let’s take a look at our filtered reads:
plotQualityProfile(filtFs[c(1, 35, 75, 90)])
plotQualityProfile(filtRs[c(1, 35, 75, 90)])

## Learn the Error Rates----
errF <- learnErrors(filtFs, multithread=TRUE)
# 100950960 total bases in 458868 reads from 24 samples will be used for learning the error rates.

errR <- learnErrors(filtRs, multithread=TRUE)
# 100666500 total bases in 915150 reads from 47 samples will be used for learning the error rates.

# visualize the estimated error rates of forward and reverse reeds
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
#The forward and reverse didn’t look too different
#The red line is what is expected based on the quality score, 
#the black line represents the estimate, and the black dots represent the observed.
#you want the observed (black dots) to track well with the estimated (black line).

## Dereplicate
# derepFastq maintains a summary of the quality information for each 
#derepF2 <- derepFastq(filtFs, verbose=TRUE)
#derepR2 <- derepFastq(filtRs, verbose=TRUE)
#
## The sample names in these objects are initially the file names of the samples, 
## This sets them to the sample names for the rest of the workflow
#names(derepF2) <- sample.names
#names(derepR2) <- sample.names
# The dereplication step is no longer listed as part of the standard dada2 tutorial, as it is performed by the dada() step if that command is given filenames.

## Inferring ASVs----
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
#Sample 1 - 14753 reads in 5454 unique sequences.
#Sample 2 - 16119 reads in 6016 unique sequences.
#Sample 3 - 18575 reads in 7262 unique sequences.

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
#Sample 1 - 14753 reads in 6242 unique sequences.
#Sample 2 - 16119 reads in 7276 unique sequences.
#Sample 3 - 18575 reads in 8715 unique sequences.

#dadaFs2 <- dada(derepF2, err=errF, multithread=TRUE, pool=TRUE)

#dadaRs2 <- dada(derepR2, err=errR, multithread=TRUE, pool=TRUE)

# Inspect the returned dada-class object:
dadaFs[[1]]
dadaRs[[1]]

## Merging forward and reverse reads/Merge paired reads----
merged_amplicons <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 45, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(merged_amplicons[[1]])

## Generating a count table/Construct sequence table----
seqtab_B1 <- makeSequenceTable(merged_amplicons)
dim(seqtab_B1)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_B1)))

# Save sequence table as RDS file
#saveRDS(seqtab_B1, "C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/Pre_processing/seqtab_B1_20240101.rds")
#seqtab_B1 <- readRDS("seqtab_B1_20240101.rds")

# The approximate size for these amplicons should be 254 bp. Amplicons of 
# sizes other than this most likely to be non-specific, so need to remove them
# Remove non-target-length sequences from the sequence table and keep only those ranges 250-256 bp
seqtab2 <- seqtab_B1[,nchar(colnames(seqtab_B1)) %in% 250:256]
#seqtab2

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab2)))

## Chimera identification/Remove chimeras----
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab2)
#  0.9896911 # in this case we lost just over 1% in terms of abundance

#sample.names <- row.names(seqtab2)  # need to check if we need this line or not?

# save the sequence table
#saveRDS(seqtab.nochim, "C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/Pre_processing/seqtab_nochimera_B1_20240101.rds")

## Track reads before, during, after processing----
#Overview of counts throughout
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merged_amplicons, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
#write.csv(track, "pre_processing_summary_B1_20240101.csv")
head(track)
#      input filtered denoisedF denoisedR merged nonchim
#1_PA1 15969    14753     13838     14235  11875   11537
#1_PB1 17501    16119     15083     15715  13185   13011
#1_PC1 20142    18575     17104     17965  14751   14666
#1_PD1 19550    18190     17144     17734  15222   14766
#1_TA1 27092    25084     24264     24657  22110   21937
#1_TD1 19447    17838     16918     17385  14855   14735

tail(track)
#          input filtered denoisedF denoisedR merged nonchim
#6_TH2       13875    12885     12010     12419  10010    9933
#PA_bl       12813    11750     11038     11426  10081    9936
#PB_bl       15342    14088     13673     13875  12855   12310
#PC_bl       12623    11331     11066     11180  10817   10630
#pcrBLANK1-1    13        7         1         1      1       1
#PD_bl       12984    11722     11362     11532  10528    9424


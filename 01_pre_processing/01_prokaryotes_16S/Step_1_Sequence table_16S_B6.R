# 02/01/2024
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
path <- "B6"

# Check if it has all the files you need
list.files (path)

# First need to set up a few variables that need to be used later on
# one holding the file names of all the forward reads
fnFs <- sort(list.files(path, pattern="_r1.fq.gz", full.names = TRUE))
# and one with the reverse
fnRs <- sort(list.files(path, pattern="_r2.fq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fq
# Create another variable with all sample names
sample.names <- sapply(strsplit(basename(fnFs), "_r"), `[`, 1)
sample.names <- sapply(strsplit(basename(sample.names), "3086_"), `[`, 2)

# Write th sample names as text file (optional)
#write(sample.names, file="C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/Pre_processing/B6_sample_names.txt")

# Assign file names and sample names to variables
#fnFs <- fnFs
#fnRs <- fnRs
#sample.names <- sample.names
#sample.names

# Inspect read quality profiles (don't forget to run the dada2 package before starting ploting)
## Inspect forward read quality ----
#we can generate plots for all samples or for a subset of them
plotQualityProfile(fnFs[c(1, 30, 50, 85, 100, 120)])

## Inspect reverse read quality ----
plotQualityProfile(fnRs[c(1, 30, 50, 85, 100, 120)])

# a quality score of 40 and a quality score of 20 is an expected error rate of 1 in 10,000 vs 1 in 100.
# Truncating at 30 is ideal. We shouldn't truncate less than 30.
# For the forward sequences, we can truncate somewhat near 30 (~210 bp) and for reverse sequences
# we can truncate somewhat near 30 (~190 bp) where the quality distribution crashes.
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
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(210, 190), #truncates reads at this defined base
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,    # maxN = maximum number of ambiguous nucleotides (reads with 1 or more ambiguous necleotides will be removed)
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE, on Mac set multithread=TRUE
# maxEE = maximum number of estimated errors allowed in our reads, c(2,2) means F and R reads with more than 2 errors will be discarded
# truncQ = truncates the read at the first nucleotide with a specific quality score. truncQ = 2 means that the probability of the base being incorrect is 63%
head(out)

#Now let’s take a look at our filtered reads:
plotQualityProfile(filtFs[c(1, 30, 50, 85, 100, 120)])
plotQualityProfile(filtRs[c(1, 30, 50, 85, 100, 120)])

## Learn the Error Rates----
errF <- learnErrors(filtFs, multithread=TRUE)
# 105921270 total bases in 504387 reads from 26 samples will be used for learning the error rates.

errR <- learnErrors(filtRs, multithread=TRUE)
# 103119840 total bases in 542736 reads from 28 samples will be used for learning the error rates.

# visualize the estimated error rates of forward and reverse reeds
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
#The forward and reverse didn’t look too different
#The red line is what is expected based on the quality score, 
#the black line represents the estimate, and the black dots represent the observed.
#you want the observed (black dots) to track well with the estimated (black line).

## Inferring ASVs----
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
#Sample 1 - 12142 reads in 4613 unique sequences.
#Sample 2 - 22641 reads in 7663 unique sequences.
#Sample 3 - 25886 reads in 8836 unique sequences.

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
#Sample 1 - 12142 reads in 4687 unique sequences.
#Sample 2 - 22641 reads in 7955 unique sequences.
#Sample 3 - 25886 reads in 9492 unique sequences.

# Inspect the returned dada-class object:
dadaFs[[1]]
dadaRs[[1]]

## Merging forward and reverse reads/Merge paired reads----
merged_amplicons <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 45, verbose=TRUE)

#while i truncate the reverse sequence at 180, merging was not working. so I truncate at 120 and now it's working
# Inspect the merger data.frame from the first sample
head(merged_amplicons[[1]])

## Generating a count table/Construct sequence table----
seqtab_B6 <- makeSequenceTable(merged_amplicons)
#dim(seqtab_B6)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_B6)))

# Save sequence table as RDS file
#saveRDS(seqtab_B6, "C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/Pre_processing/seqtab_B6_20231231.rds")

# The approximate size for these amplicons should be 254 bp. Amplicons of 
# sizes other than this most likely to be non-specific, so need to remove them
# Remove non-target-length sequences from the sequence table and keep only those ranges 250-256 bp
seqtab2 <- seqtab_B6[,nchar(colnames(seqtab_B6)) %in% 250:256]
#seqtab2

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab2)))

## Chimera identification/Remove chimeras----
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab2)
#  0.996287 # in this case we lost less than 1% in terms of abundance

#sample.names <- row.names(seqtab2)  # need to check if we need this line or not?

# save the sequence table
#saveRDS(seqtab.nochim, "C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/Pre_processing/seqtab_nochimera_B6_20240102.rds")

## Track reads before, during, after processing----
#Overview of counts throughout
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merged_amplicons, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
#write.csv(track, "pre_processing_summary_B6_20240102.csv")

head(track)
#        input filtered denoisedF denoisedR merged nonchim
#15_PA5 13282    12142     11449     11565  10130   10118
#15_PB5 24736    22641     21857     21851  20190   20146
#15_PC5 29029    25886     24987     25002  23146   23088
#15_PD5 40795    36825     35610     35633  32892   32362
#15_PE5 38917    35253     34027     34108  31585   31533
#15_PF5 14388    13045     12352     12326  10904   10902

tail(track)
#    input filtered denoisedF denoisedR merged nonchim
#T1a 20557    18427     17767     17657  16319   16292
#T1b 28535    25372     24595     24594  22983   22920
#T1c 27201    24358     23460     23569  21893   21875
#T2a 19846    18048     17589     17589  16627   16595
#T2b 26723    24149     23589     23530  22354   22319
#T2c 23085    20959     20441     20421  19406   19328


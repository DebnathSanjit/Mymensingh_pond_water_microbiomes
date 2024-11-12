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
path <- ("B7")

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
sample.names <- sapply(strsplit(basename(sample.names), "10017_"), `[`, 2)

# Write th sample names as text file (optional)
#write(sample.names, file="C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/Pre_processing/B7_sample_names.txt")

#fnFs <- fnFs
#fnRs <- fnRs
#sample.names <- sample.names
#sample.names

# Inspect read quality profiles (don't forget to run the dada2 package before starting ploting)
## Inspect forward read quality ----
#we can generate plots for all samples or for a subset of them
plotQualityProfile(fnFs[c(1, 40, 80, 120, 150, 170)])

## Inspect reverse read quality ----
plotQualityProfile(fnRs[c(1, 40, 80, 120, 150, 170)])

# a quality score of 40 and a quality score of 20 is an expected error rate of 1 in 10,000 vs 1 in 100.
# Truncating at 30 is ideal. We shouldn't truncate less than 30.
# For the forward sequences, we can truncate somewhat near 30 (~190 bp) and for reverse sequences
# we can truncate somewhat near 30 (~180 bp) where the quality distribution crashes.
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
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(190, 180), #truncates reads at this defined base
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,    # maxN = maximum number of ambiguous nucleotides (reads with 1 or more ambiguous necleotides will be removed)
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE, on Mac set multithread=TRUE
# maxEE = maximum number of estimated errors allowed in our reads, c(2,2) means F and R reads with more than 2 errors will be discarded
# truncQ = truncates the read at the first nucleotide with a specific quality score. truncQ = 2 means that the probability of the base being incorrect is 63%
head(out)
# 10017_20_TC2_rep_r2.fq this sequence has some problem, so this pair is not filtering, so has to remove from this batch

#Now let’s take a look at our filtered reads:
plotQualityProfile(filtFs[c(1, 40, 80, 120, 150, 170)])
plotQualityProfile(filtRs[c(1, 40, 80, 120, 150, 170)])

## Learn the Error Rates----
errF <- learnErrors(filtFs, multithread=TRUE)
#103907960 total bases in 546884 reads from 38 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, multithread=TRUE)
# 101734200 total bases in 565190 reads from 39 samples will be used for learning the error rates.

# visualize the estimated error rates of forward and reverse reeds
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
#The forward and reverse didn’t look too different
#The red line is what is expected based on the quality score, 
#the black line represents the estimate, and the black dots represent the observed.
#you want the observed (black dots) to track well with the estimated (black line).

## Inferring ASVs----
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
#Sample 1 - 14851 reads in 5095 unique sequences.
#Sample 2 - 5896 reads in 2176 unique sequences.
#Sample 3 - 8897 reads in 3429 unique sequences.

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
#Sample 1 - 14851 reads in 5119 unique sequences.
#Sample 2 - 5896 reads in 1982 unique sequences.
#Sample 3 - 8897 reads in 3965 unique sequences.

# Inspect the returned dada-class object:
dadaFs[[1]]
dadaRs[[1]]

## Merging forward and reverse reads/Merge paired reads----
merged_amplicons <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 45, verbose=TRUE)
#12452 paired-reads (in 423 unique pairings) successfully merged out of 13705 (in 718 pairings) input.
#5025 paired-reads (in 218 unique pairings) successfully merged out of 5365 (in 305 pairings) input.
#6906 paired-reads (in 336 unique pairings) successfully merged out of 7875 (in 509 pairings) input.

#while i truncate the reverse sequence at 180, merging was not working. so I truncate at 120 and now it's working
# Inspect the merger data.frame from the first sample
head(merged_amplicons[[1]])

## Generating a count table/Construct sequence table----
seqtab_B7 <- makeSequenceTable(merged_amplicons)
#dim(seqtab_B7)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_B7)))

# Save sequence table as RDS file
#saveRDS(seqtab_B7, "C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/Pre_processing/seqtab_B7_20240101.rds")

# The approximate size for these amplicons should be 254 bp. Amplicons of 
# sizes other than this most likely to be non-specific, so need to remove them
# Remove non-target-length sequences from the sequence table and keep only those ranges 250-256 bp
seqtab2 <- seqtab_B7[,nchar(colnames(seqtab_B7)) %in% 250:256]
#seqtab2

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab2)))

## Chimera identification/Remove chimeras----
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab2)
#  0.9910031 # in this case we lost less than 1% in terms of abundance

#sample.names <- row.names(seqtab2)  # need to check if we need this line or not?

# save the sequence table
#saveRDS(seqtab.nochim, "C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/Pre_processing/seqtab_nochimera_B7_20240102.rds")

## Track reads before, during, after processing----
#Overview of counts throughout
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merged_amplicons, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
#write.csv(track, "pre_processing_summary_B7_20240101.csv")

head(track)
#        input filtered denoisedF denoisedR merged nonchim
#10_PA8 15772    14851     14069     14136  12452   12383
#10_PB8  6352     5896      5513      5571   5025    4990
#10_PC8  9533     8897      8212      8173   6906    6906
#10_PD8 16106    14995     14449     14373  12829   12817
#10_PE8 16437    15581     14826     14846  13337   13244
#10_PF8 14132    13341     12525     12530  10979   10967

tail(track)
#          input filtered denoisedF denoisedR merged nonchim
#mockDNA1_2 50281    45816     45673     45687  41648   38339
#mockDNA3_1 44688    40737     40617     40530  37012   33912
#pcr_neg2_1    78       63        42        43     41      41
#pcr_neg2_2    12       10         1         1      0       0
#pcr_neg4_2    12        8         3         3      3       3
#pcr_neg6_3  1101      980       962       963    957     955




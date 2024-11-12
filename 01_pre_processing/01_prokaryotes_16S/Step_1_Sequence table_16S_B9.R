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
path <- ("B9")

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
#write(sample.names, file="C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/B9_sample_names.txt")

#instead of subsetting, I will try all sequence for Pangasisus
#fnFs <- fnFs
#fnRs <- fnRs
#sample.names <- sample.names
#sample.names

# Inspect read quality profiles (don't forget to run the dada2 package before starting ploting)
## Inspect forward read quality ----
#we can generate plots for all samples or for a subset of them
plotQualityProfile(fnFs[c(1, 20, 70, 100, 223, 228)])

## Inspect reverse read quality ----
plotQualityProfile(fnRs[c(1, 20, 70, 100, 223, 228)])

#a quality score of 40 and a quality score of 20 is an expected error rate of 1 in 10,000 vs 1 in 100.
# For the forward sequences, we can truncate somewhat near 40 (~200 bp) but for reverse sequences
# we can truncate somewhat near 30 (~180 bp) where the quality distribution crashes.
# But after truncating at 180 for the reverse, mergepair function didn't work. so I truncate at 120
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
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230, 160), #truncates reads at this defined base
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,    # maxN = maximum number of ambiguous nucleotides (reads with 1 or more ambiguous necleotides will be removed)
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE, on Mac set multithread=TRUE
# maxEE = maximum number of estimated errors allowed in our reads, c(2,2) means F and R reads with more than 2 errors will be discarded
# truncQ = truncates the read at the first nucleotide with a specific quality score. truncQ = 2 means that the probability of the base being incorrect is 63%
head(out)

#Now let’s take a look at our filtered reads:
plotQualityProfile(filtFs[c(1, 20, 70, 100, 223, 228)])
plotQualityProfile(filtRs[c(1, 20, 70, 100, 223, 228)])

## Learn the Error Rates----
errF <- learnErrors(filtFs, multithread=TRUE)
# 102778030 total bases in 446861 reads from 25 samples will be used for learning the error rates.

errR <- learnErrors(filtRs, multithread=TRUE)
# 108598080 total bases in 678738 reads from 38 samples will be used for learning the error rates.

# visualize the estimated error rates of forward and reverse reeds
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
#The forward and reverse didn’t look too different
#The red line is what is expected based on the quality score, 
#the black line represents the estimate, and the black dots represent the observed.
#you want the observed (black dots) to track well with the estimated (black line).

## Inferring ASVs----
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
#Sample 1 - 9500 reads in 3348 unique sequences.
#Sample 2 - 5364 reads in 1998 unique sequences.
#Sample 3 - 6630 reads in 2760 unique sequences.

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
#Sample 1 - 9500 reads in 3158 unique sequences.
#Sample 2 - 5364 reads in 1974 unique sequences.
#Sample 3 - 6630 reads in 2756 unique sequences.

# Inspect the returned dada-class object:
dadaFs[[1]]
dadaRs[[1]]

## Merging forward and reverse reads/Merge paired reads----
merged_amplicons <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 45, verbose=TRUE)
#7746 paired-reads (in 290 unique pairings) successfully merged out of 8661 (in 472 pairings) input.
#4237 paired-reads (in 206 unique pairings) successfully merged out of 4765 (in 325 pairings) input.
#4963 paired-reads (in 251 unique pairings) successfully merged out of 5753 (in 402 pairings) input.

#while i truncate the reverse sequence at 180, merging was not working. so I truncate at 120 and now it's working
# Inspect the merger data.frame from the first sample
head(merged_amplicons[[1]])

## Generating a count table/Construct sequence table----
seqtab_B9 <- makeSequenceTable(merged_amplicons)
#dim(seqtab_B9)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_B9)))

# Save sequence table as RDS file
#saveRDS(seqtab_B9, "C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/Pre_processing/seqtab_B9_20240104.rds")

# The approximate size for these amplicons should be 254 bp. Amplicons of 
# sizes other than this most likely to be non-specific, so need to remove them
# Remove non-target-length sequences from the sequence table and keep only those ranges 250-256 bp
seqtab2 <- seqtab_B9[,nchar(colnames(seqtab_B9)) %in% 250:256]
#seqtab2

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab2)))

## Chimera identification/Remove chimeras----
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
# Identified 1072 bimeras out of 14677 input sequences.
#dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab2)
#  0.9931946 # in this case we lost less than 1% in terms of abundance

#sample.names <- row.names(seqtab2)  # need to check if we need this line or not?

# save the sequence table
#saveRDS(seqtab.nochim, "C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/Pre_processing/seqtab_nochimera_B9_20240109.rds")

## Track reads before, during, after processing----
#Overview of counts throughout
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merged_amplicons, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
#write.csv(track, "pre_processing_summary_B9_20240101.csv")

head(track)
#      input filtered denoisedF denoisedR merged nonchim
#1_PA2 10188     9500      8894      9014   7746    7656
#1_PB2  5829     5364      4938      5025   4237    4230
#1_PC2  7340     6630      5985      6113   4963    4930
#1_PD2 14000    12889     12278     12361  11162   10752
#1_TA2 13547    12185     11660     11706  10371   10321
#1_TD2  7965     7107      6580      6656   5636    5631

tail(track)
#          input filtered denoisedF denoisedR merged nonchim
#mockDNA2_2 35984    30458     30346     30347  28743   26769
#pcr_neg3_1    61       53        35        33     28      28
#pcr_neg3_2   223      178       165       159    145     145
#pcr_neg3_3    13        6         1         1      0       0
#pcr_neg3_4    44       30        19        19     19      19
#pcr_neg6_2    79       55        29        33     27      27
#

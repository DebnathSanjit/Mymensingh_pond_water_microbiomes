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
path <- "B2"

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
sample.names <- sapply(strsplit(basename(sample.names), "2716_"), `[`, 2)

# Write th sample names as text file (optional)
#write(sample.names, file="C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/Pre_processing/B2_sample_names.txt")

#fnFs <- fnFs
#fnRs <- fnRs
#sample.names <- sample.names
#sample.names

# Inspect read quality profiles (don't forget to run the dada2 package before starting ploting)
## Inspect forward read quality ----
#we can generate plots for all samples or for a subset of them
plotQualityProfile(fnFs[c(1, 25, 50, 70)])

## Inspect reverse read quality ----
plotQualityProfile(fnRs[c(1, 25, 50, 70)])
# a quality score of 40 and a quality score of 20 is an expected error rate of 1 in 10,000 vs 1 in 100.
# Truncating at 30 is ideal. We shouldn't truncate less than 30.
# For the forward sequences, we can truncate somewhat near 30 (~240 bp) and for reverse sequences
# we can truncate somewhat near 30 (~200 bp) where the quality distribution crashes.
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
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,200), #truncates reads at this defined base
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,    # maxN = maximum number of ambiguous nucleotides (reads with 1 or more ambiguous necleotides will be removed)
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE, on Mac set multithread=TRUE
# maxEE = maximum number of estimated errors allowed in our reads, c(2,2) means F and R reads with more than 2 errors will be discarded
# truncQ = truncates the read at the first nucleotide with a specific quality score. truncQ = 2 means that the probability of the base being incorrect is 63%
head(out)
#                      reads.in reads.out
#2716_10_PA5_r1.fq.gz    20564     18488
#2716_10_PB5_r1.fq.gz     9796      8970
#2716_10_PC5_r1.fq.gz    10909      9882
#2716_10_PD5_r1.fq.gz    23166     20916
#2716_10_PE5_r1.fq.gz    20380     18492
#2716_10_PF5_r1.fq.gz    17559     16168

#Now let’s take a look at our filtered reads:
plotQualityProfile(filtFs[c(1, 25, 50, 70)])
plotQualityProfile(filtRs[c(1, 25, 50, 70)])

## Learn the Error Rates----
errF <- learnErrors(filtFs, multithread=TRUE)
# 100949520 total bases in 420623 reads from 28 samples will be used for learning the error rates.

errR <- learnErrors(filtRs, multithread=TRUE)
# 104638600 total bases in 523193 reads from 32 samples will be used for learning the error rates.

# visualize the estimated error rates of forward and reverse reeds
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
#The forward and reverse didn’t look too different
#The red line is what is expected based on the quality score, 
#the black line represents the estimate, and the black dots represent the observed.
#you want the observed (black dots) to track well with the estimated (black line).

## Inferring ASVs----
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
#Sample 1 - 18488 reads in 6677 unique sequences.
#Sample 2 - 8970 reads in 2714 unique sequences.
#Sample 3 - 9882 reads in 3402 unique sequences.

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
#Sample 1 - 18488 reads in 6597 unique sequences.
#Sample 2 - 8970 reads in 2696 unique sequences.
#Sample 3 - 9882 reads in 3628 unique sequences.

# Inspect the returned dada-class object:
dadaFs[[1]]
dadaRs[[1]]

## Merging forward and reverse reads/Merge paired reads----
merged_amplicons <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 45, verbose=TRUE)

#while i truncate the reverse sequence at 180, merging was not working. so I truncate at 120 and now it's working
# Inspect the merger data.frame from the first sample
head(merged_amplicons[[1]])

## Generating a count table/Construct sequence table----
seqtab_B2 <- makeSequenceTable(merged_amplicons)
dim(seqtab_B2)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_B2)))

# Save sequence table as RDS file
#saveRDS(seqtab_B2, "C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/Pre_processing/seqtab_B2_20231231.rds")

# The approximate size for these amplicons should be 254 bp. Amplicons of 
# sizes other than this most likely to be non-specific, so need to remove them
# Remove non-target-length sequences from the sequence table and keep only those ranges 250-256 bp
seqtab2 <- seqtab_B2[,nchar(colnames(seqtab_B2)) %in% 250:256]
seqtab2

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab2)))

## Chimera identification/Remove chimeras----
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab2)
#  0.9947242 # in this case we lost less than 1% in terms of abundance

sample.names <- row.names(seqtab2)  # need to check if we need this line or not?

# save the sequence table
#saveRDS(seqtab.nochim, "C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/Pre_processing/seqtab_nochimera_B2_20240102.rds")

## Track reads before, during, after processing----
#Overview of counts throughout
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merged_amplicons, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
#write.csv(track, "pre_processing_summary_B2_20240102.csv")

head(track)
#        input filtered denoisedF denoisedR merged nonchim
#10_PA5 20564    18488     17360     17349  15838   15652
#10_PB5  9796     8970      8546      8513   7915    7887
#10_PC5 10909     9882      9232      9204   8327    8311
#10_PD5 23166    20916     20152     20217  19148   19095
#10_PE5 20380    18492     17668     17746  16712   16611
#10_PF5 17559    16168     15178     15171  13933   13877

tail(track)
#            input filtered denoisedF denoisedR merged nonchim
#9_TD8       12598    11406     10903     10844  10122   10083
#9_TE5        8676     7758      7304      7203   6463    6440
#9_TF5       15820    13161     12479     12419  11003   11000
#9_TG5       17116    15084     14245     14250  13017   12958
#9_TH5       21995    19972     19186     19240  17715   17678
#extr_neg2-1   128      108        78        80     78      74

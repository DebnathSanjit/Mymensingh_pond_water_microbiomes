# 01/04/2024
# Sanjit Debnath

#Processing with DADA2 in R
# Removing primers----
# https://benjjneb.github.io/dada2/ITS_workflow.html

## Load dada2 package----
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")

## Setup working dictionary first----
setwd("C:/Users/Scd_Research_Data/Data/From_Dominique/My_interest/All_species/Sequences")
getwd ()

# batch 2----
## Show the path where the sequencing files are----
path <- paste0(getwd(), "/18S/B2")

# Check if it has all the files you need
list.files (path)

# First need to set up a few variables that need to be used later on
# one holding the file names of all the forward reads
fnFs <- sort(list.files(path, pattern="_r1.fq.gz", full.names = TRUE))
# and one with the reverse
fnRs <- sort(list.files(path, pattern="_r2.fq.gz", full.names = TRUE))

# Identify primers
FWD <- "TATCGCCGTTCGGTACACACCGCCCGTC"
REV <- "AGTCAGTCAGCATGATCCTTCTGCAGGTTCACCTAC"

# to verify the presence and orientation of these primers in the data
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), 
               Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

#“pre-filter” the sequences just to remove those with Ns, but perform no other filtering.
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = FALSE)

#We are now ready to count the number of times the primers appear in the forward and reverse read
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#                    Forward Complement Reverse RevComp
#FWD.ForwardReads       0          0       0       0
#FWD.ReverseReads       0          0       0       4475
#REV.ForwardReads       0          0       0       1
#REV.ReverseReads       0          0       0       0

#Remove Primers----
#Install cutadapat. To install cutadapt, first need to install "Anaconda Promt (Miniconda3)" and follow the code below
# Open the command line (cmd.exe) and run py -m pip install cutadapt.
#Test whether it worked by running py -m cutadapt --version. You should see the version number of Cutadapt.
# https://cutadapt.readthedocs.io/en/stable/installation.html#installation-on-windows
#After installing cutadapt, we need to tell R the path to the cutadapt command.
#set the cutadapt path (https://github.com/benjjneb/dada2/issues/1874)
cutadapt <- "C:/Users/scd226/AppData/Local/Programs/Python/Python312/Scripts/cutadapt.exe"
system2(cutadapt, args = "--version") # Run shell commands from R
# 4.6 # if R found cutadapt, it will show the cutadapt version

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# As a sanity check, we will count the presence of primers in the first cutadapt-ed sample
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

#This is cutadapt 4.6 with Python 3.12.1
#Command line parameters: -g TATCGCCGTTCGGTACACACCGCCCGTC -a GTAGGTGAACCTGCAGAAGGATCATGCTGACTGACT -G AGTCAGTCAGCATGATCCTTCTGCAGGTTCACCTAC -A GACGGGCGGTGTGTACCGAACGGCGATA -n 2 -o C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/Sequences/18S/B1/test/cutadapt/2715_1_PA1_r1.fq.gz -p C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/Sequences/18S/B1/test/cutadapt/2715_1_PA1_r2.fq.gz C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/Sequences/18S/B1/test/filtN/2715_1_PA1_r1.fq.gz C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/Sequences/18S/B1/test/filtN/2715_1_PA1_r2.fq.gz
#Run "cutadapt --help" to see command-line options.
#See https://cutadapt.readthedocs.io/ for full documentation.
#
#cutadapt: error: unrecognized arguments: - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/Sequences/18S/B1/test/cutadapt/2715_1_PA1_r2.fq.gz C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/Sequences/18S/B1/test/filtN/2715_1_PA1_r1.fq.gz C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/Sequences/18S/B1/test/filtN/2715_1_PA1_r2.fq.gz
#This is cutadapt 4.6 with Python 3.12.1
# to solve this use, need to use such working dictionary where it doesn't contain any space or parentheses. 
# see this: https://github.com/benjjneb/dada2/issues/984
#                    Forward Complement Reverse RevComp
#FWD.ForwardReads       0          0       0       0
#FWD.ReverseReads       0          0       0       0
#REV.ForwardReads       0          0       0       0
#REV.ReverseReads       0          0       0       0

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_r1.fq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_r2.fq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
#get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
#sample.names <- unname(sapply(cutFs, get.sample.name))
#head(sample.names)

# First need to set up a few variables that need to be used later on
# one holding the file names of all the forward reads
#fnFs <- sort(list.files(path, pattern="_r1.fq.gz", full.names = TRUE))
# and one with the reverse
#fnRs <- sort(list.files(path, pattern="_r2.fq.gz", full.names = TRUE))

## extract sample names----
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fq
# Create another variable with all sample names
sample.names <- sapply(strsplit(basename(cutFs), "_r"), `[`, 1)
# Then remove the sequencing run number "2715_" and anything preceding it
sample.names <- sapply(strsplit(basename(sample.names), "2725_"), `[`, 2)

# Write th sample names as text file (optional)
write(sample.names, file="C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/1.From_Dominique/R_scripts/pre_processing_18S/B1_sample_names_20240401.txt")

# Inspect read quality profiles (don't forget to run the dada2 package before starting ploting)
## Inspect forward read quality ----
#we can generate plots for all samples or for a subset of them
plotQualityProfile(fnFs[c(1, 65, 90, 180)]) 

## Inspect reverse read quality ----
plotQualityProfile(fnRs[c(1, 65, 90, 180)])

#a quality score of 40 and a quality score of 20 is an expected error rate of 1 in 10,000 vs 1 in 100.
# For the forward sequences, we can truncate somewhat near 40 (~200 bp) but for reverse sequences
# we can truncate somewhat near 30 (~100 bp) where the quality distribution crashes.
# Although it seems we are removing most of our sequences, but it will cover the complete V4 with higher quality bases

## Quality Filter and trim----
# Place filtered files in filtered/ subdirectory
# Assign the filenames for the filtered fq.gz files.
filtFs <- file.path(path, "18S_B2_filtered", paste0(sample.names, "_f_filt.fq.gz"))
filtRs <- file.path(path, "18S_B2_filtered", paste0(sample.names, "_r_filt.fq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
#sample.names

# Quality filtering is done by 'filterAndTrim()' function
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(130,130), # Jamie trimed at 100 bp position
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)
# On Windows set multithread=FALSE, on Mac set multithread=TRUE
head(out)

#Now let’s take a look at our filtered reads:
plotQualityProfile(filtFs[c(1, 65, 90, 180)])
plotQualityProfile(filtRs[c(1, 65, 90, 180)])

## Learn the Error Rates----
errF <- learnErrors(filtFs, multithread=TRUE)
#101282350 total bases in 779095 reads from 81 samples will be used for learning the error rates.

errR <- learnErrors(filtRs, multithread=TRUE)
# 101282350 total bases in 779095 reads from 81 samples will be used for learning the error rates.

# visualize the estimated error rates of forward and reverse reeds
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
#The forward and reverse didn’t look too different
#The red line is what is expected based on the quality score, 
#the black line represents the estimate, and the black dots represent the observed.
#you want the observed (black dots) to track well with the estimated (black line).

## Inferring ASVs----
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
#Sample 1 - 12081 reads in 3086 unique sequences.
#Sample 2 - 15433 reads in 3735 unique sequences.
#Sample 3 - 11809 reads in 2978 unique sequences.
#Sample 4 - 9470 reads in 1973 unique sequences.

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
#Sample 1 - 12081 reads in 2552 unique sequences.
#Sample 2 - 15433 reads in 3166 unique sequences.
#Sample 3 - 11809 reads in 2699 unique sequences.
#Sample 4 - 9470 reads in 1667 unique sequences.

# Inspect the returned dada-class object:
dadaFs[[1]]
dadaRs[[5]]

## Merging forward and reverse reads/Merge paired reads----
merged_amplicons <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 30, verbose=TRUE)
#11152 paired-reads (in 254 unique pairings) successfully merged out of 11547 (in 355 pairings) input.
#14204 paired-reads (in 281 unique pairings) successfully merged out of 14689 (in 409 pairings) input.
#10805 paired-reads (in 248 unique pairings) successfully merged out of 11188 (in 358 pairings) input.

# Inspect the merger data.frame from the first sample
head(merged_amplicons[[1]])

## Generating a count table/Construct sequence table----
seqtab_B2 <- makeSequenceTable(merged_amplicons)
dim(seqtab_B2)
#  195 4313

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_B2)))

# Save sequence table as RDS file
#saveRDS(seqtab_B2, "C:/Users/Scd_Research_Data/Data/From_Dominique/My_interest/All_species/R_scripts/Scripts_for_18S/primer_removed/seqtab_B2_18S_20240401.rds")
seqtab_B2 <- readRDS("C:/Users/Scd_Research_Data/Data/From_Dominique/My_interest/All_species/R_scripts/Scripts_for_18S/primer_removed/seqtab_B2_18S_20240401.rds")

# The approximate size for these amplicons should be 90-150 bp. Amplicons of 
# sizes other than this most likely to be non-specific, so need to remove them
# Remove non-target-length sequences from the sequence table and keep only those ranges 90-150 bp
seqtab2 <- seqtab_B2[,nchar(colnames(seqtab_B2)) %in% 90:150]

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab2)))
#130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 
#78 522 243 540 392 237 151 179 185 174 145 131 111  81  67  65  50  40  51  29  16

## Chimera identification/Remove chimeras----
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
# 195 3482

sum(seqtab.nochim)/sum(seqtab2)
#  0.9996529 # in this case we lost less than 1% in terms of abundance

# save the sequence table
saveRDS(seqtab.nochim, "C:/Users/Scd_Research_Data/Data/From_Dominique/My_interest/All_species/R_scripts/Scripts_for_18S/primer_removed/seqtab_nochimera_B2_18S_20240401.rds")

## Track reads before, during, after processing----
#Overview of counts throughout
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merged_amplicons, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
#        input filtered denoisedF denoisedR merged nonchim
#10_PA2 12215    12081     11748     11673  11152   10553
#10_PA5 15606    15433     14916     14938  14204   12937
#10_PA8 11956    11809     11372     11411  10805   10117
#10_PB2  9567     9470      9219      9265   8862    8455
#10_PB5 10529    10416     10270     10225   9935    9588
#10_PB8 13882    13740     13502     13487  13134   12439

tail(track)
#              input filtered denoisedF denoisedR merged nonchim
#extr_neg2-1      9        7         1         3      0       0
#mock_dna2-1   5942     5828      5813      5827   4689    3216
#mock_dna2-2   8551     8423      8417      8420   7168    6407
#mock_extr2-1  7453     7343      7341      7330   5647    5153
#pcr_neg2-1       6        3         2         2      2       2
#pcr_neg2-2       3        1         1         1      0       0

write.csv(track, "C:/Users/Scd_Research_Data/Data/From_Dominique/My_interest/All_species/R_scripts/Scripts_for_18S/primer_removed/B2_18S_dada2_summary_20240401.csv")


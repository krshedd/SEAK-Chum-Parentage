# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2014 CM 38, Project 37 Chum SNP selection for the AHRP
## 380 fish (95 from each of the 4 pedigree streams), all 188 available chum markers
## Get all 2013 adults from Fish Creek to boost sample size
## Split all hatchery fish into their own silly (5 sillys total, 4 for the streams and 1 for DIPAC)
## Kyle Shedd Wed Jul 08 15:24:10 2015
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


ls()
rm(list=ls(all=TRUE))
search()
getwd()
setwd("V:/WORK/Chum/AHRP")

## save.image("Phase1_CM38_SNPselection_newer.RData")
## load("Phase1_CM38_SNPselection_newer.RData")

#This sources all of the new GCL functions to this workspace
source("V:/DATA/R_GEN/GCL Source Scripts/Functions.GCL.R")
source("V:/WORK/Kyle/R Source Scripts/Functions.GCL_KS.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get unaltered .gcl objects for Admiralty, Prospect, and Sawmill Creek (CM 38)
CMADMCR13.gcl <- dget(file = "V:/WORK/Chum/AHRP/Raw genotypes/Phase 1_QCed/CMADMCR13.txt")
CMPROSCR13.gcl <- dget(file = "V:/WORK/Chum/AHRP/Raw genotypes/Phase 1_QCed/CMPROSCR13.txt")
CMSAWCR13.gcl <- dget(file = "V:/WORK/Chum/AHRP/Raw genotypes/Phase 1_QCed/CMSAWCR13.txt")

# Check structure of attributes table
str(CMADMCR13.gcl$attributes)
names(CMADMCR13.gcl$attributes)

# Get rid of "added" attributes (can't pool collections if attributes tables are different)
CMADMCR13.gcl$attributes <- CMADMCR13.gcl$attributes[, -c(25:29)]
CMPROSCR13.gcl$attributes <- CMPROSCR13.gcl$attributes[, -c(25:29)]
CMSAWCR13.gcl$attributes <- CMSAWCR13.gcl$attributes[, -c(25:29)]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### VIC/FAM Invstigation from Biomark.csv files ####
## Make sure that VIC/FAM are correct between CM38 and 40
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Panel A 1:96 SNPs
file.CM38.1.2.A = "V:\\DATA\\All SNP data\\Chum\\Project CM038 - PBT Phase 1\\Biomark\\Data\\CM38A\\Combined Chips\\Done\\CM38Phase1_1A-2A\\CM38Phase1_1A-2A Combined Chip Run.csv"
file.CM40.1.2.A = "V:\\DATA\\All SNP data\\Chum\\Project CM040 - PBT Phase 2\\Project CM40 BioMark\\Data\\CM40A\\Combined Chips\\Done\\CM40A 1-2 after correction\\CM40A 01-02 after correction Combined Chip Run.csv"

biomark.CM38.1.2.A <- read.csv(file = file.CM38.1.2.A, skip = 15, as.is = TRUE)
biomark.CM40.1.2.A <- read.csv(file = file.CM40.1.2.A, skip = 15, as.is = TRUE)

str(biomark.CM38.1.2.A)
str(biomark.CM40.1.2.A)

CM40.A.alleles <- biomark.CM40.1.2.A[1:96, 4:6]
CM38.A.alleles <- biomark.CM38.1.2.A[1:96, 4:6]

# Are loci in same order? YES
table(CM38.A.alleles$Assay == CM40.A.alleles$Assay)

# Are VIC/FAM switched? YES, 2 markers
table(CM38.A.alleles$Allele.Y == CM40.A.alleles$Allele.Y)

# Which markers are switched in panel A? "Oke_RPN1-80"  "Oke_ROA1-209"
SNP.issues.PanelA <- CM38.A.alleles$Assay[!CM38.A.alleles$Allele.Y == CM40.A.alleles$Allele.Y]


## Panel B 97:188 SNPs
file.CM38.1.2.B = "V:\\DATA\\All SNP data\\Chum\\Project CM038 - PBT Phase 1\\Biomark\\Data\\CM38B\\Combined Chips\\CM38Phase1_1B-2B\\CM38Phase1_1B-2B Combined Chip Run.csv"
file.CM40.1.2.B = "V:\\DATA\\All SNP data\\Chum\\Project CM040 - PBT Phase 2\\Project CM40 BioMark\\Data\\CM40B\\Combined Chips\\Done\\CM40B 5-6\\CM40B 5-6 Combined Chip Run.csv"

biomark.CM38.1.2.B <- read.csv(file = file.CM38.1.2.B, skip = 15, as.is = TRUE)
biomark.CM40.1.2.B <- read.csv(file = file.CM40.1.2.B, skip = 15, as.is = TRUE)

str(biomark.CM38.1.2.B)
str(biomark.CM40.1.2.B)

unique(biomark.CM38.1.2.B$Assay)
length(unique(biomark.CM38.1.2.B$Assay))

CM40.B.alleles <- biomark.CM40.1.2.B[1:93, 4:6]
CM38.B.alleles <- biomark.CM38.1.2.B[1:93, 4:6]

# Are loci in same order? YES
table(CM38.B.alleles$Assay == CM40.B.alleles$Assay)

# Are VIC/FAM switched? YES, 1 marker
table(CM38.B.alleles$Allele.Y == CM40.B.alleles$Allele.Y)

# Which markers are switched in panel B? "Oke_U2025-86"
SNP.issues.PanelB <- CM38.B.alleles$Assay[!CM38.B.alleles$Allele.Y == CM40.B.alleles$Allele.Y]

# Master list
SNP.issues <- c(SNP.issues.PanelA, SNP.issues.PanelB)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Flip Allele 1/2 for "Oke_ROA1-209", "Oke_RPN1-80", and "Oke_U2025-86" ####
require(qdap) # will do gsub on vectors

### Scores
## Admiralty
# Oke_ROA1-209
table(CMADMCR13.gcl$scores[,"Oke_ROA1-209",])
CMADMCR13.gcl$scores[,"Oke_ROA1-209",] <- mgsub(pattern = c("A", "G"), replacement = c("V", "F"), text.var = CMADMCR13.gcl$scores[,"Oke_ROA1-209",])
CMADMCR13.gcl$scores[,"Oke_ROA1-209",] <- mgsub(pattern = c("V", "F"), replacement = c("G", "A"), text.var = CMADMCR13.gcl$scores[,"Oke_ROA1-209",])
table(CMADMCR13.gcl$scores[,"Oke_ROA1-209",])

# Oke_RPN1-80
table(CMADMCR13.gcl$scores[,"Oke_RPN1-80",])
CMADMCR13.gcl$scores[,"Oke_RPN1-80",] <- mgsub(pattern = c("A", "G"), replacement = c("V", "F"), text.var = CMADMCR13.gcl$scores[,"Oke_RPN1-80",])
CMADMCR13.gcl$scores[,"Oke_RPN1-80",] <- mgsub(pattern = c("V", "F"), replacement = c("G", "A"), text.var = CMADMCR13.gcl$scores[,"Oke_RPN1-80",])
table(CMADMCR13.gcl$scores[,"Oke_RPN1-80",])

# Oke_U2025-86
table(CMADMCR13.gcl$scores[,"Oke_U2025-86",])
CMADMCR13.gcl$scores[,"Oke_U2025-86",] <- mgsub(pattern = c("A", "G"), replacement = c("V", "F"), text.var = CMADMCR13.gcl$scores[,"Oke_U2025-86",])
CMADMCR13.gcl$scores[,"Oke_U2025-86",] <- mgsub(pattern = c("V", "F"), replacement = c("G", "A"), text.var = CMADMCR13.gcl$scores[,"Oke_U2025-86",])
table(CMADMCR13.gcl$scores[,"Oke_U2025-86",])


## Prospect
# Oke_ROA1-209
table(CMPROSCR13.gcl$scores[,"Oke_ROA1-209",])
CMPROSCR13.gcl$scores[,"Oke_ROA1-209",] <- mgsub(pattern = c("A", "G"), replacement = c("V", "F"), text.var = CMPROSCR13.gcl$scores[,"Oke_ROA1-209",])
CMPROSCR13.gcl$scores[,"Oke_ROA1-209",] <- mgsub(pattern = c("V", "F"), replacement = c("G", "A"), text.var = CMPROSCR13.gcl$scores[,"Oke_ROA1-209",])
table(CMPROSCR13.gcl$scores[,"Oke_ROA1-209",])

# Oke_RPN1-80
table(CMPROSCR13.gcl$scores[,"Oke_RPN1-80",])
CMPROSCR13.gcl$scores[,"Oke_RPN1-80",] <- mgsub(pattern = c("A", "G"), replacement = c("V", "F"), text.var = CMPROSCR13.gcl$scores[,"Oke_RPN1-80",])
CMPROSCR13.gcl$scores[,"Oke_RPN1-80",] <- mgsub(pattern = c("V", "F"), replacement = c("G", "A"), text.var = CMPROSCR13.gcl$scores[,"Oke_RPN1-80",])
table(CMPROSCR13.gcl$scores[,"Oke_RPN1-80",])

# Oke_U2025-86
table(CMPROSCR13.gcl$scores[,"Oke_U2025-86",])
CMPROSCR13.gcl$scores[,"Oke_U2025-86",] <- mgsub(pattern = c("A", "G"), replacement = c("V", "F"), text.var = CMPROSCR13.gcl$scores[,"Oke_U2025-86",])
CMPROSCR13.gcl$scores[,"Oke_U2025-86",] <- mgsub(pattern = c("V", "F"), replacement = c("G", "A"), text.var = CMPROSCR13.gcl$scores[,"Oke_U2025-86",])
table(CMPROSCR13.gcl$scores[,"Oke_U2025-86",])


## Sawmill
# Oke_ROA1-209
table(CMSAWCR13.gcl$scores[,"Oke_ROA1-209",])
CMSAWCR13.gcl$scores[,"Oke_ROA1-209",] <- mgsub(pattern = c("A", "G"), replacement = c("V", "F"), text.var = CMSAWCR13.gcl$scores[,"Oke_ROA1-209",])
CMSAWCR13.gcl$scores[,"Oke_ROA1-209",] <- mgsub(pattern = c("V", "F"), replacement = c("G", "A"), text.var = CMSAWCR13.gcl$scores[,"Oke_ROA1-209",])
table(CMSAWCR13.gcl$scores[,"Oke_ROA1-209",])

# Oke_RPN1-80
table(CMSAWCR13.gcl$scores[,"Oke_RPN1-80",])
CMSAWCR13.gcl$scores[,"Oke_RPN1-80",] <- mgsub(pattern = c("A", "G"), replacement = c("V", "F"), text.var = CMSAWCR13.gcl$scores[,"Oke_RPN1-80",])
CMSAWCR13.gcl$scores[,"Oke_RPN1-80",] <- mgsub(pattern = c("V", "F"), replacement = c("G", "A"), text.var = CMSAWCR13.gcl$scores[,"Oke_RPN1-80",])
table(CMSAWCR13.gcl$scores[,"Oke_RPN1-80",])

# Oke_U2025-86
table(CMSAWCR13.gcl$scores[,"Oke_U2025-86",])
CMSAWCR13.gcl$scores[,"Oke_U2025-86",] <- mgsub(pattern = c("A", "G"), replacement = c("V", "F"), text.var = CMSAWCR13.gcl$scores[,"Oke_U2025-86",])
CMSAWCR13.gcl$scores[,"Oke_U2025-86",] <- mgsub(pattern = c("V", "F"), replacement = c("G", "A"), text.var = CMSAWCR13.gcl$scores[,"Oke_U2025-86",])
table(CMSAWCR13.gcl$scores[,"Oke_U2025-86",])



### Counts
## Admiralty
apply(CMADMCR13.gcl$counts[, SNP.issues[1], c("Allele 1", "Allele 2")], 2, function(allele) {sum(allele, na.rm = TRUE)})
CMADMCR13.gcl$counts[, SNP.issues, ] <- CMADMCR13.gcl$counts[, SNP.issues, c(2, 1, 3)]
apply(CMADMCR13.gcl$counts[, SNP.issues[1], c("Allele 1", "Allele 2")], 2, function(allele) {sum(allele, na.rm = TRUE)})

## Prospect
apply(CMPROSCR13.gcl$counts[, SNP.issues[1], c("Allele 1", "Allele 2")], 2, function(allele) {sum(allele, na.rm = TRUE)})
CMPROSCR13.gcl$counts[, SNP.issues, ] <- CMPROSCR13.gcl$counts[, SNP.issues, c(2, 1, 3)]
apply(CMPROSCR13.gcl$counts[, SNP.issues[1], c("Allele 1", "Allele 2")], 2, function(allele) {sum(allele, na.rm = TRUE)})

## Sawmill
apply(CMSAWCR13.gcl$counts[, SNP.issues[1], c("Allele 1", "Allele 2")], 2, function(allele) {sum(allele, na.rm = TRUE)})
CMSAWCR13.gcl$counts[, SNP.issues, ] <- CMSAWCR13.gcl$counts[, SNP.issues, c(2, 1, 3)]
apply(CMSAWCR13.gcl$counts[, SNP.issues[1], c("Allele 1", "Allele 2")], 2, function(allele) {sum(allele, na.rm = TRUE)})


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get unaltered .gcl object for Fish Creek (CM 40)
CMFISHCR13.gcl <- dget(file = "V:/WORK/Chum/AHRP/Raw genotypes/Phase 2_QCed/CMFISHCR13.txt")
#CMFISHCRT13.gcl <- dget(file = "V:/WORK/Chum/AHRP/Raw genotypes/Phase 2_QCed/CMFISHCRT13.txt") # NOTE: these will not be used since we do not know ancestral origin (no otoliths!)

# Check structure of attributes table
str(CMFISHCR13.gcl$attributes)
names(CMFISHCR13.gcl$attributes)

# Get rid of "added" attributes (can't pool collections if attributes tables are different)
CMFISHCR13.gcl$attributes <- CMFISHCR13.gcl$attributes[, -c(23:30)]


# Confirming that the two loci that had the wrong allele calls (Oke_U1022-114 and Oke_IL8r-272) are OKAY
table(CMADMCR13.gcl$scores[,"Oke_IL8r-272",]); table(CMFISHCR13.gcl$scores[,"Oke_IL8r-272",]) # Looks okay

table(CMADMCR13.gcl$scores[,"Oke_U1022-114",]); table(CMFISHCR13.gcl$scores[,"Oke_U1022-114",]) # Looks okay


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get LocusControl from CM40
LocusControl <- dget(file = "V:/WORK/Chum/AHRP/Objects/Phase 2/OriginalLocusControl.txt")
loci188 <- LocusControl$locusnames



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Need to pull in Finsight otolith-origin data Fri May 08 13:15:52 2015 ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
metadata <- read.csv(file = "Phase 1_StreamSpecimens_20150508.csv", skip = 1)
head(metadata)
str(metadata)

# Add KEY (Individual = Barcode_Well)
metadata$KEY <- apply(cbind(metadata$Sample.Tray.Id, metadata$Sample.Cell), 1, function(row) {paste(row[1:2], collapse = "_")})

# Create "KEY" in attributes table for DNATRAY_WELLNO
CMADMCR13.gcl$attributes$KEY <- apply(cbind(CMADMCR13.gcl$attributes$DNA_TRAY_CODE, CMADMCR13.gcl$attributes$DNA_TRAY_WELL_CODE), 1, function(row) {paste(row[1:2], collapse = "_")})
CMFISHCR13.gcl$attributes$KEY <- apply(cbind(CMFISHCR13.gcl$attributes$DNA_TRAY_CODE, CMFISHCR13.gcl$attributes$DNA_TRAY_WELL_CODE), 1, function(row) {paste(row[1:2], collapse = "_")})
CMPROSCR13.gcl$attributes$KEY <- apply(cbind(CMPROSCR13.gcl$attributes$DNA_TRAY_CODE, CMPROSCR13.gcl$attributes$DNA_TRAY_WELL_CODE), 1, function(row) {paste(row[1:2], collapse = "_")})
CMSAWCR13.gcl$attributes$KEY <- apply(cbind(CMSAWCR13.gcl$attributes$DNA_TRAY_CODE, CMSAWCR13.gcl$attributes$DNA_TRAY_WELL_CODE), 1, function(row) {paste(row[1:2], collapse = "_")})

# Add otolith present
CMADMCR13.gcl$attributes$OTOLITH_MARK_PRESENT <- as.character(metadata$Mark.Present[match(CMADMCR13.gcl$attributes$KEY, metadata$KEY)])
CMFISHCR13.gcl$attributes$OTOLITH_MARK_PRESENT <- as.character(metadata$Mark.Present[match(CMFISHCR13.gcl$attributes$KEY, metadata$KEY)])
CMPROSCR13.gcl$attributes$OTOLITH_MARK_PRESENT <- as.character(metadata$Mark.Present[match(CMPROSCR13.gcl$attributes$KEY, metadata$KEY)])
CMSAWCR13.gcl$attributes$OTOLITH_MARK_PRESENT <- as.character(metadata$Mark.Present[match(CMSAWCR13.gcl$attributes$KEY, metadata$KEY)])

# Add otolith ID
CMADMCR13.gcl$attributes$OTOLITH_MARK_ID <- as.character(metadata$Mark.Id[match(CMADMCR13.gcl$attributes$KEY, metadata$KEY)])
CMFISHCR13.gcl$attributes$OTOLITH_MARK_ID <- as.character(metadata$Mark.Id[match(CMFISHCR13.gcl$attributes$KEY, metadata$KEY)])
CMPROSCR13.gcl$attributes$OTOLITH_MARK_ID <- as.character(metadata$Mark.Id[match(CMPROSCR13.gcl$attributes$KEY, metadata$KEY)])
CMSAWCR13.gcl$attributes$OTOLITH_MARK_ID <- as.character(metadata$Mark.Id[match(CMSAWCR13.gcl$attributes$KEY, metadata$KEY)])

# Add scale age
CMADMCR13.gcl$attributes$SCALE_AGE <- as.character(metadata$Age[match(CMADMCR13.gcl$attributes$KEY, metadata$KEY)])
CMFISHCR13.gcl$attributes$SCALE_AGE <- as.character(metadata$Age[match(CMFISHCR13.gcl$attributes$KEY, metadata$KEY)])
CMPROSCR13.gcl$attributes$SCALE_AGE <- as.character(metadata$Age[match(CMPROSCR13.gcl$attributes$KEY, metadata$KEY)])
CMSAWCR13.gcl$attributes$SCALE_AGE <- as.character(metadata$Age[match(CMSAWCR13.gcl$attributes$KEY, metadata$KEY)])

# Summary stats of metadata
StreamSILLYs <- sapply(objects(pattern = "\\.gcl"), function(silly) {unlist(strsplit(x = silly, split = ".gcl"))[1]})
names(StreamSILLYs) <- StreamSILLYs

sapply(StreamSILLYs, function(silly) {table(get(paste(silly, ".gcl", sep = ""))$attributes$OTOLITH_MARK_PRESENT)})
sapply(StreamSILLYs, function(silly) {table(get(paste(silly, ".gcl", sep = ""))$attributes$OTOLITH_MARK_ID)})
sapply(StreamSILLYs, function(silly) {table(get(paste(silly, ".gcl", sep = ""))$attributes$SCALE_AGE)})
sapply(StreamSILLYs, function(silly) {table(get(paste(silly, ".gcl", sep = ""))$attributes$OTOLITH_MARK_PRESENT, get(paste(silly, ".gcl", sep = ""))$attributes$SCALE_AGE)})



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create new .gcl objects w/ only Natural-origin fish ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get FK_FISH_IDs
Natural_IDs <- sapply(StreamSILLYs, function(silly) {AttributesToIDs.GCL(silly = silly, attribute = "OTOLITH_MARK_PRESENT", matching = "N")})

# Create SILLYs
sapply(StreamSILLYs, function(silly) {PoolCollections.GCL(collections = silly, loci = loci188, IDs = Natural_IDs[silly], newname = paste(silly, "_natural", sep = ""))})

NATsamples <- unlist(strsplit(objects(pattern = "_natural.gcl"), split = ".gcl"))

# Check sample sizes
sapply(NATsamples, function(silly) {get(paste(silly, ".gcl", sep = ""))$n})
sapply(StreamSILLYs, function(silly) {table(get(paste(silly, ".gcl", sep = ""))$attributes$OTOLITH_MARK_PRESENT)}, USE.NAMES = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Create a separate SILLY for DIPAC hatchery fish
DIPAC.BYs <- grep(pattern = "DIPAC", x = unique(unlist(sapply(StreamSILLYs, function(silly) {unique(get(paste(silly, ".gcl", sep = ""))$attributes$OTOLITH_MARK_ID)}, simplify = TRUE), use.names = FALSE)), value = TRUE)
Hatchery_IDs <- sapply(StreamSILLYs, function(silly) {AttributesToIDs.GCL(silly = silly, attribute = "OTOLITH_MARK_ID", matching = DIPAC.BYs)})

PoolCollections.GCL(collections = StreamSILLYs, loci = loci188, IDs = Hatchery_IDs, newname = "DIPAC_strays")
str(DIPAC_strays.gcl)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## All 5 SILLYs of interest (natural + DIPAC)
Phase1samples <- c(NATsamples, "DIPAC_strays")
names(Phase1samples) <- Phase1samples

## Get sample size by locus
OriginalPhase1SampleSizebyLocus <- SampSizeByLocus.GCL(Phase1samples, loci188)
min(OriginalPhase1SampleSizebyLocus) # 37
apply(OriginalPhase1SampleSizebyLocus, 1, min) / apply(OriginalPhase1SampleSizebyLocus, 1, max) # 0.65 for Sawmill due to rotten tissue

## Get number of individuals per silly before removing missing loci individuals
OriginalPhase1ColSize <- sapply(paste(Phase1samples, ".gcl", sep = ''), function(x) get(x)$n)

## Remove individuals with >20% missing data
Phase1MissLoci <- RemoveIndMissLoci.GCL(sillyvec = Phase1samples, loci = loci188, proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSizePhase1PostMissLoci <- sapply(paste(Phase1samples, ".gcl", sep = ''), function(x) get(x)$n)

Phase1SampleSizes <- matrix(data = NA, nrow = 5, ncol = 4, dimnames = list(Phase1samples, c("Genotyped", "Missing", "Duplicate", "Final")))
Phase1SampleSizes[, 1] <- OriginalPhase1ColSize
Phase1SampleSizes[, 2] <- OriginalPhase1ColSize - ColSizePhase1PostMissLoci

## Check within collections for duplicate individuals.
Phase1DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = Phase1samples, loci = loci188, quantile = NULL, minproportion = 0.95)

## Remove duplicate individuals
Phase1RemovedDups <- RemoveDups.GCL(Phase1DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSizePhase1PostDuplicate <- sapply(paste(Phase1samples, ".gcl", sep = ''), function(x) get(x)$n)

Phase1SampleSizes[, 3] <- ColSizePhase1PostMissLoci - ColSizePhase1PostDuplicate
Phase1SampleSizes[, 4] <- ColSizePhase1PostDuplicate

require(xlsx)
write.xlsx(Phase1SampleSizes,file="Data/Phase1SampleSizes.xlsx")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### MAF ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get allele counts/frequencies
x <- FreqPop.GCL(sillyvec = Phase1samples, loci = loci188)
str(x)
p <- x[, , "Allele 1"] / (x[, , "Allele 1"] + x[, , "Allele 2"])

# Create Freq Plots
dir.create(path = "FreqPlots")
pdf(file = paste(getwd(), "/FreqPlots/SEAKChum5Pops185SNPs_FreqPlots.pdf", sep = ""), width = 11, height = 8.5, family = "Times", pointsize = 20)   
for(loci in sort(loci188)) {
  plot(p[, loci], ylim = c(0, 1), ylab = "Frequency", xlab = "Population", axes = FALSE, pch = 16, col = 1:5, cex = 2, main = loci)
  lines(supsmu(x = 1:5, y = p[, loci]))
  axis(side = 2)
  abline(h = 0)
  text(x = 1:5, y = rep(0, 5), labels = c("ADM", "FISH", "PROS", "SAW", "DIPAC"), pos = 1, offset = 1, srt = 45, xpd = TRUE)
  text(x = 1:5, y = rep(0, 5), labels = paste("n = ", Phase1SampleSizes[, "Final"], sep = ""), pos = 1, offset = 2.5, xpd = TRUE)
}; rm(loci)
dev.off()






LocusControl$ploidy["Oke_ROA1-209"]
x[, "Oke_ROA1-209", ]

#### LOCUS CONTROL ISSUES ####
#### I checked the Fluidigm plots and VIC/FAM were flipped for CM38 and CM40
#### I tried to re-read in from ReadLOKI.gcl, but it bombed on Sawmill Creek and also had funny sample sizes
#### This will require further investigation before proceeding

freqs <- array(data = NA, dim = c(length(Phase1samples), length(loci188), 3), dimnames = list(Phase1samples, loci188, c("mean", "5%", "95%")))
for(silly in Phase1samples){
  for(loci in loci188){
    mod <- prop.test(x = x[silly, loci, "Allele 1"], n = x[silly, loci, "Allele 1"] + x[silly, loci, "Allele 2"], conf.level = 0.95, correct = FALSE)
    freqs[silly, loci, ] <- c(mod$estimate, mod$conf.int)
  }
}



pdf(file = paste(getwd(), "/FreqPlots/SEAKChum5Pops185SNPs_FreqPlots95CI.pdf", sep = ""), width = 11, height = 8.5, family = "Times", pointsize = 20)   
for(loci in sort(loci188)) {
  plot(freqs[, loci, "mean"], ylim = c(0, 1), ylab = "Frequency", xlab = "Population", axes = FALSE, pch = 16, col = 1:5, cex = 1.5, main = loci)
  arrows(x0 = 1:5, y0 = freqs[, loci, "5%"], x1 = 1:5, y1 = freqs[, loci, "95%"], angle = 90, code = 3, lwd = 2)
  points(x = 1:5, y = freqs[, loci, "mean"], pch = 16, col = 1:5, cex = 1.5)
  lines(supsmu(x = 1:5, y = p[, loci]))
  axis(side = 2)
  abline(h = 0)
  text(x = 1:5, y = rep(0, 5), labels = c("ADM", "FISH", "PROS", "SAW", "DIPAC"), pos = 1, offset = 1, srt = 45, xpd = TRUE)
  text(x = 1:5, y = rep(0, 5), labels = paste("n = ", Phase1SampleSizes[, "Final"], sep = ""), pos = 1, offset = 2.5, xpd = TRUE)
}; rm(loci)
dev.off()


dir.create("FSTAT")
Phase1.5Pops.188LociSummaryStats <- HoFisFstTable.GCL(sillyvec = Phase1samples, loci = loci188, dir = "FSTAT")

dir.create("Output")
write.table(Phase1.5Pops.188LociSummaryStats, file = "Output/Phase1.5Pops.188LociSummaryStats.xls", sep = "\t", col.names = c("Ho", "Fis", "Fst"), row.names = TRUE)

hist(Phase1.5Pops.188LociSummaryStats$Fst, col = 8)
Phase1.5Pops.188LociSummaryStats["Oke_gdh1-191", "Fst"]
max(Phase1.5Pops.188LociSummaryStats$Fst, na.rm = TRUE)





### Just a sec, need to make sure to call Allele 1 or 2 the right one

# MAF
maf <- pmin(p, 1-p)

# Sample sizes
NATSampleSizes[,4]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### HWE ####
## Calculate HWE for each "Pop"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gcl2Genepop.GCL(sillyvec = Phase1samples, loci = loci188[!LocusControl$ploidy == 1], path = "Data/Phase1samples_185nuclearloci.gen")

# NOTE: DO NOT USE ADEGENET FOR HWE,  USE GENEPOP
# Calculate HWE only with Natural-origin fish and DIPAC strays
# From Tech Doc 2, a marker fails if 1) overall p<0.01, or 2) any pop p<0.05)

source("V:/DATA/R_GEN/tempGCL Source Scripts/ReadGenepopHWEKS.GCL.R")
HWE <- ReadGenepopHWEKS.GCL(file = "Data/Phase1samples_185nuclearloci.txt.P")
str(HWE)
HWE$SummaryPValues
colnames(HWE$SummaryPValues) <- c(Phase1samples, "Overall Pops")
HWE$SummaryPValues <- cbind(HWE$SummaryPValues, "n Pops < 0.05" = apply(HWE$SummaryPValues[, Phase1samples], 1, function(row) {sum(row < 0.05, na.rm = TRUE)}))

# view histogram of HWE p-values
apply(HWE$SummaryPValues[-dim(HWE$SummaryPValues)[1], Phase1samples], 2, function(silly) {hist(silly, col = 8, breaks = seq(from = 0, to = 1, by = 0.05))})


# How many loci with Fischer across pops p < 0.01?
table(HWE$SummaryPValues[-dim(HWE$SummaryPValues)[1], "Overall Pops"] < 0.01)
HWE_across_pops <- rownames(HWE$SummaryPValues)[HWE$SummaryPValues[-dim(HWE$SummaryPValues)[1], "Overall Pops"] < 0.01 & !is.na(HWE$SummaryPValues[-dim(HWE$SummaryPValues)[1], "Overall Pops"])]
HWE_across_pops_mat <- HWE$SummaryPValues[HWE_across_pops, ]


### These were out of HWE in WASSIP but not here...
HWE$SummaryPValues["Oke_U2038-32", ]
maf[, "Oke_U2038-32"] # very low MAF

HWE$SummaryPValues["Oke_U2040-77", ]
maf[, "Oke_U2040-77"] # hrmm



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### LD ####
## Calculate HWE for each "Pop"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


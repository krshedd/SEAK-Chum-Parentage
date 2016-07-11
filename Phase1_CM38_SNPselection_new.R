# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2014 CM 38, Project 37 Chum SNP selection for the AHRP
## 380 fish (95 from each of the 4 pedigree streams), all 188 available chum markers
## Kyle Shedd Tuesday Wed May 27 13:44:19 2015
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


ls()
rm(list=ls(all=TRUE))
search()
getwd()
setwd("V:/WORK/Chum/AHRP")

## save.image("Phase1_CM38_SNPselection_new.RData")
## load("Phase1_CM38_SNPselection_new.RData")

## save.image("V:/WORK/Chum/AHRP/Phase1_CM38_SNPselection.RData")
## load("V:/WORK/Chum/AHRP/Phase1_CM38_SNPselection.RData")

#This sources all of the new GCL functions to this workspace
source("V:/DATA/R_GEN/GCL Source Scripts/Functions.GCL.R")
source("V:/WORK/Kyle/R Source Scripts/Functions.GCL_KS.R")


## Pull all data for each silly code and create .gcl objects for each
AHRPsamples <- c("CMADMCR13","CMFISHCR13","CMPROSCR13","CMSAWCR13")

# Original ====================================================================
## Get genotypes out of LOKI:
ReadLOKI.GCL(sillyvec=AHRPsamples, markersuite="Chum2014_188 SNPs")  # ReadLOKI.GCL will not work in RStudio, need to run in R
objects(pattern="*.gcl")#check objects

# Check sample sizes by SILLY and tissue type
for (collection in AHRPsamples){
  print(get(paste(collection,".gcl",sep=''))$n)
} # Why only 94 for Admirilty and Sawmill? Hans says it is due to random sampling of rows that had NTCs

for (collection in AHRPsamples){
  print(table(get(paste(collection,".gcl",sep=''))$attributes$PK_TISSUE_TYPE))
} # 1 Unkown for Admiralty and 1 Axillary for Sawmill


## Create a locus list from Locus Control object
loci188 <- LocusControl$locusnames

## Save original LocusControl
dput(x=LocusControl,file="Objects/OriginalLocusControl.txt")

## Save unaltered .gcl's as back-up:
for (collection in AHRPsamples){
  dput(x=get(paste(collection,".gcl",sep='')),file=paste("Raw genotypes/",collection,".txt",sep=''))
}

## Get LocusControl and unaltered .gcl's to start from scratch... =============
LocusControl <- dget(file = "Objects/Phase 1/OriginalLocusControl.txt")
loci188 <- LocusControl$locusnames

for(silly in AHRPsamples) {
  assign(x = paste(silly, ".gcl", sep = ""), 
         value = dget(file = paste("Raw genotypes/Phase 1_QCed/", silly, ".txt", sep="")))
}

objects(pattern = "\\.gcl")

################################################################################################################
# Brute force exercises
# determine the number of scored loci per individual
table(CMADMCR13.gcl$counts[,,1]+CMADMCR13.gcl$counts[,,2])
sum(is.na(CMADMCR13.gcl$counts[,,1]+CMADMCR13.gcl$counts[,,2]))

scores_fail_indv=NULL
for(i in 1:94){
  scores_fail_indv[i]=sum(is.na(CMADMCR13.gcl$counts[i,,1]+CMADMCR13.gcl$counts[i,,2]))
};rm(i)
scores_indv=188-scores_fail_indv
scores_indv_percent=round(scores_indv/188*100,2)
hist(scores_indv_percent,breaks=100)
sum(scores_indv_percent<80)

# determine the number of individuals scored per loci
str(CMADMCR13.gcl$counts)
table(CMADMCR13.gcl$counts[,,1]+CMADMCR13.gcl$counts[,,2])
sum(is.na(CMADMCR13.gcl$counts[,,1]+CMADMCR13.gcl$counts[,,2]))

scores_fail_loci=NULL
for(i in 1:188){
  scores_fail_loci[i]=sum(is.na(CMADMCR13.gcl$counts[,i,1]+CMADMCR13.gcl$counts[,i,2]))
}
scores_loci=94-scores_fail_loci
scores_loci_percent=round(scores_loci/94*100,2)
hist(scores_loci_percent,breaks=20)
sum(scores_loci_percent<80)





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PHASE 1: ####
# unranked, veto of markers not in HWE, in LD, or <80% fish genotyped (relative to best marker)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

phase1=array(dim=c(188,6,6),dimnames=list(loci188,c(AHRPsamples,"Overall","Keep"),c("HWE","LD","Success Rate (%)","MAF","sdMAF","Score")))

# Need to pull in Finsight otolith-origin data Fri May 08 13:15:52 2015 =======
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


# Determine Hatchery/Natural numbers
sapply(AHRPsamples, function(silly) {table(get(paste(silly,".gcl",sep=""))$attributes$OTOLITH_MARK_PRESENT)}, USE.NAMES = TRUE)

## Create new .gcl objects w/ only Natural-origin fish ========================
Natural_IDs <- sapply(AHRPsamples, function(silly) {AttributesToIDs.GCL(silly = silly, attribute = "OTOLITH_MARK_PRESENT", matching = "N")})

sapply(AHRPsamples, function(silly) {PoolCollections.GCL(collections = silly, loci = loci188, IDs = Natural_IDs[silly], newname = paste(silly, "_natural", sep = ""))})

NATsamples <- unlist(strsplit(objects(pattern = "_natural.gcl"), split = ".gcl"))

# Check sample sizes
sapply(NATsamples, function(silly) {get(paste(silly, ".gcl", sep = ""))$n})
sapply(AHRPsamples, function(silly) {table(get(paste(silly,".gcl",sep=""))$attributes$OTOLITH_MARK_PRESENT)}, USE.NAMES = TRUE)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data QC/Massage ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### Natural-only ====
## Get sample size by locus
OriginalNATSampleSizebyLocus <- SampSizeByLocus.GCL(NATsamples, loci188)
min(OriginalNATSampleSizebyLocus) ### 29
apply(OriginalNATSampleSizebyLocus,1,min)/apply(OriginalNATSampleSizebyLocus,1,max) ### 0.65 for Sawmill

## Get number of individuals per silly before removing missing loci individuals
OriginalNATColSize <- sapply(paste(NATsamples,".gcl",sep=''), function(x) get(x)$n)

## Remove individuals with >20% missing data
NATMissLoci <- RemoveIndMissLoci.GCL(sillyvec=NATsamples,loci=loci188,proportion=0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSizeNATPostMissLoci <- sapply(paste(NATsamples,".gcl",sep=''), function(x) get(x)$n)

NATSampleSizes <- matrix(data=NA,nrow=4,ncol=4,dimnames=list(NATsamples,c("Genotyped","Missing","Duplicate","Final")))
NATSampleSizes[,1] <- OriginalNATColSize
NATSampleSizes[,2] <- OriginalNATColSize-ColSizeNATPostMissLoci

## Check within collections for duplicate individuals.
NATDuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec=NATsamples,loci=loci188,quantile=NULL,minproportion=0.95)

## Remove duplicate individuals
NATRemovedDups <- RemoveDups.GCL(NATDuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSizeNATPostDuplicate <- sapply(paste(NATsamples,".gcl",sep=''), function(x) get(x)$n)

NATSampleSizes[,3] <- ColSizeNATPostMissLoci-ColSizeNATPostDuplicate
NATSampleSizes[,4] <- ColSizeNATPostDuplicate

require(xlsx)
write.xlsx(NATSampleSizes,file="Data/NATSampleSizes_new.xlsx")

#### Hatchery & Natural ####
## Get sample size by locus
OriginalAHRPSampleSizebyLocus <- SampSizeByLocus.GCL(AHRPsamples, loci188)
min(OriginalAHRPSampleSizebyLocus) ### 56
apply(OriginalAHRPSampleSizebyLocus,1,min)/apply(OriginalAHRPSampleSizebyLocus,1,max) ### 0.60 for Sawmill

## Get number of individuals per silly before removing missing loci individuals
OriginalAHRPColSize <- sapply(paste(AHRPsamples,".gcl",sep=''), function(x) get(x)$n)

## Remove individuals with >20% missing data
AHRPMissLoci <- RemoveIndMissLoci.GCL(sillyvec=AHRPsamples,loci=loci188,proportion=0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSizeAHRPPostMissLoci <- sapply(paste(AHRPsamples,".gcl",sep=''), function(x) get(x)$n)

AHRPSampleSizes <- matrix(data=NA,nrow=4,ncol=4,dimnames=list(AHRPsamples,c("Genotyped","Missing","Duplicate","Final")))
AHRPSampleSizes[,1] <- OriginalAHRPColSize
AHRPSampleSizes[,2] <- OriginalAHRPColSize-ColSizeAHRPPostMissLoci

## Check within collections for duplicate individuals.
AHRPDuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec=AHRPsamples,loci=loci188,quantile=NULL,minproportion=0.95)

## Remove duplicate individuals
AHRPRemovedDups <- RemoveDups.GCL(AHRPDuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSizeAHRPPostDuplicate <- sapply(paste(AHRPsamples,".gcl",sep=''), function(x) get(x)$n)

AHRPSampleSizes[,3] <- ColSizeAHRPPostMissLoci-ColSizeAHRPPostDuplicate
AHRPSampleSizes[,4] <- ColSizeAHRPPostDuplicate

require(xlsx)
write.xlsx(AHRPSampleSizes,file="Data/AHRPSampleSizes_new.xlsx")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## HWE ####
## Calculate HWE only with Natural-origin fish
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gcl2Genepop.GCL(sillyvec=NATsamples,loci=loci188,path="Data/NATsamples_new.gen")

# NOTE: DO NOT USE ADEGENET FOR HWE,  USE GENEPOP
# Calculate HWE only with Natural-origin fish
# From Tech Doc 2, a marker fails if 1) overall p<0.01, or 2) any pop p<0.05)

source("V:/DATA/R_GEN/tempGCL Source Scripts/ReadGenepopHWEKS.GCL.R")
HWE <- ReadGenepopHWEKS.GCL(file = "Data/NATsamples_new.txt.P")
str(HWE)
HWE$SummaryPValues
colnames(HWE$SummaryPValues) <- c(AHRPsamples, "Overall Pops")
HWE$SummaryPValues <- cbind(HWE$SummaryPValues, "n Pops < 0.05" = apply(HWE$SummaryPValues[, 1:4], 1, function(row) {sum(row < 0.05)}))

# How many loci with Fischer across pops p < 0.01?
table(HWE$SummaryPValues[-dim(HWE$SummaryPValues)[1], "Overall Pops"] < 0.01)
HWE_across_pops <- rownames(HWE$SummaryPValues)[HWE$SummaryPValues[-dim(HWE$SummaryPValues)[1], "Overall Pops"] < 0.01 & !is.na(HWE$SummaryPValues[-dim(HWE$SummaryPValues)[1], "Overall Pops"])]
HWE_across_pops_mat <- HWE$SummaryPValues[HWE_across_pops, ]


### These were out of HWE in WASSIP but not here...
HWE$SummaryPValues["Oke_U2038-32", ]
maf[, "Oke_U2038-32"] # very low MAF

HWE$SummaryPValues["Oke_U2040-77", ]
maf[, "Oke_U2040-77"] # hrmm


### This marker was out of HWE here, but not in WASSIP (not in ANY pops)
HWE$SummaryPValues["Oke_U1012-60", ]
maf[, "Oke_U1012-60"] # very high MAF
# Look into why Oke_U1012-60 was out of HWE?
HWE$DataByPop[HWE$DataByPop$Locus == "Oke_U1012-60", ] # very negative Fis for Admiralty might indicate an excess of hets scoring error!!!



# How many loci with n pops p < 0.05
table(HWE$SummaryPValues[-dim(HWE$SummaryPValues)[1], "n Pops < 0.05"])
HWE_within_pops <- rownames(HWE$SummaryPValues)[HWE$SummaryPValues[-dim(HWE$SummaryPValues)[1], "n Pops < 0.05"] >= 2 & !is.na(HWE$SummaryPValues[-dim(HWE$SummaryPValues)[1], "n Pops < 0.05"])]
HWE$SummaryPValues[HWE_within_pops, ]

# Output summary for conditional formatting in Excel
write.csv(x = HWE$SummaryPValues, file = "Data/NATsamples_new_HWE.csv")

# Get rid of markers out of HWE in > 2 pops or overall < 0.01
failed_loci_HWE <- unique(c(HWE_across_pops, HWE_within_pops))
HWE$SummaryPValues[failed_loci_HWE, ]

# Look into why Oke_U1012-60 was out of HWE?
HWE$DataByPop[HWE$DataByPop$Locus == "Oke_U1012-60", ]
# Check others for out of whack Fis
table(HWE$DataByPop[, "WC Fis"] < -0.5)
hist(HWE$DataByPop[, "WC Fis"], breaks = seq(from = -1, to = 1, by = 0.1), col = 8, xlab = "Fis", main = "Histogram of Fis")

# I looked into how this marker was scored and given:
#  1) the .bml genotyping plot,
#  2) no other collections were out of HWE, and
#  3) that this marker was not out of HWE for ANY of the 29 WASSIP populations
# I am not going to remove this marker for parentage

failed_loci_HWE <- failed_loci_HWE[-3]
HWE$SummaryPValues[failed_loci_HWE, ]


# IMPORTANT Decision Point ====================================================
# Earlier on Mon May 11 16:58:21 2015 Chris Habicht and I (Kyle) decided that
# the way to deal with HWE and LD for this project given the limited sample 
# size (didn't account for the fact that hatchery fish would be removed) and
# some poor quality genotyping was to reject and loci that are out of HWE or in
# LD for this project OR WASSIP, thus keeping the most conservative markers

HWE$SummaryPValues[failed_loci_HWE, ]

failed_loci_HWE_WASSIP <- c("Oke_RS27-94", "Oke_U1103-150", "Oke_U2038-32", "Oke_U2040-77", "Oke_XBP1-82")
failed_loci_HWE_WASSIP %in% loci188 # all properly spelled

failed_loci_HWE_final <- unique(c(failed_loci_HWE, failed_loci_HWE_WASSIP))
# 6 loci fail for HWE

# Look at HWE for all AHRP samples, not just natural-origin fish. See if it changes anything. KS on Wed Jul 08 14:55:43 2015
date()

source("V:/DATA/R_GEN/tempGCL Source Scripts/ReadGenepopHWEKS.GCL.R")
HWE <- ReadGenepopHWEKS.GCL(file = "Data/OLD/AHRPsamples.txt.P")
str(HWE)
HWE$SummaryPValues
colnames(HWE$SummaryPValues) <- c(AHRPsamples, "Overall Pops")
HWE$SummaryPValues <- cbind(HWE$SummaryPValues, "n Pops < 0.05" = apply(HWE$SummaryPValues[, 1:4], 1, function(row) {sum(row < 0.05)}))

# How many loci with Fischer across pops p < 0.01?
table(HWE$SummaryPValues[-dim(HWE$SummaryPValues)[1], "Overall Pops"] < 0.01)
HWE_across_pops <- rownames(HWE$SummaryPValues)[HWE$SummaryPValues[-dim(HWE$SummaryPValues)[1], "Overall Pops"] < 0.01 & !is.na(HWE$SummaryPValues[-dim(HWE$SummaryPValues)[1], "Overall Pops"])]
HWE_across_pops_mat <- HWE$SummaryPValues[HWE_across_pops, ]


### These were out of HWE in WASSIP but not here...
HWE$SummaryPValues["Oke_U2038-32", ]
maf[, "Oke_U2038-32"] # very low MAF

HWE$SummaryPValues["Oke_U2040-77", ]
maf[, "Oke_U2040-77"] # hrmm


### This marker was out of HWE here, but not in WASSIP (not in ANY pops)
HWE$SummaryPValues["Oke_U1012-60", ]
maf[, "Oke_U1012-60"] # very high MAF
# Look into why Oke_U1012-60 was out of HWE?
HWE$DataByPop[HWE$DataByPop$Locus == "Oke_U1012-60", ] # very negative Fis for Admiralty might indicate an excess of hets scoring error!!!



# How many loci with n pops p < 0.05
table(HWE$SummaryPValues[-dim(HWE$SummaryPValues)[1], "n Pops < 0.05"])
HWE_within_pops <- rownames(HWE$SummaryPValues)[HWE$SummaryPValues[-dim(HWE$SummaryPValues)[1], "n Pops < 0.05"] >= 2 & !is.na(HWE$SummaryPValues[-dim(HWE$SummaryPValues)[1], "n Pops < 0.05"])]
HWE$SummaryPValues[HWE_within_pops, ]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## LD ####
## Calculate LD only with Natural-origin fish
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# My adaptation to Andy's function, gives a data.frame of p-values for Pop and Overall
source("V:/DATA/R_GEN/tempGCL Source Scripts/ReadGenepopDisKS.GCL.R")
LD <- ReadGenepopDisKS.GCL(file = "Data/NATsamples_new.txt.DIS") 
str(LD)
head(LD); tail(LD)

# Overall p-values
table(LD$Overall < 0.001)
LD[which(LD$Overall < 0.001), ]

# Which markers are p < 0.05 at > half pops
LD_within_pops <- LD[which(apply(LD[, 3:6], 1, function(locuspair) {sum(locuspair < 0.05) > 2})), ]
LD_within_pops

# Sorted table
t(t(sort(table(unlist(LD[which(apply(LD[, 3:6], 1, function(locuspair) {sum(locuspair < 0.05) > 2})), 1:2]))))) # FUUUUUCK

write.csv(x = LD[which(apply(LD[, 3:6], 1, function(locuspair) {sum(locuspair < 0.05) > 2})), ], file = "Data/NATsamples_new_DIS.csv")

## NOTE: NATsamples_new_DIS.xlsx has the association blocks ====
# Compare these pairs of loci in LD to what was cound in WASSIP
# RIR.5J.2012.25 Table 4

# Published name vs. unpublished
LocusControl$locusnames[which(LocusControl$Publishedlocusnames == "Oke_U401-143")]
LocusControl$locusnames[which(LocusControl$Publishedlocusnames == "Oke_U401-220")]

# This was signfiicant in WASSIP
LD[which(LD$Locus1 == "Oke_U1021-102" & LD$Locus2 == "Oke_U1022-139"), ]
# This wasn't but "should" have been
LD[which(LD$Locus1 == "Oke_U1021-102" & LD$Locus2 == "Oke_U1022-114"), ]

# This was significant in this, but NOT in WASSIP
LD[which(LD$Locus1 == "Oke_U2048-91" & LD$Locus2 == "Oke_NUPR1-70"), ]


failed_loci_LD <- sort(unique(c(LD_within_pops$Locus1, LD_within_pops$Locus2)))
failed_loci_LD_WASSIP <- sort(c("Oke_ATP5L-105", "Oke_ATP5L-248", "Oke_CCT3-143", "Oke_CCT3-220", "Oke_DBLOH-79", "Oke_gdh1-191", "Oke_gdh1-234", "Oke_gdh1-62", "Oke_IL8r-272", "Oke_IL8r2-406", "Oke_NHERF-123", "Oke_NHERF-54", "Oke_PDIA3-475", "Oke_PDIA3-82", "Oke_pgap-111", "Oke_pgap-92",  "Oke_psmd9-188", "Oke_psmd9-57", "Oke_TCTA-202", "Oke_TCTA-99", "Oke_u0602-244", "Oke_U1001-79", "Oke_U1002-165", "Oke_U1002-262", "Oke_U1010-251", "Oke_U1012-241", "Oke_U1012-60", "Oke_U1022-114", "Oke_U1022-139", "Oke_U507-286", "Oke_U507-87", "Oke_UBA3-245", "Oke_U1021-102"))

failed_loci_LD_tentative <- unique(c(failed_loci_LD, failed_loci_LD_WASSIP))


#### Determine which to keep based on which marker has higher MAF score =======
## Natural-origin only!
# Check sample sizes
NATSampleSizes[,4]

# Get allele counts/frequencies
x <- FreqPop.GCL(sillyvec=NATsamples,loci=loci188)
str(x)
p <- x[, , "Allele 1"] / (x[, , "Allele 1"] + x[, , "Allele 2"])

# MAF
maf <- pmin(p, 1-p)

# Sample sizes
NATSampleSizes[,4]

# Package for weighted mean and sd, weighting by sample size (number of fish)
require(SDMTools)
wt.maf.NAT <- apply(maf, 2, function(loci) {wt.mean(x = loci, wt = NATSampleSizes[,4])})
wt.sd.NAT <- apply(maf, 2, function(loci) {wt.sd(x = loci, wt = NATSampleSizes[,4])})

score <- (2 * wt.maf.NAT) / (1 + wt.sd.NAT)

hist(score)

# Get scores for LD loci
failed_loci_LD_tentative_scores <- score[names(score) %in% failed_loci_LD_tentative]

t(t(failed_loci_LD_tentative_scores))

# Add scores to df
LD_within_pops$Locus1.Score <- failed_loci_LD_tentative_scores[match(LD_within_pops$Locus1, names(failed_loci_LD_tentative_scores))]
LD_within_pops$Locus2.Score <- failed_loci_LD_tentative_scores[match(LD_within_pops$Locus2, names(failed_loci_LD_tentative_scores))]

# Which loci has highest score?
str(LD_within_pops)
LD_within_pops$LocusWinner <- apply(LD_within_pops, 1, function(comp) {comp[which.max(comp[8:9])]})


# Create a list of "association blocks"
association <- list(data.frame("loci" = failed_loci_LD_tentative[c(23, 5, 22)]),
                    data.frame("loci" = failed_loci_LD_tentative[c(29, 30, 35)]),
                    data.frame("loci" = failed_loci_LD_tentative[c(6:8)]),
                    data.frame("loci" = failed_loci_LD_tentative[c(24:25)]),
                    data.frame("loci" = failed_loci_LD_tentative[c(27:28)]),
                    data.frame("loci" = failed_loci_LD_tentative[c(32:33)]),
                    data.frame("loci" = failed_loci_LD_tentative[c(34, 26)]),
                    data.frame("loci" = failed_loci_LD_tentative[c(3:4)]),
                    data.frame("loci" = failed_loci_LD_tentative[c(14:15)]),
                    data.frame("loci" = failed_loci_LD_tentative[c(20:21)]),
                    data.frame("loci" = failed_loci_LD_tentative[c(1:2)]),
                    data.frame("loci" = failed_loci_LD_tentative[c(11:12)]),
                    data.frame("loci" = failed_loci_LD_tentative[c(18:19)]),
                    data.frame("loci" = failed_loci_LD_tentative[c(16:17)]),
                    data.frame("loci" = failed_loci_LD_tentative[c(9:10)]),
                    data.frame("loci" = failed_loci_LD_tentative[c(31, 13)])
)
# Assign scores
association_scores <- lapply(association, function(block) {failed_loci_LD_tentative_scores[match(block$loci, names(failed_loci_LD_tentative_scores))]})
# Keep loci with best score
keep_loci_LD_score <- unlist(lapply(association_scores, function(block) {names(which.max(block))}))



#### Try with just straight weighted MAF from the SDMTools package ============

# Weighted mean p weighted by number of allele calls per locus (locus specific sample size) for comparison
w.maf <- apply(maf * (x[, , "Allele 1"] + x[, , "Allele 2"]), 2, sum) / apply(x, 2, function(loci) sum(loci, na.rm = TRUE))
# SDMTools vs. by hand (this is the difference between locus specific sample size and number of fish sample size)
table(round(w.maf - wt.maf.NAT, 3)) # very minimal differences

# Get wt.maf for LD loci
failed_loci_LD_tentative_wt.maf <- wt.maf.NAT[names(wt.maf.NAT) %in% failed_loci_LD_tentative]

# Assign wt.maf
association_wt.maf <- lapply(association, function(block) {failed_loci_LD_tentative_wt.maf[match(block$loci, names(failed_loci_LD_tentative_wt.maf))]})
# Keep loci with best score
keep_loci_LD_wt.maf <- unlist(lapply(association_wt.maf, function(block) {names(which.max(block))}))

## Compare score vs. wt.maf
cbind(keep_loci_LD_score, keep_loci_LD_wt.maf)
outersect(keep_loci_LD_score, keep_loci_LD_wt.maf)
# No diffrences between keeping loci by "score" or straight weighted MAF (wt.maf)


# Which loci did WASSIP keep?
keep_loci_LD_WASSIP <- list(failed_loci_LD_WASSIP[c(5)],
                            failed_loci_LD_WASSIP[c(28, 30)],
                            failed_loci_LD_WASSIP[c(6, 8)],
                            failed_loci_LD_WASSIP[c(24)],
                            failed_loci_LD_WASSIP[c(26)],
                            failed_loci_LD_WASSIP[c(31)],
                            failed_loci_LD_WASSIP[c(25)],
                            failed_loci_LD_WASSIP[c(3)],
                            failed_loci_LD_WASSIP[c(13)],
                            failed_loci_LD_WASSIP[c(19)],
                            failed_loci_LD_WASSIP[c(1)],
                            failed_loci_LD_WASSIP[c(12)],
                            failed_loci_LD_WASSIP[c(18)],
                            failed_loci_LD_WASSIP[c(15:16)],
                            failed_loci_LD_WASSIP[c(10)])
)

keep_loci_LD_score %in% unlist(keep_loci_LD_WASSIP)


### Drop loci due to LD ====
failed_loci_LD_final <- outersect(failed_loci_LD_tentative, keep_loci_LD_score)


## All loci that fail the "gating" measures ====
length(failed_loci_HWE_final); length(failed_loci_LD_final)
failed_loci_ALL <- unique(c(failed_loci_HWE_final, failed_loci_LD_final))


## Loci that passed the "gating" measures ====
loci163 <- loci188[!loci188 %in% failed_loci_ALL]
dput(x = loci163, file = "V:/WORK/Chum/AHRP/Objects/loci163.txt")


# MAF of loci kept
hist(wt.maf.NAT[loci188])
hist(wt.maf.NAT[loci163], add = TRUE, col = 8)

# Ploidy of loci163
table(LocusControl$ploidy[loci163])



## How many of the 188 should we order? Which ones? ===========================
# None of the ones that failed at the gating measures
# Only the ones with "reasonably high"  MAF, but what is that?

# How many of the 163 have an MAF > 0.1?
sum(wt.maf.NAT[loci163] > 0.1)
# > 0.05?
sum(wt.maf.NAT[loci163] > 0.05)

# Think about exclusions rates by MAF...
MAF <- seq(from = 0.01, to = 0.5, by = 0.01)

## Diploid
# MAF and power for single parents (can only exclude when find alternate homozygotes for parent and offspring, hets are worthless)
# 2 * p^2 * q^2
plot(2 * (MAF^2) * (1 - MAF)^2, type = "l", lwd = 5, bty = "n", axes = FALSE, cex.lab = 2, xlab = "MAF", ylab = "Probability")
axis(side = 1, labels = seq(from = 0, to = 0.5, by = 0.1), at = seq(from = 0, to = 50, by = 10), cex.axis = 2)
axis(side = 2, cex.axis = 1.5)
abline(h = 0.005, lwd = 5, col = 2) # avg genotyping error rate for SNPs

## Haploid (only good for Females)
# MAF and power for single parents (can only exclude when find alternate homozygotes for parent and offspring, hets are worthless)
# There are no hets for haploid markers, so here we have 2pq (divided by 2 to penalize for only females) (2 * p * q) / 2
plot((MAF) * (1 - MAF), type = "l", lwd = 5, bty = "n", axes = FALSE, cex.lab = 2, xlab = "MAF", ylab = "Probability")
axis(side = 1, labels = seq(from = 0, to = 0.5, by = 0.1), at = seq(from = 0, to = 50, by = 10), cex.axis = 2)
axis(side = 2, cex.axis = 1.5)
abline(h = 0.005, lwd = 5, col = 2) # avg genotyping error rate for SNPs


## Both together
par(mar = c(5.1, 5.1, 2.1, 2.1))
plot(2 * (MAF) * (1 - MAF) / 2, type = "l", lwd = 5, lty = 2, bty = "n", axes = FALSE, cex.lab = 2, xlab = "MAF", ylab = "Naive Exclusion Probability", ylim = c(0, 0.25)) # / 2 is a penalty, because this only helps with Females, not Males
points(2 * (MAF^2) * (1 - MAF)^2, type = "l", lwd = 5)
axis(side = 1, labels = seq(from = 0, to = 0.5, by = 0.1), at = seq(from = 0, to = 50, by = 10), cex.axis = 2, pos = 0, lwd = 3)
axis(side = 2, cex.axis = 1.5, at = seq(from = 0, to = 0.25, by = 0.05), pos = 0, lwd = 3)
abline(h = 0.005, lwd = 5, col = 2) # avg genotyping error rate for SNPs, conservatively very high given that we care about homo-homo conflicts or miscoring mitochondrial markers
legend("topleft", bty = "n", legend = c("Haploid", "Diploid"), lty = c(2, 1), lwd = 5, cex = 2.5)




# If you have a marker panel where ALL SNPs are of a certain MAF, how likely would it be to find no exclusions for a false parent?
nSNPs <- 96
plot((1 - (2 * (MAF^2) * (1 - MAF)^2))^nSNPs, type = "l", lwd = 5, bty = "n", axes = FALSE, cex.lab = 2, xlab = "MAF", ylab = "Probability")
axis(side = 1, labels = seq(from = 0, to = 0.5, by = 0.1), at = seq(from = 0, to = 50, by = 10), cex.axis = 2)
axis(side = 2, cex.axis = 1.5)

# Power of of the 129 vs. 163 marker sets
prod(1 - ((wt.maf.NAT[loci163]^2) * (1 - wt.maf.NAT[loci163])^2))
prod(1 - ((wt.maf.NAT[loci129]^2) * (1 - wt.maf.NAT[loci129])^2))


# Likelihood of no mismatches by chance by number of markers (high-grading)
plot(
  sapply(seq_along(loci163), function(nloci){
    loci <- sort(wt.maf.NAT[loci163], decreasing = TRUE)[1:nloci]
    prod(1 - ((loci^2) * (1 - loci)^2))
    }),
type = "l", lwd = 5, xlab = "nSNPs", ylab = "Probability")



# Met with Chris and decided to keep markers with a MAF of > 0.05 for now =====
loci129 <- loci163[which(wt.maf.NAT[loci163] > 0.05)]
dput(x = loci129, file = "V:/WORK/Chum/AHRP/Objects/loci129.txt")


# Double check that all are Kosher
sort(HWE$SummaryPValues[loci129, "Overall Pops"])
sort(wt.maf.NAT[loci129])

# Check MAF in ALL Fish Creek adults (2013)
maf_Fish <- dget(file = "V:/WORK/Chum/AHRP/Objects/Phase 2/FishCreekAllParents_MAF.txt")
cbind(wt.maf.NAT[loci129], maf["CMFISHCR13_natural", loci129], maf_Fish[loci129])
hist(maf_Fish[loci129]); min(maf_Fish[loci129])

# Okay, send it out!
cbind(common = loci129, published = LocusControl$Publishedlocusnames[which(LocusControl$locusnames %in% loci129)])

#### What about haploid mitochondrial markers !!!!=============================
# Chat with Jim Jasper about mitochondrial markers made me re-consider their utility in parentage analysis
# Check and see if I threw out any mitochondrial markers

table(LocusControl$ploidy[loci188]) # started with 3
table(LocusControl$ploidy[loci163]) # didn't lose any to HWE or LD
table(LocusControl$ploidy[loci129]) # none
table(LocusControl$ploidy[outersect(loci129, loci163)]) # all 3 were eliminated due to low MAF

loci188[which(LocusControl$ploidy == 1)]

wt.maf.NAT[loci188[which(LocusControl$ploidy == 1)]] # appears fixed
maf[, loci188[which(LocusControl$ploidy == 1)]] # appears fixed
x[, loci188[which(LocusControl$ploidy == 1)], ] # yup, it is fixed

# Check in all the Fish Creek data to confirm that these three are fixed
maf_Fish <- dget(file = "V:/WORK/Chum/AHRP/Objects/Phase 2/FishCreekAllParents_MAF.txt")
maf_Fish[loci188[which(LocusControl$ploidy == 1)]] # appears fixed

## Check in WASSIP
# Ran this code in V:\WORK\WASSIP\Chum\Baseline\WASSIPchumBaseline.R
# Data is in load("V:/Work/WASSIP/chum/Baseline/WASSIPchumBaseline.RData")

# Admiralty
WASSIP432Freq[grep(pattern = "ADM", x = WASSIP436Collections, value = TRUE), c("Oke_Cr30", "Oke_Cr386", "Oke_ND3???69"), ] # fixed
# Fish Creek
WASSIP432Freq[grep(pattern = "CMFCJUN06S", x = WASSIP436Collections, value = TRUE), c("Oke_Cr30", "Oke_Cr386", "Oke_ND3???69"), ] # fixed
# DIPAC
WASSIP432Freq[grep(pattern = "CMDIPAC06", x = WASSIP436Collections, value = TRUE), c("Oke_Cr30", "Oke_Cr386", "Oke_ND3???69"), ] # fixed
# Prospect
WASSIP432Freq[grep(pattern = "PRO", x = WASSIP436Collections, value = TRUE), c("Oke_Cr30", "Oke_Cr386", "Oke_ND3???69"), ]
# Sawmill
WASSIP432Freq[grep(pattern = "SAW", x = WASSIP436Collections, value = TRUE), c("Oke_Cr30", "Oke_Cr386", "Oke_ND3???69"), ]

# Well, nothing to worry about then, MAF = 0, thus, no power

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stopped Here ================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






# How many pops in LD for each pair
LD$fail=apply(LD[,3:6]<0.05,1,function(row){sum(row)})
table(LD$fail)

# Bonferroni correct
table(LD[,"Overall"]<(0.05/dim(LD)[1]))

# Which markers were involved in any LD for each pop
phase1[,1,"LD"]=rownames(phase1)%in%unique(c(as.character(LD[which(LD$CMADMCR13<0.05),1]),as.character(LD[which(LD$CMADMCR13<0.05),2]))) # 0 = not in any LD, 1 = involved in some LD
phase1[,2,"LD"]=rownames(phase1)%in%unique(c(as.character(LD[which(LD$CMFISHCR13<0.05),1]),as.character(LD[which(LD$CMFISHCR13<0.05),2]))) # 0 = not in any LD, 1 = involved in some LD
phase1[,3,"LD"]=rownames(phase1)%in%unique(c(as.character(LD[which(LD$CMPROSCR13<0.05),1]),as.character(LD[which(LD$CMPROSCR13<0.05),2]))) # 0 = not in any LD, 1 = involved in some LD
phase1[,4,"LD"]=rownames(phase1)%in%unique(c(as.character(LD[which(LD$CMSAWCR13<0.05),1]),as.character(LD[which(LD$CMSAWCR13<0.05),2]))) # 0 = not in any LD, 1 = involved in some LD
phase1[,5,"LD"]=rownames(phase1)%in%unique(c(as.character(LD[which(LD$Overall<0.05),1]),as.character(LD[which(LD$Overall<0.05),2]))) # 0 = not in any LD, 1 = involved in some LD
phase1[,,2]

# How many pairs in LD for 3 pops
failat3=LD[which(LD$fail>=3),]
unique(c(as.character(failat3[,1]),as.character(failat3[,2])))

locifailat3=data.frame(Locus1=as.character(failat3[,1]),Locus2=as.character(failat3[,2]))
locifailat3=data.frame(lapply(locifailat3,as.character),stringsAsFactors=FALSE)

y=FreqPop.GCL(sillyvec=NATsamples,loci=loci188)
py=y[,,1]/(y[,,1]+y[,,2])
NATMAF=t(pmin(py,1-py))
NATMAF=cbind(NATMAF,apply(NATMAF,1,function(row){mean(row,na.rm=TRUE)}))

locifailat3$MAF1=NATMAF[match(locifailat3$Locus1,rownames(NATMAF)),5]
locifailat3$MAF2=NATMAF[match(locifailat3$Locus2,rownames(NATMAF)),5]
locifailat3$fail=ifelse(locifailat3$MAF1>locifailat3$MAF2,locifailat3$Locus2,locifailat3$Locus1)

# Markers that fail at 3 marker LD cutoff
unique(locifailat3$fail) # 18

# How many pairs in LD for 2 pops
failat2=LD2[which(LD2$fail>=2),]
unique(c(as.character(failat2[,1]),as.character(failat2[,2])))

locifailat2=data.frame(Locus1=as.character(failat2[,1]),Locus2=as.character(failat2[,2]))
locifailat2=data.frame(lapply(locifailat2,as.character),stringsAsFactors=FALSE)

locifailat2$MAF1=NATMAF[match(locifailat2$Locus1,rownames(NATMAF)),5]
locifailat2$MAF2=NATMAF[match(locifailat2$Locus2,rownames(NATMAF)),5]
locifailat2$fail=ifelse(locifailat2$MAF1>locifailat2$MAF2,locifailat2$Locus2,locifailat2$Locus1)

phase1[,6,"LD"]=rownames(phase1)%in%locifailat2$fail # 0 = not in any LD, 1 = involved in some LD

# Markers that fail at 2 marker LD cutoff
unique(locifailat2$fail) # 99

## Sticking with 3 marker LD cutoff
phase1[,6,"LD"]=rownames(phase1)%in%locifailat3$fail # 0 = not in any LD, 1 = involved in some LD
phase1[,,2]

table(apply(phase1[,6,1:2],1,function(row){sum(row)}))

keeploci=rownames(phase1)[which(apply(phase1[,6,1:2],1,function(row){sum(row)})==0)] # keep 83 if use 2 marker LD cutoff, keep 149 if use 3 marker LD cutoff

# MAF for 83 loci kept with 2 marker LD cutoff
phase1[,6,"LD"]=rownames(phase1)%in%locifailat2$fail # 0 = not in any LD, 1 = involved in some LD
keeploci=rownames(phase1)[which(apply(phase1[,6,1:2],1,function(row){sum(row)})==0)] # keep 83 if use 2 marker LD cutoff, keep 149 if use 3 marker LD cutoff
mean(phase1[match(keeploci,rownames(phase1)),5,4]) # 0.237

# MAF for 149 loci kept with 3 marker LD cutoff
phase1[,6,"LD"]=rownames(phase1)%in%locifailat3$fail # 0 = not in any LD, 1 = involved in some LD
keeploci=rownames(phase1)[which(apply(phase1[,6,1:2],1,function(row){sum(row)})==0)] # keep 83 if use 2 marker LD cutoff, keep 149 if use 3 marker LD cutoff
mean(phase1[match(keeploci,rownames(phase1)),5,4]) # 0.235, but this is all 149, not the "best" 96
mean(sort(phase1[match(keeploci,rownames(phase1)),5,4],decreasing=TRUE)[1:96]) # 0.337 for naive "best" 96 (ranked by MAF, not "score")

## Check HWE with all fish (natural and hatchery)
gcl2Genepop.GCL(sillyvec=AHRPsamples,loci=loci188,path="Data/AHRPsamples.gen")
LD2all=ReadGenepopDis.GCL(file="Data/NATsamples.txt.DIS") 
# How many pops in LD for each pair
LD2all$fail=apply(LD2all[,3:6]<0.05,1,function(row){sum(row)})
table(LDall2$fail)


## Success Rate ###############################################################################
## Calculate with ALL fish (Hatchery and Natural) #############################################

## Get sample size by locus
OriginalAHRPSamplesSampleSizebyLocus <- SampSizeByLocus.GCL(AHRPsamples,loci188)
min(OriginalAHRPSamplesSampleSizebyLocus) ## 0, we had 5 markers fail to load on chip, have been run off chip, but waiting on data to be loaded into LOKI, loaded as of 8/6/14 09:18:42
apply(OriginalAHRPSamplesSampleSizebyLocus,1,min)/apply(OriginalAHRPSamplesSampleSizebyLocus,1,max) ## gives the worst marker for each SILLY; ADM, FISH, and SAW had marker load failures

## Get number of individuals per silly before removing missing loci individuals
OriginalAHRPSamplesColSize <- sapply(paste(AHRPsamples,".gcl",sep=''), function(x) get(x)$n)
OriginalAHRPPercentbyLocus=round(OriginalAHRPSamplesSampleSizebyLocus/OriginalAHRPSamplesColSize*100,2)
# Transpose to show loci as row
t(OriginalAHRPPercentbyLocus)

# Number of loci that were <80% success by silly
apply(OriginalAHRPPercentbyLocus,1,function(x) sum(x<80)) # only Sawmill, due to rotten fish / poor tissue quality (see Tech Doc 3)
# Number of loci that were <80% success by silly, excluding 5 markers that failed on chip
apply(OriginalAHRPPercentbyLocus,1,function(x) sum(x<80 & x>0)) # This was to exclude failed markers

# Average marker (loci) success across silly
OriginalAHRPPercentbyLocusOverall=round(matrix(apply(OriginalAHRPSamplesSampleSizebyLocus,2,sum),dimnames=list(names(OriginalAHRPSamplesSampleSizebyLocus),"%"))/sum(OriginalAHRPSamplesColSize)*100,2)
# Number of markers that failed, <80% success across all fish sampled
sum(OriginalAHRPPercentbyLocusOverall<80) # HURRAY, all markers were good enough
summary(OriginalAHRPPercentbyLocusOverall) # mean success = 92.15
# Which 5 markers failed?
rownames(OriginalAHRPPercentbyLocusOverall)[which(OriginalAHRPPercentbyLocusOverall<80)]
# ADM - U1012.241; FISH - rab5a.117, NHERF.54, & RS9.379; & SAW - txnrdl.74
x=apply(t(OriginalAHRPPercentbyLocus), 1, function(row) all(row !=0 ))
names(x)[which(t(t(x))=="FALSE")];rm(x)
# These are the 5 loci that failed to load on chip and have been re-run off chip
# ADM - U1012.241; FISH - rab5a.117, NHERF.54, & RS9.379; & SAW - txnrdl.74
# These data were loaded into LOKI as of 8/6/14 9:24:56

# Add Success Rate to big Phase 1 array
phase1[,1:5,3]=cbind(t(OriginalAHRPPercentbyLocus),OriginalAHRPPercentbyLocusOverall)
phase1[,6,3]=ifelse(phase1[,5,3]<80,1,0)
phase1[,,1:3]
phase1[,6,]
keeploci=rownames(phase1)[which(apply(phase1[,6,1:3],1,function(row){sum(row)})==0)]





#######################################################################################################################################################################################################################
PHASE 1: ranked, by MAF, sd MAF, & score= (2*MAF)/(1+SD MAF) ###################################################################################################################################################
#######################################################################################################################################################################################################################

str(LocusControl)
mitoloci=names(which(LocusControl$ploidy==1))
str(x)
x[,mitoloci,] # mitochondrial loci are all fixed for all 4 pops

## MAF ########################################################################################

## All samples (hatchery + natural)
x=FreqPop.GCL(sillyvec=AHRPsamples,loci=loci188)
str(x)
p=x[,,1]/(x[,,1]+x[,,2])

# Plot allele frequencies
maf <- pmin(p, 1-p)

for(loci in loci188){
  plot(maf[, loci], pch = 16, ylim = c(0, 0.5), xlab = "", ylab = "MAF", bty = "n", axes = FALSE, cex = 3, main = loci, cex.main = 2)
  axis(side = 2)
  axis(side = 1, labels = NA, at = 1:4)
  text(x = 1:4, y = rep(-0.07, 4), labels = AHRPsamples, srt = 45, pos = 1, xpd = TRUE)
  abline(h = 0.2, col = 2, lwd = 6)
}

## Natural-origin only!
# Check sample sizes
sapply(NATsamples, function(silly) {get(paste(silly, ".gcl", sep = ""))$n})

# Get allele frequencies
x=FreqPop.GCL(sillyvec=NATsamples,loci=loci188)
str(x)
p=x[,,1]/(x[,,1]+x[,,2])

# Plot allele frequencies
maf <- pmin(p, 1-p)

# All
for(loci in loci188){
  plot(maf[, loci], pch = 16, ylim = c(0, 0.5), xlab = "", ylab = "MAF", bty = "n", axes = FALSE, cex = 3, main = loci, cex.main = 2)
  axis(side = 2)
  axis(side = 1, labels = NA, at = 1:4)
  text(x = 1:4, y = rep(-0.07, 4), labels = AHRPsamples, srt = 45, pos = 1, xpd = TRUE)
  abline(h = 0.2, col = 2, lwd = 6)
}


# Specific
loci <- "Oke_RFC2"

plot(p[, loci], pch = 16, ylim = c(0, 1), xlab = "", ylab = "MAF", bty = "n", axes = FALSE, cex = 3, main = loci, cex.main = 2)
axis(side = 2)
axis(side = 1, labels = NA, at = 1:4)
text(x = 1:4, y = rep(-0.07, 4), labels = AHRPsamples, srt = 45, pos = 1, xpd = TRUE)
abline(h = 0.2, col = 2, lwd = 6)

loci188[which(LocusControl$Publishedlocusnames == loci)]






phase1[,1:4,4]=t(pmin(p,1-p))
phase1[,5,4]=apply(phase1[,1:4,4],1,function(x) mean(x,na.rm=TRUE))

phase1[match(keeploci,rownames(phase1)),5,4]


## SD MAF ##################################################################################

phase1[,5,5]=apply(phase1[,1:4,4],1,function(x) sd(x,na.rm=TRUE))


## SCORE ###################################################################################

phase1[,5,6]=(2*phase1[,5,4])/(1+phase1[,5,5])

#View a summary of the results
phase1[,5,]
phase1[,5,4:6]=round(phase1[,5,4:6],3)

# Histogram of MAF, all 188 markers
hist(phase1[,5,4])
mean(phase1[,5,4]) # 0.243

# MAF of top 96 markers (inluding those that fail unranked portion!!!)
t(t(t(t(sort(phase1[,5,4],decreasing=TRUE)))[1:96,1]))
# mean MAF of top 96 markers (inluding those that fail unranked portion!!!)
mean(t(t(t(t(sort(phase1[,5,4],decreasing=TRUE)))[1:96,1]))) #0.376
plot(sort(phase1[,5,4]),type="h")

## WHICH MARKERS PASS UNRANKED

rownames(phase1)%in%keeploci
phase1passed=phase1[match(keeploci,rownames(phase1)),,]
str(phase1passed)

# MAF for 149 loci that weren't vetoed by unranked measures (HWE, LD, and Success Rate)
mean(phase1passed[,5,4]) # 0.235

# Best XX loci by score for 149 loci that passed veto by unranked measures (HWE, LD, and Success Rate)
topscoreloci24=names(sort(phase1passed[,5,6],decreasing=TRUE)[1:24])
topscoreloci48=names(sort(phase1passed[,5,6],decreasing=TRUE)[1:48])
topscoreloci72=names(sort(phase1passed[,5,6],decreasing=TRUE)[1:72])
topscoreloci96=names(sort(phase1passed[,5,6],decreasing=TRUE)[1:96])
topscoreloci120=names(sort(phase1passed[,5,6],decreasing=TRUE)[1:120])
topscoreloci144=names(sort(phase1passed[,5,6],decreasing=TRUE)[1:144])

# MAF for best XX loci
mean(phase1passed[match(topscoreloci24,rownames(phase1passed)),5,4]) # 0.455
hist(phase1passed[match(topscoreloci24,rownames(phase1passed)),5,4],breaks=20)
mean(phase1passed[match(topscoreloci48,rownames(phase1passed)),5,4]) # 0.426
hist(phase1passed[match(topscoreloci48,rownames(phase1passed)),5,4],breaks=20)
mean(phase1passed[match(topscoreloci96,rownames(phase1passed)),5,4]) # 0.337
hist(phase1passed[match(topscoreloci96,rownames(phase1passed)),5,4],breaks=20)
mean(phase1passed[match(topscoreloci120,rownames(phase1passed)),5,4]) # 0.288
hist(phase1passed[match(topscoreloci120,rownames(phase1passed)),5,4],breaks=20)
mean(phase1passed[match(topscoreloci144,rownames(phase1passed)),5,4]) # 0.243
hist(phase1passed[match(topscoreloci144,rownames(phase1passed)),5,4],breaks=20)

getwd()
dput(x=topscoreloci24, file="V:/WORK/Chum/AHRP/Objects/Phase 1/topscoreloci24.txt")
dput(x=topscoreloci48, file="V:/WORK/Chum/AHRP/Objects/Phase 1/topscoreloci48.txt")
dput(x=topscoreloci72, file="V:/WORK/Chum/AHRP/Objects/Phase 1/topscoreloci72.txt")
dput(x=topscoreloci96, file="V:/WORK/Chum/AHRP/Objects/Phase 1/topscoreloci96.txt")
dput(x=topscoreloci120, file="V:/WORK/Chum/AHRP/Objects/Phase 1/topscoreloci120.txt")
dput(x=topscoreloci144, file="V:/WORK/Chum/AHRP/Objects/Phase 1/topscoreloci144.txt")


###############################################################################################################################################################
STOPPED HERE
###############################################################################################################################################################
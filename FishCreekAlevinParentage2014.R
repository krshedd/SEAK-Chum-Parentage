#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2014 Fish Creek Parentage Analysis ####
# Kyle Shedd Mon Jul 11 11:14:11 2016
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
date()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Introduction ####
# The goal of this script is to perform parentage analysis on chum salmon
# adults (2013) and alevin (2014). Parentage analysis will be performed with
# different sets of SNPs to select a final marker set.
# 1) Gating measures (HWE and LD)
# 3) Ranking measures (MAF)
# 2) Parentage analysis (subsets of SNPs)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Specific Objectives for Parentage Analysis ####
# This script will:
# 1) Import adult/alevin data
# 2) Define potential parents and offspring
# 3) Perform a data QC on mixtures
# 4) Prepare FRANz input files
# 5) Summarize FRANz results
# 6) Run fullsnplings for alevin
# 7) Generate plots and tables of results

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Initial Setup ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))

options(java.parameters = "-Xmx100g")

setwd("V:/Analysis/5_Coastwide/Multispecies/Alaska Hatchery Research Program/SEAK Chum")
source("H:/R Source Scripts/Functions.GCL_KS.R")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
username <- "krshedd"
password <- "********"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get Data from LOKI ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get collection SILLYs
SEAKChumSillys <- c("CMFISHCR14a", "CMFISHCR13", "CMFISHCRT13", "CMADMCR13", "CMSAWCR13", "CMPROSCR13")

## Create Locus Control
CreateLocusControl.GCL(markersuite = "ChumParentage2015_284SNPs", username = username, password = password)

## Save original LocusControl
loci284 <- LocusControl$locusnames
mito.loci <- which(LocusControl$ploidy == 1)

dir.create("Objects")
dput(x = LocusControl, file = "Objects/OriginalLocusControl.txt")
dput(x = loci284, file = "Objects/loci284.txt")
dput(x = mito.loci, file = "Objects/mito.loci.txt")

#~~~~~~~~~~~~~~~~~~
## Pull all data for each silly code and create .gcl objects for each
# sillyvec = SEAKChumSillys; username = username; password = password
LOKI2R.GCL(sillyvec = SEAKChumSillys, username = username, password = password)  # Had to bust open the function and run line by line with `options(java.parameters = "-Xmx100g")` as opposed to `10g`, otherwise hit GC overhead and run out of heap space
rm(username, password)
objects(pattern = "\\.gcl")

## Save unaltered .gcl's as back-up:
dir.create("Raw genotypes")
dir.create("Raw genotypes/OriginalCollections")
invisible(sapply(SEAKChumSillys, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections/" , silly, ".txt", sep = ''))} )); beep(8)

## Original sample sizes by SILLY
collection.size.original <- sapply(SEAKChumSillys, function(silly) get(paste(silly, ".gcl", sep = ""))$n)



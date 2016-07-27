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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Add Metatdata from OceanAK ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Add bottle data from OceanAK - CMFISHCR14a
oceanak.CMFISHCR14a.dt <- fread(input = "OceanAK/Bulk Tissue Inventory 20160726.txt")
oceanak.CMFISHCR14a.df <- data.frame(oceanak.CMFISHCR14a.dt)
str(oceanak.CMFISHCR14a.df)

bottles <- oceanak.CMFISHCR14a.df$Bottle.Name  # bottle names
number.fish.bottle <- apply(oceanak.CMFISHCR14a.df, 1, function(btl) {length(btl["First.Vial"]:btl["Last.Vial"])} )  # n fish per bottle
bottle.id <- rep(bottles, number.fish.bottle)  # vector of bottle names
length(bottle.id) == CMFISHCR14a.gcl$n  # confirm equal
CMFISHCR14a.gcl$attributes$BOTTLE_ID <- bottle.id  # add bottle ID to each individual


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sillys to filter OceanAK Salmon Biological Fact for
writeClipboard(paste(SEAKChumSillys, collapse = ";"))

# Only have data for CMADMCR13, CMFISHCR13, CMPROSCR13, and CMSAWCR13 (not CMFISHCRT13 since they were tagged alive)
require(data.table)
oceanak.dt <- fread(input = "OceanAK/Salmon Biological Data 20160726.txt")  # freaky fast
str(oceanak.dt)
# Convert to data.frame
oceanak.df <- data.frame(oceanak.dt)
str(oceanak.df)

# Create data keys for both (barcode + position)
oceanak.df$Key <- paste(oceanak.df$DNA.Tray.Code, oceanak.df$DNA.Tray.Well.Code, sep = "_")

CMFISHCR13.gcl$attributes$Key <- paste(CMFISHCR13.gcl$attributes$DNA_TRAY_CODE, CMFISHCR13.gcl$attributes$DNA_TRAY_WELL_CODE, sep = "_")
CMFISHCRT13.gcl$attributes$Key <- paste(CMFISHCRT13.gcl$attributes$DNA_TRAY_CODE, CMFISHCRT13.gcl$attributes$DNA_TRAY_WELL_CODE, sep = "_")
CMADMCR13.gcl$attributes$Key <- paste(CMADMCR13.gcl$attributes$DNA_TRAY_CODE, CMADMCR13.gcl$attributes$DNA_TRAY_WELL_CODE, sep = "_")
CMPROSCR13.gcl$attributes$Key <- paste(CMPROSCR13.gcl$attributes$DNA_TRAY_CODE, CMPROSCR13.gcl$attributes$DNA_TRAY_WELL_CODE, sep = "_")
CMSAWCR13.gcl$attributes$Key <- paste(CMSAWCR13.gcl$attributes$DNA_TRAY_CODE, CMSAWCR13.gcl$attributes$DNA_TRAY_WELL_CODE, sep = "_")

# Note that floy tag individuals from CMFISHCRT13 were not pulled into CMFISHCR13 by LOKI2R
table(CMFISHCR13.gcl$attributes$PK_TISSUE_TYPE)

# Which floy tag individuals were recapured and had their otoliths read? (CMFISHCR13 individuals with tissue = floy tag)
oceanak.CMFISHCR13.dt <- fread(input = "OceanAK/GEN_SAMPLED_FISH_TISSUE 20160726.txt")
oceanak.CMFISHCR13.df <- data.frame(oceanak.CMFISHCR13.dt)
str(oceanak.CMFISHCR13.df)

oceanak.CMFISHCR13.df$Key <- paste(oceanak.CMFISHCR13.df$DNA_TRAY_CODE, oceanak.CMFISHCR13.df$DNA_TRAY_WELL_CODE, sep = "_")

CMFISHCR13T.recaptures <- CMFISHCRT13.gcl$attributes$SillySource[match(oceanak.CMFISHCR13.df$CAPTURE_LOCATION, CMFISHCRT13.gcl$attributes$CAPTURE_LOCATION)]


match(CMFISHCR13.gcl$attributes$Key, oceanak.df$Key)
match(oceanak.CMFISHCR13.df$Key, oceanak.df$Key)
match(CMADMCR13.gcl$attributes$Key, oceanak.df$Key)
match(CMPROSCR13.gcl$attributes$Key, oceanak.df$Key)
match(CMSAWCR13.gcl$attributes$Key, oceanak.df$Key)















# Determin DWP tray and position ID to get olotith metadata ===================
str(CMFISHCR13.gcl)

ADULTsamples = c("CMFISHCR13")

ChumSamples=data.frame()
for (collection in ADULTsamples){
  ChumSamples=rbind(ChumSamples,get(paste(collection,".gcl",sep=""))$attributes[c(4,20:22)])
}
write.xlsx(ChumSamples,"Data/CM40_SamplesforParentage.xlsx")

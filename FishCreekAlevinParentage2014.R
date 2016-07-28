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
dput(x = SEAKChumSillys, file = "Objects/SEAKChumSillys.txt")

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
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/5_Coastwide/Multispecies/Alaska Hatchery Research Program/SEAK Chum")
# This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

## Get objects
LocusControl <- dget(file = "Objects/OriginalLocusControl.txt")

SEAKobjects <- list.files(path = "Objects", recursive = FALSE)
SEAKobjects <- SEAKobjects[!SEAKobjects %in% c("OriginalLocusControl.txt")]
SEAKobjects

invisible(sapply(SEAKobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)


## Get un-altered mixtures
invisible(sapply(SEAKChumSillys, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/OriginalCollections/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Add Metatdata from OceanAK ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(data.table)

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


#~~~~~~~~~~~~~~~~~~
# Note that floy tag individuals from CMFISHCRT13 were not pulled into CMFISHCR13 by LOKI2R
table(CMFISHCR13.gcl$attributes$PK_TISSUE_TYPE)

# Which floy tag individuals were recapured and had their otoliths read? (CMFISHCR13 individuals with tissue = floy tag)
oceanak.CMFISHCR13.dt <- fread(input = "OceanAK/GEN_SAMPLED_FISH_TISSUE 20160726.txt")
oceanak.CMFISHCR13.df <- data.frame(oceanak.CMFISHCR13.dt)
str(oceanak.CMFISHCR13.df); dim(oceanak.CMFISHCR13.df)

oceanak.CMFISHCR13.df$Key <- paste(oceanak.CMFISHCR13.df$DNA_TRAY_CODE, oceanak.CMFISHCR13.df$DNA_TRAY_WELL_CODE, sep = "_")

# CMFISHCRT13.recaptures.Key <- CMFISHCRT13.gcl$attributes$Key[match(oceanak.CMFISHCR13.df$CAPTURE_LOCATION, CMFISHCRT13.gcl$attributes$CAPTURE_LOCATION)]


# Pool CMFISHCRT13 fish that were recaptured
CMFISHCRT13.recapture.IDs <- list(CMFISHCRT13 = na.omit(as.character(CMFISHCRT13.gcl$attributes$FK_FISH_ID[
  match(oceanak.CMFISHCR13.df$CAPTURE_LOCATION, CMFISHCRT13.gcl$attributes$CAPTURE_LOCATION)]
  )))
PoolCollections.GCL(collections = "CMFISHCRT13", loci = loci284, CMFISHCRT13.recapture.IDs, newname = "CMFISHCRT13_recapture")
CMFISHCRT13_recapture.gcl$attributes$Key <- paste(CMFISHCRT13_recapture.gcl$attributes$DNA_TRAY_CODE, CMFISHCRT13_recapture.gcl$attributes$DNA_TRAY_WELL_CODE, sep = "_")
str(CMFISHCRT13_recapture.gcl)

# Note that the Key codes for the axillaries do not match the otoliths
table(oceanak.CMFISHCR13.df$Key %in% oceanak.df$Key)
table(CMFISHCRT13_recapture.gcl$attributes$Key %in% oceanak.df$Key)

# Need to replace the dimnames for counts/scores with "true" fish numbers and also replace attributes
dimnames(CMFISHCRT13_recapture.gcl$counts)[[1]] <- oceanak.CMFISHCR13.df$FK_FISH_ID
dimnames(CMFISHCRT13_recapture.gcl$scores)[[1]] <- oceanak.CMFISHCR13.df$FK_FISH_ID
names(CMFISHCRT13_recapture.gcl$attributes)
names(oceanak.CMFISHCR13.df)

str(CMFISHCRT13_recapture.gcl$attributes)

CMFISHCRT13_recapture.gcl$attributes[, c("FK_FISH_ID", "DNA_TRAY_CODE", "DNA_TRAY_WELL_CODE", "DNA_TRAY_WELL_POS", "Key")] <- 
  oceanak.CMFISHCR13.df[, c("FK_FISH_ID", "DNA_TRAY_CODE", "DNA_TRAY_WELL_CODE", "DNA_TRAY_WELL_POS", "Key")]

str(CMFISHCRT13_recapture.gcl)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Matching
SEAKChumSillys.otoltihs <- c("CMFISHCR13", "CMFISHCRT13_recapture", "CMADMCR13", "CMPROSCR13", "CMSAWCR13")

# Create list of matches by silly
SEAKChumSillys.otoltihs.match <- sapply(SEAKChumSillys.otoltihs, function(silly) {
  match(get(paste(silly, ".gcl", sep = ''))$attributes$Key, oceanak.df$Key)
}, simplify = FALSE)
str(SEAKChumSillys.otoltihs.match)

str(oceanak.df)

save.image(file = "FishCreekAlevinParentage2014.RData")

# Add field/otolith metadata using match
require(qdap)
sapply(SEAKChumSillys.otoltihs, function(silly) {
  my.gcl <- get(paste(silly, ".gcl", sep = ''))
  atts <- c("Sex", "Length.Mm", "Otolith.Mark.Present", "Otolith.Mark.ID", "Otolith.Mark.Status.Code")
  my.gcl$attributes <- cbind(my.gcl$attributes, oceanak.df[SEAKChumSillys.otoltihs.match[[silly]], atts])
  my.gcl$attributes$Natural.Hatchery <- mgsub(pattern = c("NO", "YES"), replacement = c("N", "H"), text.var = my.gcl$attributes$Otolith.Mark.Present)
  my.gcl$attributes$Natural.Hatchery[my.gcl$attributes$Natural.Hatchery == ""] = "U"
  assign(x = paste(silly, ".gcl", sep = ''), value = my.gcl, pos = 1)
})
str(CMFISHCR13.gcl)
table(CMFISHCR13.gcl$attributes$Otolith.Mark.Present); table(CMFISHCR13.gcl$attributes$Natural.Hatchery)
unique(CMFISHCR13.gcl$attributes$Otolith.Mark.Present)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save post-match, pre-QC .gcl objects
dir.create("Raw genotypes/PostMetadataPreQC")
invisible(sapply(SEAKChumSillys.otoltihs, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/PostMetadataPreQC/" , silly, ".txt", sep = ''))} )); beep(8)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SEAKChumSillys.all <- sapply(objects(pattern = "\\.gcl"), function(gclobject) {unlist(strsplit(x = gclobject, split = ".gcl"))} )
samp.size <- sapply(SEAKChumSillys.all, function(silly) {get(paste(silly, ".gcl", sep = ''))$n} )

require(xlsx); require(beepr)

SEAKChumSillys.all

SEAKChumSillys.all_SampleSizes <- matrix(data = NA, nrow = length(SEAKChumSillys.all), ncol = 5, 
                                        dimnames = list(SEAKChumSillys.all, c("Genotyped", "Alternate", "Missing", "Duplicate", "Final")))

#### Check loci
## Get sample size by locus
Original_SEAKChumSillys.all_SampleSizebyLocus <- SampSizeByLocus.GCL(SEAKChumSillys.all, loci284)
min(Original_SEAKChumSillys.all_SampleSizebyLocus)  # 0

sort(apply(Original_SEAKChumSillys.all_SampleSizebyLocus,1,min)/apply(Original_SEAKChumSillys.all_SampleSizebyLocus,1,max))  # Several under 0.8
table(apply(Original_SEAKChumSillys.all_SampleSizebyLocus,1,min)/apply(Original_SEAKChumSillys.all_SampleSizebyLocus,1,max) < 0.8)  # all 7 SILLY's with at least one locus fail

# Remove loci that failed for all individuals
loci2remove <- loci284[apply(Original_SEAKChumSillys.all_SampleSizebyLocus, 2, sum) == 0]
loci277 <- loci284[!loci284 %in% loci2remove]

# Percent of individuals genotyped per locus
Original_SEAKChumSillys.all_SampleSizebyLocus <- SampSizeByLocus.GCL(SEAKChumSillys.all, loci277)
SEAKChumSillys.all.percent.per.locus <- apply(Original_SEAKChumSillys.all_SampleSizebyLocus, 2, function(locus) {locus / samp.size} )
require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(SEAKChumSillys.all.percent.per.locus, col.regions = new.colors, xlab = "SILLY", ylab = "Locus", at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 90)), aspect = "fill")  # aspect = "iso" will make squares

# View histogram of number of individuals genotyped per loci
sapply(SEAKChumSillys.all, function(silly) {hist(as.numeric(Original_SEAKChumSillys.all_SampleSizebyLocus[silly, ]), main = silly, col = 8, xlab = "Sample size per locus")} )



#### Check individuals
### Initial
## Get number of individuals per silly before removing missing loci individuals
Original_SEAKChumSillys.all_ColSize <- sapply(paste(SEAKChumSillys.all, ".gcl", sep = ''), function(x) get(x)$n)
SEAKChumSillys.all_SampleSizes[, "Genotyped"] <- Original_SEAKChumSillys.all_ColSize


### Alternate
## Indentify alternate species individuals
ptm <- proc.time()
SEAKChumSillys.all_Alternate <- FindAlternateSpecies.GCL(sillyvec = SEAKChumSillys.all, species = "chum"); beep(8)
proc.time() - ptm

## Remove Alternate species individuals
RemoveAlternateSpecies.GCL(AlternateSpeciesReport = SEAKChumSillys.all_Alternate, AlternateCutOff = 0.5, FailedCutOff = 0.5); beep(2)

## Get number of individuals per silly after removing alternate species individuals
ColSize_SEAKChumSillys.all_PostAlternate <- sapply(paste(SEAKChumSillys.all, ".gcl", sep = ''), function(x) get(x)$n)
SEAKChumSillys.all_SampleSizes[, "Alternate"] <- Original_SEAKChumSillys.all_ColSize-ColSize_SEAKChumSillys.all_PostAlternate


### Missing
## Remove individuals with >20% missing data
SEAKChumSillys.all_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = SEAKChumSillys.all, proportion = 0.8); beep(8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_SEAKChumSillys.all_PostMissLoci <- sapply(paste(SEAKChumSillys.all, ".gcl", sep = ''), function(x) get(x)$n)
SEAKChumSillys.all_SampleSizes[, "Missing"] <- ColSize_SEAKChumSillys.all_PostAlternate-ColSize_SEAKChumSillys.all_PostMissLoci


### Duplicate
## Check within collections for duplicate individuals.
SEAKChumSillys.all_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = SEAKChumSillys.all, loci = loci277, quantile = NULL, minproportion = 0.95); beep(8)
SEAKChumSillys.all_DuplicateCheckReportSummary <- sapply(SEAKChumSillys.all, function(x) SEAKChumSillys.all_DuplicateCheck95MinProportion[[x]]$report)

## Remove duplicate individuals
SEAKChumSillys.all_RemovedDups <- RemoveDups.GCL(SEAKChumSillys.all_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_SEAKChumSillys.all_PostDuplicate <- sapply(paste(SEAKChumSillys.all, ".gcl", sep = ''), function(x) get(x)$n)
SEAKChumSillys.all_SampleSizes[, "Duplicate"] <- ColSize_SEAKChumSillys.all_PostMissLoci-ColSize_SEAKChumSillys.all_PostDuplicate


### Final
SEAKChumSillys.all_SampleSizes[, "Final"] <- ColSize_SEAKChumSillys.all_PostDuplicate
SEAKChumSillys.all_SampleSizes

dir.create("Output")
write.xlsx(SEAKChumSillys.all_SampleSizes, file = "Output/SEAKChumSillys.all_SampleSizes.xlsx")

## Save post-QC .gcl's as back-up:
dir.create(path = "Raw genotypes/PostQCCollections")
invisible(sapply(SEAKChumSillys.all, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/PostQCCollections/" , silly, ".txt", sep = ''))} )); beep(8)


# Percent of individuals genotyped per locus
Original_SEAKChumSillys.all_postQC_SampleSizebyLocus <- SampSizeByLocus.GCL(SEAKChumSillys.all[-7], loci277)
samp.size <- sapply(SEAKChumSillys.all[-7], function(silly) {get(paste(silly, ".gcl", sep = ''))$n} )
SEAKChumSillys.all.percent.per.locus <- apply(Original_SEAKChumSillys.all_postQC_SampleSizebyLocus, 2, function(locus) {locus / samp.size} )
require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(SEAKChumSillys.all.percent.per.locus, col.regions = new.colors, xlab = "SILLY", ylab = "Locus", at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 90)), aspect = "fill")  # aspect = "iso" will make squares


# Hatchery vs. Natural per silly
sapply(SEAKChumSillys.all, function(silly) {table(get(paste(silly, ".gcl", sep = ''))$attributes$Natural.Hatchery)} )


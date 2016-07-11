## 2014 CM 40 Chum SNP selection for the AHRP
## Parentage analsysis of Fish Creek 2013 adults and 2014 alevin
## Kyle Shedd Tuesday December 9 2014 09:08:25
#==============================================================================

ls()
rm(list=ls(all=TRUE))
search()
getwd()
setwd("V:/WORK/Chum/AHRP")
#setwd("F:/2015 AHRP April Seattle")


# Load/Save ===================================================================
## save.image("V:/WORK/Chum/AHRP/Phase2_CM40_SNPselection_AlevinParentage_QCed_FRANz.RData")
## load("V:/WORK/Chum/AHRP/Phase2_CM40_SNPselection_AlevinParentage_QCed_FRANz.RData")

## save.image("F:/2015 AHRP April Seattle/Phase2_CM40_SNPselection_AlevinParentage_QCed_FRANz.RData")
## load("F:/2015 AHRP April Seattle/Phase2_CM40_SNPselection_AlevinParentage_QCed_FRANz.RData")


# Source scripts ==============================================================
#This sources all of the GCL functions to this workspace
#source("V:/DATA/R_GEN/GCL Source Scripts/Functions.GCL.R")                      # GCL functions
#source("V:/WORK/Kyle/R Source Scripts/Functions.GCL_KS.R")                      # Kyle's functions
#source("V:/WORK/Pink/AHRG/Parentage simulations/R Functions/HW_Functions.R")    # Kyle's functions for parentage analysis

source("F:/2015 AHRP April Seattle/R Functions/HW_Functions.R")
source("F:/2015 AHRP April Seattle/GCL Source Scripts/Functions.GCL.R")


## Pull all data for each silly code and create .gcl objects for each
PEDsamples = c("CMFISHCR13","CMFISHCRT13","CMFISHCR14a")

## Get genotypes out of LOKI:
ReadLOKI.GCL(sillyvec=PEDsamples, markersuite="Chum2014_188 SNPs")  # ReadLOKI.GCL will not work in RStudio, need to run in R
objects(pattern="*.gcl")#check objects
objects(pattern="\\.gcl")

# Check sample sizes by SILLY and tissue type
sampsize=sapply(paste(PEDsamples,".gcl",sep=''), function(SILLY) {get(SILLY)$n})

sapply(paste(PEDsamples,".gcl",sep=''), function(SILLY) {table(get(SILLY)$attributes$PK_TISSUE_TYPE)})
#### NOTE THAT CMFISHCR13 ONLY CONTAINS FISH WITH GENOTYPES, SO THAT IS WHY THERE ARE 0 FLOYTAG FISH

## Create a locus list from Locus Control object
loci188 <- LocusControl$locusnames

## Save original LocusControl
dput(x=LocusControl,file="Objects/Phase 2/OriginalLocusControl_QCed.txt")

## Save unaltered .gcl's as back-up:
for (collection in PEDsamples){
  dput(x=get(paste(collection,".gcl",sep='')),file=paste("Raw genotypes/Phase 2_QCed/",collection,".txt",sep=''))
}


## Add Metatdata ==============================================================
# Add bottle data from OceanAK - CMFISHCR14a ==================================
str(CMFISHCR14a.gcl$attributes)
bottle=readClipboard() # read list of bottle B1-B69
number=as.numeric(readClipboard()) # read list of lab count per bottle
dput(x=rep(bottle,number),file=paste("Objects/Phase 2/BottleID.txt"))

bottle_id=dget(file="Objects/Phase 2/BottleID.txt")
CMFISHCR14a.gcl$attributes$BOTTLE_ID=bottle_id # add bottle ID to each individual
str(CMFISHCR14a.gcl$attributes)


# Determin DWP tray and position ID to get olotith metadata ===================
str(CMFISHCR13.gcl)

ADULTsamples = c("CMFISHCR13")

ChumSamples=data.frame()
for (collection in ADULTsamples){
  ChumSamples=rbind(ChumSamples,get(paste(collection,".gcl",sep=""))$attributes[c(4,20:22)])
}
write.xlsx(ChumSamples,"Data/CM40_SamplesforParentage.xlsx")
#### NOTE THAT THIS LIST DOES NOT INCLUDE FLOY TAG INDIVIDUALS, AS THEY DO NOT HAVE A GENOTYPE FOR CMFISHCR13, AND THUS WERE NOT PULLED WITH ReadLOKI.GCL
## Instead I pulled a tissue report from OceanAK and sent that to Eric. I sent him a list of fish with barcode # and position # (i.e. 1300001382_1, 1300001382_2, etc.)


# Add field data such as SEX and ORIGIN designation from OTOLITH data =========
# Data source is "V:/WORK/Chum/AHRP/Data/CM40_oto_data_20141209_SQLpull_KS.xlsx" coppied ID; OTOLITH_MARK_PRESENT; and OTOLITH_MARK_ID columns
# Note that ID is derived from DNA_TRAY_CODE + _ + DNA_TRAY_WELL_CODE

#### DATA KEY =================================================================
## Create FISH_BARCODE attribute which is = DNATRAY_POSITION and is the data key
ParentsSILLY=c("CMFISHCR13","CMFISHCRT13")

SILLY=paste(ParentsSILLY[1],".gcl",sep='')
CMFISHCR13.gcl$attributes$FISH_BARCODE=unlist(apply(cbind(get(SILLY)$attributes$DNA_TRAY_CODE,get(SILLY)$attributes$DNA_TRAY_WELL_CODE),1,function(row){paste(row[1:2],collapse="_")}))
str(CMFISHCR13.gcl$attributes)

SILLY=paste(ParentsSILLY[2],".gcl",sep='')
CMFISHCRT13.gcl$attributes$FISH_BARCODE=unlist(apply(cbind(get(SILLY)$attributes$DNA_TRAY_CODE,get(SILLY)$attributes$DNA_TRAY_WELL_CODE),1,function(row){paste(row[1:2],collapse="_")}))
str(CMFISHCRT13.gcl$attributes)

#### OTOLITH ==================================================================
## Load data
# Get OTOLITH_ORIGIN data from field from Eric L. SQL querry on 12-09-2014
require(xlsx)
CMFISHCR13_oto=read.xlsx(file="Data/CM40_oto_data_20141209_SQLpull_KS.xlsx",sheetName="Kyle")
str(CMFISHCR13_oto);head(CMFISHCR13_oto)

## Pair otolith data
# Add OTOLITH_ORIGIN data to CMFISHCR13.gcl attributes by matching FISH_BARCODE
match(CMFISHCR13.gcl$attributes$FISH_BARCODE,as.character(CMFISHCR13_oto$FISH_BARCODE))
CMFISHCR13.gcl$attributes$OTOLITH_ORIGIN=as.character(CMFISHCR13_oto$OTOLITH_MARK_PRESENT[match(CMFISHCR13.gcl$attributes$FISH_BARCODE,as.character(CMFISHCR13_oto$FISH_BARCODE))])
str(CMFISHCR13.gcl$attributes)
CMFISHCR13.gcl$attributes$OTOLITH_ORIGIN

CMFISHCR13.gcl$attributes$OTOLITH_ORIGIN[is.na(CMFISHCR13.gcl$attributes$OTOLITH_ORIGIN)]="U"
CMFISHCR13.gcl$attributes$OTOLITH_ORIGIN[which(CMFISHCR13.gcl$attributes$OTOLITH_ORIGIN=="YES")]="H"
CMFISHCR13.gcl$attributes$OTOLITH_ORIGIN[which(CMFISHCR13.gcl$attributes$OTOLITH_ORIGIN=="NO")]="N"
table(CMFISHCR13.gcl$attributes$OTOLITH_ORIGIN)

# Add OTOLITH_ORIGIN data to CMFISHCRT13.gcl attributes by matching FISH_BARCODE
# Matched duplicate FISH_BARCODES by floy tag number to assign origin data to CMFISHCRT13.gcl, because CMFISHCR13.gcl doesn't contain floy tag fish genotypes, and thus are not in the .gcl object
# Get field data for CMFISHCRT13.gcl
CMFISHCRT13_field=read.xlsx(file="Data/CMFISHCRT13 field data.xlsx",sheetName="Kyle")
str(CMFISHCRT13_field)

# Tag FISH_BARCODES
as.character(CMFISHCRT13_field$FISH_CODE[!is.na(CMFISHCRT13_field$CMFISHCR13_FISH_CODE)])

# Of the floy tag recoveries, where do the paired data go?
match(as.character(CMFISHCRT13_field$FISH_CODE[!is.na(CMFISHCRT13_field$CMFISHCR13_FISH_CODE)]),CMFISHCRT13.gcl$attributes$FISH_BARCODE)

# Carcass FISH_BARCODES
as.character(CMFISHCRT13_field$CMFISHCR13_FISH_CODE[!is.na(CMFISHCRT13_field$CMFISHCR13_FISH_CODE)])

# Match carcass FISH_BARCODES with otolith data
match(as.character(CMFISHCRT13_field$CMFISHCR13_FISH_CODE[!is.na(CMFISHCRT13_field$CMFISHCR13_FISH_CODE)]),
      as.character(CMFISHCR13_oto$FISH_BARCODE))

as.character(CMFISHCR13_oto$OTOLITH_MARK_PRESENT[match(as.character(CMFISHCRT13_field$CMFISHCR13_FISH_CODE[!is.na(CMFISHCRT13_field$CMFISHCR13_FISH_CODE)]),
                                                      as.character(CMFISHCR13_oto$FISH_BARCODE))])

str(CMFISHCRT13.gcl$attributes)
CMFISHCRT13.gcl$attributes$OTOLITH_ORIGIN=NA
CMFISHCRT13.gcl$attributes$OTOLITH_ORIGIN[match(as.character(CMFISHCRT13_field$FISH_CODE[!is.na(CMFISHCRT13_field$CMFISHCR13_FISH_CODE)]),CMFISHCRT13.gcl$attributes$FISH_BARCODE)]=
  as.character(CMFISHCR13_oto$OTOLITH_MARK_PRESENT[match(as.character(CMFISHCRT13_field$CMFISHCR13_FISH_CODE[!is.na(CMFISHCRT13_field$CMFISHCR13_FISH_CODE)]),
                                                         as.character(CMFISHCR13_oto$FISH_BARCODE))])

CMFISHCRT13.gcl$attributes$OTOLITH_ORIGIN[is.na(CMFISHCRT13.gcl$attributes$OTOLITH_ORIGIN)]="U"
CMFISHCRT13.gcl$attributes$OTOLITH_ORIGIN[which(CMFISHCRT13.gcl$attributes$OTOLITH_ORIGIN=="YES")]="H"
CMFISHCRT13.gcl$attributes$OTOLITH_ORIGIN[which(CMFISHCRT13.gcl$attributes$OTOLITH_ORIGIN=="NO")]="N"
table(CMFISHCRT13.gcl$attributes$OTOLITH_ORIGIN)

#### QC =======================================================================
# QC in excel
CMFISHCRT13.gcl$attributes$OTOLITH_ORIGIN[which(CMFISHCRT13.gcl$attributes$OTOLITH_ORIGIN!="U")]
CMFISHCRT13.gcl$attributes$FISH_BARCODE[which(CMFISHCRT13.gcl$attributes$OTOLITH_ORIGIN!="U")]
# Checks out!

#### OTOLITH AGE ==============================================================
# Add OTOLITH_AGE data to CMFISHCR13.gcl attributes by matching FISH_BARCODE
CMFISHCR13_oto$OTOLITH_AGE=13-1-as.numeric(sapply(as.character(CMFISHCR13_oto$OTOLITH_MARK_ID),function(markID){unlist(strsplit(markID,split="0"))[2]})) # DOUBLE CHECK; Andrew agrees, check with Lorna
table(CMFISHCR13_oto$OTOLITH_AGE)

CMFISHCR13.gcl$attributes$OTOLITH_AGE=CMFISHCR13_oto$OTOLITH_AGE[match(CMFISHCR13.gcl$attributes$FISH_BARCODE,as.character(CMFISHCR13_oto$FISH_BARCODE))]
str(CMFISHCR13.gcl$attributes)
cbind(CMFISHCR13.gcl$attributes$OTOLITH_AGE,CMFISHCR13.gcl$attributes$OTOLITH_ORIGIN)

# Add OTOLITH_AGE data to CMFISHCRT13.gcl attributes by matching FISH_BARCODE
str(CMFISHCRT13.gcl$attributes)
CMFISHCRT13.gcl$attributes$OTOLITH_AGE=NA
CMFISHCRT13.gcl$attributes$OTOLITH_AGE[match(as.character(CMFISHCRT13_field$FISH_CODE[!is.na(CMFISHCRT13_field$CMFISHCR13_FISH_CODE)]),CMFISHCRT13.gcl$attributes$FISH_BARCODE)]=
  as.character(CMFISHCR13_oto$OTOLITH_MARK_ID[match(as.character(CMFISHCRT13_field$CMFISHCR13_FISH_CODE[!is.na(CMFISHCRT13_field$CMFISHCR13_FISH_CODE)]),
                                                         as.character(CMFISHCR13_oto$FISH_BARCODE))])
cbind(CMFISHCRT13.gcl$attributes$FISH_BARCODE,CMFISHCRT13.gcl$attributes$OTOLITH_AGE,CMFISHCRT13.gcl$attributes$OTOLITH_ORIGIN)
CMFISHCRT13.gcl$attributes$OTOLITH_AGE=13-1-as.numeric(sapply(as.character(CMFISHCRT13.gcl$attributes$OTOLITH_AGE),function(markID){unlist(strsplit(markID,split="0"))[2]}))

#### SEX ======================================================================
# Get SEX data from field data worksheets
CMFISHCRT13_field=read.xlsx(file="Data/CMFISHCRT13 field data.xlsx",sheetName="Kyle")
str(CMFISHCRT13_field)

CMFISHCR13_field=read.xlsx(file="Data/CMFISHCR13 field data.xlsx",sheetName="Kyle")
str(CMFISHCR13_field)


## Pair SEX data
# Add SEX data to .gcl attributes by matching FISH_BARCODE
match(CMFISHCRT13.gcl$attributes$FISH_BARCODE,as.character(CMFISHCRT13_field$FISH_CODE))
CMFISHCRT13.gcl$attributes$SEX=as.character(CMFISHCRT13_field$SEX[match(CMFISHCRT13.gcl$attributes$FISH_BARCODE,as.character(CMFISHCRT13_field$FISH_CODE))])

match(CMFISHCR13.gcl$attributes$FISH_BARCODE,as.character(CMFISHCR13_field$FISH_CODE))
CMFISHCR13.gcl$attributes$SEX=as.character(CMFISHCR13_field$SEX[match(CMFISHCR13.gcl$attributes$FISH_BARCODE,as.character(CMFISHCR13_field$FISH_CODE))])

#### QC =======================================================================
# QC in excel with field data sheets and otolith SQL pull
str(CMFISHCR13.gcl$attributes)
cbind(CMFISHCR13.gcl$attributes$FISH_BARCODE,CMFISHCR13.gcl$attributes$SEX,CMFISHCR13.gcl$attributes$OTOLITH_ORIGIN,CMFISHCR13.gcl$attributes$OTOLITH_AGE)

str(CMFISHCRT13.gcl$attributes)
cbind(CMFISHCRT13.gcl$attributes$FISH_BARCODE,CMFISHCRT13.gcl$attributes$SEX,CMFISHCRT13.gcl$attributes$OTOLITH_ORIGIN,CMFISHCRT13.gcl$attributes$OTOLITH_AGE)


# Brute force exercises #######################################################
# determine the number of scored loci per individual
PEDsamples
SILLY=get(paste(PEDsamples[3],".gcl",sep=''))

table(SILLY$counts[,,1]+SILLY$counts[,,2])
sum(is.na(SILLY$counts[,,1]+SILLY$counts[,,2]))

scores_fail_indv=NULL
for(i in 1:SILLY$n){
  scores_fail_indv[i]=sum(is.na(SILLY$counts[i,,1]+SILLY$counts[i,,2]))
};rm(i)
scores_indv=188-scores_fail_indv
scores_indv_percent=round(scores_indv/188*100,2)
hist(scores_indv_percent,breaks=100,col=8)
sum(scores_indv_percent<80) # CMFISHCR13 - 14 individuals; CMFISHCRT13 - 2 individuals; CMFISHCR14a - 14 individuals

# determine the number of individuals scored per loci
str(SILLY$counts)
table(SILLY$counts[,,1]+SILLY$counts[,,2])
sum(is.na(SILLY$counts[,,1]+SILLY$counts[,,2]))

scores_fail_loci=NULL
for(i in 1:188){
  scores_fail_loci[i]=sum(is.na(SILLY$counts[,i,1]+SILLY$counts[,i,2]))
};rm(i)
scores_loci=SILLY$n-scores_fail_loci
scores_loci_percent=round(scores_loci/SILLY$n*100,2)
hist(scores_loci_percent,breaks=20,col=8)
sum(scores_loci_percent<80) # 1 maker didn't do well
loci188[which(scores_loci_percent<80)] # CMFISHCR13 - 1 loci: Oke_mcfd2-86; CMFISHCRT13 - 2 loci: Oke_hmgb1-66, OKEGnRH-373; CMFISHCR14a - 8 loci: 


# Data QC/Massage for INDIVIDUALS =============================================

## Get sample size by locus
OriginalPEDSampleSizebyLocus <- SampSizeByLocus.GCL(PEDsamples,loci188)
round(100*apply(OriginalPEDSampleSizebyLocus,2,function(col){col/sampsize}),2) # Percent of samples w/in SILLY w/ genotype by locus
# Percent of individuals w/in SILLY genotyped per locus
OriginalPEDPercentbyLocus=round(100*apply(OriginalPEDSampleSizebyLocus,2,function(col){col/sampsize}),2)
hist(OriginalPEDPercentbyLocus[1,],breaks=20,xlab="% individuals genotyped",col=8)
hist(OriginalPEDPercentbyLocus[2,],breaks=20,xlab="% individuals genotyped",col=8)
hist(OriginalPEDPercentbyLocus[3,],breaks=20,xlab="% individuals genotyped",col=8) # whoa, this had some bad markers!
apply(OriginalPEDSampleSizebyLocus,1,min)/apply(OriginalPEDSampleSizebyLocus,1,max) ### 0.30 for Alevin

## Get number of individuals per silly before removing missing loci individuals
OriginalPEDColSize <- sapply(paste(PEDsamples,".gcl",sep=''), function(x) get(x)$n)

## Remove individuals with >20% missing data
PEDMissLoci <- RemoveIndMissLoci.GCL(sillyvec=PEDsamples,loci=loci188,proportion=0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSizePEDPostMissLoci <- sapply(paste(PEDsamples,".gcl",sep=''), function(x) get(x)$n)

PEDSampleSizes <- matrix(data=NA,nrow=length(PEDsamples),ncol=4,dimnames=list(PEDsamples,c("Genotyped","Missing","Duplicate","Final")))
PEDSampleSizes[,1] <- OriginalPEDColSize
PEDSampleSizes[,2] <- OriginalPEDColSize-ColSizePEDPostMissLoci

## Check within collections for duplicate individuals.
PEDDuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec=PEDsamples,loci=loci188,quantile=NULL,minproportion=0.95)

## Remove duplicate individuals
PEDRemovedDups <- RemoveDups.GCL(PEDDuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSizePEDPostDuplicate <- sapply(paste(PEDsamples,".gcl",sep=''), function(x) get(x)$n)

PEDSampleSizes[,3] <- ColSizePEDPostMissLoci-ColSizePEDPostDuplicate
PEDSampleSizes[,4] <- ColSizePEDPostDuplicate
PEDSampleSizes

write.xlsx(PEDSampleSizes,file="Data/PEDSampleSizes_QCed.xlsx")



# Data QC/Massage for LOCI ====================================================
# Remove mitochondrial markers
table(LocusControl$ploidy)

assign(paste("loci",188-table(LocusControl$ploidy)[1],sep=''),loci188[(which(LocusControl$ploidy==2))])

# Remove markers that are genotyped at less for less than 80% individuals
str(OriginalPEDPercentbyLocus)

# List of loci with < 80% genotypes by SILLY
LociRemove=apply(OriginalPEDPercentbyLocus,1,function(row){names(which(row<80))})
str(LociRemove)

# Position of loci with < 80% genotypes within loci185
LociRemoveNum=unique(sapply(unlist(LociRemove),function(locus){which(gsub(pattern="-",replacement=".",loci185)==locus)}))
185-length(LociRemoveNum)

assign(paste("loci",185-length(LociRemoveNum),sep=''),loci185[-c(LociRemoveNum)])

objects(patter="*loci")

# USE THIS MARKER SET???
loci176



# Create Input Files for Parents and Offspring ================================
#Compare sample size for collected vs. successfully genotyped and up to QC standards
sampsize;ColSizePEDPostDuplicate
sum(sampsize);sum(ColSizePEDPostDuplicate)

str(CMFISHCR13.gcl)
str(CMFISHCRT13.gcl)
str(CMFISHCR14a.gcl)

table(CMFISHCR13.gcl$attributes$CAPTURE_LOCATION)
table(CMFISHCR13.gcl$attributes$OTOLITH_ORIGIN)


#### NOTE #####
# KS 12-9-14 I created SOLOMON files for first run (non QCed genotypes and all 188 loci), but will not use SOLOMON anymore for computational and data reasons
# SOLOMON =====================================================================
### Create PARENT input files

## Old function, individual names are "SillySource"
#gcl2SOLOMON.HW(CMFISHCR13.gcl,"SOLOMON/CMFISHCR13.gcl.txt")
#gcl2SOLOMON.HW(CMFISHCRT13.gcl,"SOLOMON/CMFISHCRT13.gcl.txt")
## New function, individual names are "DNA_TRAY_CODE"_"DNA_TRAY_WELL_CODE"
gcl2SOLOMON.GCL(sillyvec=c("CMFISHCR13","CMFISHCRT13"),loci=loci188,output="SOLOMON/Parents.txt")

### Create OFFSPRING input file

#gcl2SOLOMON.HW(CMFISHCR14a.gcl,"SOLOMON/CMFISHCR14a.gcl.txt")
gcl2SOLOMON.GCL(sillyvec=c("CMFISHCR14a"),loci=loci188,output="SOLOMON/Offspring.txt")


# GO RUN IN SOLOMON!!!!


objects(pattern="*loci")

# Mykiss ======================================================================
## Can run Mykiss, but need to get sex data from HWI view in OceanAK
# "CMFISHCR13 field data.xlsx" and "CMFISHCRT13 field data.xlsx"
gcl2Mykiss.HW

str(CMFISHCR13.gcl$attributes)
str(CMFISHCRT13.gcl$attributes)


# Now we need to get into the Mykiss input file format of a Genepop file with 3 POPs: Dam, Sire, Offspring
# Test
gcl2Genepop.GCL(sillyvec=ParentsSILLY,loci=loci188,path="Mykiss/Test.txt")


########## Pool collections by SEX
##### Split each collection by SEX
### CMFISHCRT13.gcl 
# Dams?
unique(CMFISHCRT13.gcl$attributes$SEX)[]
CMFISHCRT13_Dam_IDs <- AttributesToIDs.GCL(silly="CMFISHCRT13",attribute="SEX",matching="F")

CMFISHCRT13_Dam_IDs <- list(as.numeric(na.omit(CMFISHCRT13_Dam_IDs)))
names(CMFISHCRT13_Dam_IDs) <- "CMFISHCRT13"

PoolCollections.GCL("CMFISHCRT13",loci=loci188,IDs=CMFISHCRT13_Dam_IDs,newname="CMFISHCRT13_Dam")
CMFISHCRT13_Dam.gcl$n ## 97

# Sires?
CMFISHCRT13_Sire_IDs <- AttributesToIDs.GCL(silly="CMFISHCRT13",attribute="SEX",matching="M")

CMFISHCRT13_Sire_IDs <- list(as.numeric(na.omit(CMFISHCRT13_Sire_IDs)))
names(CMFISHCRT13_Sire_IDs) <- "CMFISHCRT13"

PoolCollections.GCL("CMFISHCRT13",loci=loci188,IDs=CMFISHCRT13_Sire_IDs,newname="CMFISHCRT13_Sire")
CMFISHCRT13_Sire.gcl$n ## 142

table(CMFISHCRT13.gcl$attributes$SEX) # checks out



### CMFISHCR13.gcl 
# Dams?
unique(CMFISHCR13.gcl$attributes$SEX)[]
CMFISHCR13_Dam_IDs <- AttributesToIDs.GCL(silly="CMFISHCR13",attribute="SEX",matching="F")

CMFISHCR13_Dam_IDs <- list(as.numeric(na.omit(CMFISHCR13_Dam_IDs)))
names(CMFISHCR13_Dam_IDs) <- "CMFISHCR13"

PoolCollections.GCL("CMFISHCR13",loci=loci188,IDs=CMFISHCR13_Dam_IDs,newname="CMFISHCR13_Dam")
CMFISHCR13_Dam.gcl$n ## 347

# Sires?
CMFISHCR13_Sire_IDs <- AttributesToIDs.GCL(silly="CMFISHCR13",attribute="SEX",matching="M")

CMFISHCR13_Sire_IDs <- list(as.numeric(na.omit(CMFISHCR13_Sire_IDs)))
names(CMFISHCR13_Sire_IDs) <- "CMFISHCR13"

PoolCollections.GCL("CMFISHCR13",loci=loci188,IDs=CMFISHCR13_Sire_IDs,newname="CMFISHCR13_Sire")
CMFISHCR13_Sire.gcl$n ## 372

table(CMFISHCR13.gcl$attributes$SEX)



##### Pool two adult colletions together by SEX
#Mykiss_Names=c("Dams","Sires")
Potential_Parents=list(c("CMFISHCR13_Dam","CMFISHCRT13_Dam"),c("CMFISHCR13_Sire","CMFISHCRT13_Sire"))
#Potential_Parents_IDs=list(Dams=c(CMFISHCR13_Dam.gcl$attributes$FISH_BARCODE,CMFISHCRT13_Dam.gcl$attributes$FISH_BARCODE),Sires=c(CMFISHCR13_Sire.gcl$attributes$FISH_BARCODE,CMFISHCRT13_Sire.gcl$attributes$FISH_BARCODE))

lapply(1:2,function(x){PoolCollections.GCL(Potential_Parents[[x]], loci=loci188, IDs = NULL,newname = paste(Potential_Parents[[x]],collapse = "."))})

str(CMFISHCR13_Dam.CMFISHCRT13_Dam.gcl)
  as.character(CMFISHCR13_Dam.CMFISHCRT13_Dam.gcl$attributes$FISH_BARCODE)
str(CMFISHCR13_Sire.CMFISHCRT13_Sire.gcl)
# All attributes are now factors...kind of annoying but here we are


##### Create Mykiss input file
Mykiss_input=c("CMFISHCR13_Dam.CMFISHCRT13_Dam","CMFISHCR13_Sire.CMFISHCRT13_Sire","CMFISHCR14a")

# Run 1: Loci 188
gcl2Genepop.GCL(sillyvec=Mykiss_input,loci=loci188,path="Mykiss/Run1_input.txt")

# Run 2: Loci 185, get rid of mitochondrial markers
dir.create("Mykiss/Run2")
gcl2Genepop.GCL(sillyvec=Mykiss_input,loci=loci185,path="Mykiss/Run2/Run2_input.txt")

# Run 3: Loci 176, get rid of mitochondrial markers and other loci w/ < 80% individuals genotyped w/in SILLY
# Genotypes have been QCed!!!
dir.create("Mykiss/Run3")
gcl2Genepop.GCL(sillyvec=Mykiss_input,loci=loci176,path="Mykiss/Run3/Run3_input.txt")

# FRANz =======================================================================
##### Creat FRANz input file ==================================================
# Input parameters
dir.create("FRANz/Run14_loci160_QCed")
path <- "FRANz/Run14_loci160_QCed/Run14_input.dat"
sillyvec <- c("CMFISHCR13_Dam.CMFISHCRT13_Dam","CMFISHCR13_Sire.CMFISHCRT13_Sire","CMFISHCR14a")
loci <- dget(file = "V:/WORK/Chum/AHRP/Objects/loci163.txt")
loci <- loci[which(LocusControl$ploidy[loci] == 2)]
parents <- 1:2
offspring <- 3

nlocation <- 1
nloci <- length(loci)
delim <- "/"
title <- "Alevin"

ngenotypes <- sum(sapply(sillyvec, function(silly) {get(paste(silly, ".gcl", sep = ''))$n} ))
location <- "FishCreek"

# Make file!
file <- paste(nlocation, nloci, delim, title, collapse = ' ')
file <- rbind(file, paste(ngenotypes, location, collapse = ' '))
alleles <- LocusControl$alleles[loci]
ploidy <- LocusControl$ploidy[loci]
nalleles <- LocusControl$nalleles[loci]
maxchar <- nchar(nalleles)+1
maxchar_indv_name <- 10

for(silly in sillyvec) {    
  my.gcl <- get(paste(silly,".gcl",sep = ""), pos = 1)
  n <- my.gcl$n
    # individual name
  invisible(ifelse(silly%in%sillyvec[parents], assign("fishIDs", sapply(as.character(my.gcl$attributes$FISH_BARCODE), function(ID) {unlist(strsplit(ID, split = "130000"))[2]} )), assign("fishIDs", as.character(my.gcl$attributes$FK_FISH_ID))))
    # make it identifiable, but only 10 characters, take up extra with " "
  indv_IDs <- sapply(fishIDs, function(ID) {c(ID, rep(" ", times = maxchar_indv_name-nchar(ID)))} )
  indv_IDs <- sapply(indv_IDs, function(i) {paste(i, collapse = "")} )
    # assign 2000 to F0 generation and 2001 to F1 generation
  invisible(ifelse(silly%in%sillyvec[parents], assign("generation", rep(2000, times = n)), assign("generation", rep(2001, n))))
    # assign sex to F0 generation and "?" to F1 generation
  invisible(ifelse(length(as.character(my.gcl$attributes$SEX)) == 0, assign("sex", rep("?", times = n)), assign("sex", as.character(my.gcl$attributes$SEX))))
    # IDs to look up genotypes
  IDs <- dimnames(my.gcl$scores)[[1]]
    # raw scores
  scores <- my.gcl$scores[, loci, ]
    # genotypes in the format 00/00 where ?  =  NA
  counts <- sapply(IDs, function(ID) {paste(sapply(loci, function(locus) {paste(sapply(1:ploidy[locus], function(ploid) {ifelse(is.na(scores[ID, locus, ploid]), paste("?", collapse = ""), paste(c(rep(0, maxchar[locus]-nchar(match(scores[ID, locus, ploid], alleles[[locus]]))), match(scores[ID, locus, ploid], alleles[[locus]])), collapse = ""))} ), collapse = "/")} ), collapse = " ")} )
    # create indivdual genotype data by row  
  indvdata <- as.character(sapply(1:n, function(ID) {paste(indv_IDs[ID], "1", generation[ID], "?", sex[ID], counts[ID], collapse = " ")} ))

  file <- rbind(file,matrix(indvdata))
}

write.table(file, path, quote = FALSE, row.names = FALSE, col.names = FALSE)


# Note: 106 ["Oke_U1010-154"] is fixed and uniformative =======================
# Also, I should run 1st w/o --fullsibtest ====================================

## Script to run FRANz from the command prompt ================================
# Run 3 176 loci, QCed genotypes, most likely LOD <FRANz --Nmax 2000 --femrepro 1:1 --malerepro 1:1 --updatefreqs V:\WORK\Chum\AHRP\FRANz\Run3_loci176_QCed\Run3_input.dat>
# Run 4 176 loci, QCed genotypes, most likely LOD, fullsibs <FRANz --Nmax 2000 --femrepro 1:1 --malerepro 1:1 --updatefreqs --fullsibtest V:\WORK\Chum\AHRP\FRANz\Run3_loci176_QCed\Run3_input.dat>
# Run 5 176 loci, QCed genotypes, positive LOD <FRANz --Nmax 2000 --femrepro 1:1 --malerepro 1:1 --updatefreqs --poutformat 2 V:\WORK\Chum\AHRP\FRANz\Run3_loci176_QCed\Run3_input.dat>

# Run 6 <FRANz --Nmmax 2000 --Nfmax 2000 --femrepro 1:1 --malerepro 1:1 --typingerror 0.005 --updatefreqs --poutformat 2 --cervusgenotypeout Cervusgenotye.txt 
# --cervusoffspringout Cervusoffspring.txt V:\WORK\Chum\AHRP\FRANz\Run3_loci176_QCed\Run3_input.dat>

# Run 7 <FRANz --Nmmax 1000 --Nfmax 1000 --femrepro 1:1 --malerepro 1:1 --typingerror 0.005 --updatefreqs --poutformat 2 --cervusgenotypeout Cervusgenotye.txt 
# --cervusoffspringout Cervusoffspring.txt V:\WORK\Chum\AHRP\FRANz\Run3_loci176_QCed\Run3_input.dat>

# New =========================================================================
# Run 7alt <FRANz --Nmmax 1000 --Nfmax 1000 --femrepro 1:1 --malerepro 1:1 --typingerror 0.005 --updatefreqs --poutformat 2 --cervusgenotypeout Cervusgenotye.txt 
# --cervusoffspringout Cervusoffspring.txt --fullsibtest V:\WORK\Chum\AHRP\FRANz\Run3_loci176_QCed\Run3_input.dat>
# NO QUALITATIVE DIFFERENCE BETWEEN USING FULLSIBLING INFORMATION OR NOT
# =============================================================================

# Use system2 function to pass commands to command prompt from R!!!
setwd("V:/WORK/Chum/AHRP/FRANz/Run14_loci160_QCed")
system("FRANz --Nmmax 1000 --Nfmax 1000 --femrepro 1:1 --malerepro 1:1 --typingerror 0.005 --updatefreqs --poutformat 2 --cervusgenotypeout Cervusgenotye.txt --cervusoffspringout Cervusoffspring.txt --fullsibtest V:\\WORK\\Chum\\AHRP\\FRANz\\Run14_loci160_QCed\\Run14_input.dat") #, stdout = "V:\\WORK\\Chum\\AHRP\\FRANz\\Run3_loci176_QCed", stderr = "V:\\WORK\\Chum\\AHRP\\FRANz\\Run3_loci176_QCed"
system("FRANz --Nmmax 1000 --Nfmax 1000 --femrepro 1:1 --malerepro 1:1 --typingerror 0.005 --updatefreqs --poutformat 2 --cervusgenotypeout Cervusgenotye.txt --cervusoffspringout Cervusoffspring.txt V:\\WORK\\Chum\\AHRP\\FRANz\\Run14_loci160_QCed\\Run14_input.dat") #, stdout = "V:\\WORK\\Chum\\AHRP\\FRANz\\Run3_loci176_QCed", stderr = "V:\\WORK\\Chum\\AHRP\\FRANz\\Run3_loci176_QCed"

## Also try SOLIDMON???

### Will need to figure out how to match up the otolith read data from Eric. L. with the fish, because most will be in CMFISHCR13, but there will be the few from CMFISHCR13T
### The floy tag number is under CAPTURE_LOCATION for both CMFISHCRT13 and CMFISHCR13
### Eric L. will send back data based on barcode_position, and I will use that to match up H/N/U designation for CMFISHCR13
### For the 24 floy tag recoveries...I will match up the CMFISHCR13 (has otolith data) and CMFISHCRT (has genotypes) barcodes by the floy tag number (CAPTURE_LOCATION)
### I will likely need to do this in excel...using OceanAK



### Mykiss output from Run 3 QCed genotypes and 176 loci
## Look at distribution of the number of MM for the "best" Dam and Sire
# Single parent
numMisDam=as.numeric(readClipboard())
table(numMisDam);hist(numMisDam,col=8,breaks=50,xlim=c(0,10),ylim=c(0,200))

numMisSire=as.numeric(readClipboard())
table(numMisSire);hist(numMisSire,col=8,breaks=50,xlim=c(0,10),ylim=c(0,200))




#Parent Pair
numMisPair=as.numeric(readClipboard())
hist(numMisPair,breaks=20)





## FRANz output from run 2, trial 2 using all assignments (regardless of posterior)
parents=readClipboard()
length(unique(parents))
table(table(parents))
plot(rev(sort(table(parents))),type="h")
par(mar=c(5.1,5.1,4.1,2.1))
plot(table(table(parents)),type="h",ylab="Frequency",xlab="Number of offpsring",main="Family Size",lwd=15,cex.lab=2,cex.main=2.5)

length(parents)


CMFISHCR13.gcl$attributes$FK_FISH_ID[which(CMFISHCR13.gcl$attributes$FISH_BARCODE=="1300001389_45")]
CMFISHCR13.gcl$attributes$FK_FISH_ID[which(CMFISHCR13.gcl$attributes$FISH_BARCODE=="1300001390_29")]


#### Sibling analysis =========================================================
# FK_FISH_ID for alevin with fullsib relationship from siblings_kyle.xlsx column X (Genotype I)
offspring_w_1plus_sib=readClipboard()
# FK_FISH_ID for alevin with at least one fullsib relationship
length(unique(offspring_w_1plus_sib))

# Plot a "histogram" of the  number of alevin per sample bottle
plot(table(table(CMFISHCR14a.gcl$attributes$BOTTLE_ID)),type="h",ylab="Frequency",xlab="Number alvein / bottle",main="Alevin sampling",lwd=15,cex.lab=2,cex.main=2.5,axes=FALSE)
  axis(side=1,c(1,seq(5,25,by=5)),cex.axis=2);axis(side=2,cex.axis=2)

# Number of offspring with at least one fullsib relationship per bottle
table(CMFISHCR14a.gcl$attributes$BOTTLE_ID[(CMFISHCR14a.gcl$attributes$FK_FISH_ID%in%offspring_w_1plus_sib)])

# Names of bottles that have at least one fullsib relationship per bottle
names(table(CMFISHCR14a.gcl$attributes$BOTTLE_ID[(CMFISHCR14a.gcl$attributes$FK_FISH_ID%in%offspring_w_1plus_sib)]))

require(scales) # allows colors to be transparent

# add "bars" for the bottles that have at least one fullsib relationship
points(table(table(CMFISHCR14a.gcl$attributes$BOTTLE_ID)[names(table(CMFISHCR14a.gcl$attributes$BOTTLE_ID[(CMFISHCR14a.gcl$attributes$FK_FISH_ID%in%offspring_w_1plus_sib)]))]),type="h",lwd=15,col=alpha(2, 0.5))
# We conclude that if at least 4 alevin are sampled in a bottle, there will be at least one full sibling relationship in the bottle


# Plot Number of alevin per bottle, sorted from smallest to largest (i.e. 1-25)
plot(sort(table(CMFISHCR14a.gcl$attributes$BOTTLE_ID)),type="h",lwd=5,axes=FALSE,xlab="",ylab="Number alevin / bottle",main="Fullsiblings",cex.lab=2,cex.main=2.5)
  axis(side=2,seq(0,25,by=5),cex.axis=2);axis(side=1,at=seq(1,66,by=1),labels=rep("",66),cex.axis=2)
  text((1:66)-1,y=rep(c(-1,-2.2),33),labels=names(sort(table(CMFISHCR14a.gcl$attributes$BOTTLE_ID))),srt=45,xpd=TRUE,cex=0.7)
  mtext(side=1,"Bottle",cex=2,line=3.5)
  abline(h=0)

names(sort(table(CMFISHCR14a.gcl$attributes$BOTTLE_ID)))

names(table(CMFISHCR14a.gcl$attributes$BOTTLE_ID[(CMFISHCR14a.gcl$attributes$FK_FISH_ID%in%offspring_w_1plus_sib)]))

match(names(table(CMFISHCR14a.gcl$attributes$BOTTLE_ID[(CMFISHCR14a.gcl$attributes$FK_FISH_ID%in%offspring_w_1plus_sib)])),names(sort(table(CMFISHCR14a.gcl$attributes$BOTTLE_ID))))



# Add "bars" for the number of those alevin that are in fullsib relationships
points(match(names(table(CMFISHCR14a.gcl$attributes$BOTTLE_ID[(CMFISHCR14a.gcl$attributes$FK_FISH_ID%in%offspring_w_1plus_sib)])),names(sort(table(CMFISHCR14a.gcl$attributes$BOTTLE_ID)))),
       table(CMFISHCR14a.gcl$attributes$BOTTLE_ID[(CMFISHCR14a.gcl$attributes$FK_FISH_ID%in%offspring_w_1plus_sib)]),type="h",lwd=5,col=2)



plot(table(table(offspring_w_1plus_sib))[-length(table(table(offspring_w_1plus_sib)))],type="h")


## Read in both genotype 1 and 2 columns
genotype1=offspring_w_1plus_sib[which(offspring_w_1plus_sib!="")]
genotype2=readClipboard()
genotype2=genotype2[which(genotype2!="")]

fullsibs=c(genotype1,genotype2)

table(fullsibs)
table(table(fullsibs))
table(table(fullsibs))/(2:)

plot(table(table(fullsibs)))



# Create Pipeline for Output files to determine relationships & RRS ===========
# Mykiss ======================================================================
# Mykiss Parent Pair output calls a mismatch if either parent is missing data and the parent with data does not have the same genotype as the offspring for a given loci!!! Even when there is no mismatch!!!
# Givne this information, we will not be using Mykiss for this project!


# FRANz =======================================================================
getwd()
setwd("V:/WORK/Chum/AHRP")
objects(pattern="*.gcl") #check objects

# Specify FRANz inputs

# Locus naming issues ====
LocusControl$locusnames
LocusControl$Publishedlocusnames
LocusControl$Publishedlocusnames = gsub(pattern=" ", replacement="", LocusControl$Publishedlocusnames)
LocusControl$Publishedlocusnames[135] = NA

# Published names for all loci
ifelse(is.na(LocusControl$Publishedlocusnames),LocusControl$locusnames,LocusControl$Publishedlocusnames)
pub188 = sort(ifelse(is.na(LocusControl$Publishedlocusnames),LocusControl$locusnames,LocusControl$Publishedlocusnames))
pub185 = readClipboard()
  pub185 = gsub("1Oke", "Oke", pub185)
  pub185 = gsub("Oke_u2020-51", "Oke_U2020-51", pub185)
  pub185 = gsub("Oke_U212-87", "Oke_u212-87", pub185)

# Which of 185 from WASSIP NOT in the 188 I have?
pub185[!pub185 %in% pub188]

# Which of the 3 188 are not in the WASSIP 185?
pub188[!pub188 %in% pub185] # "Oke_Cr30"   "Oke_Cr386"  "Oke_ND3-69"

# Full 176 SNPs ====
# Establish Posterior p-value cutoff for parentage assingment
cutoff=0.80

# Read in "parentage_kyle.csv" from FRANz
# I had to add a false column as "row" in order to get it to read in properly
# out=read.csv(file="FRANz/Run3_loci176_QCed/Run 5 posLOD/parentage_kyle.csv") #old

#setwd("F:/2015 AHRP April Seattle")
setwd("FRANz/Run3_loci176_QCed/Run 7 postLOD typingerror_Nmax")
# Read in "parentage.csv" from FRANz
# Since the last column in the "parentage.csv" output file isn't named, I skip the first row and add the col.names manually, otherwise it coerced "Offspring" as the row.names
out <- read.csv(file="parentage.csv",skip=1,header=FALSE,col.names=c("Offspring","Loci.Typed","Parent.1","Loci.Typed.1","Parent.2","Loci.Typed.2","LOD","Posterior","Common.Loci.Typed","Mismatches","n_f","n_m","Pair.LOD.Parent.1","Pair.LOD.Parent.2","Posterior.Parent.1","Posterior.Parent.2","ML"), as.is=TRUE)

str(out);head(out)

out$Parent.1=as.character(out$Parent.1);out$Parent.1[which(out$Parent.1=="")]=NA
out$Parent.2=as.character(out$Parent.2);out$Parent.2[which(out$Parent.2=="")]=NA

# Number of unique parents (recaptures)
sum(!is.na(unique(c(out$Parent.1,out$Parent.2))))
# Simple Peterson Estimator of Genetic Mark-Recapture
(959*553)/67


# Mean Posterior of offspring that have >1 Parent
mean(out$Posterior[!is.na(out$Parent.1)])

# Total number of offspring assigned >=1 parent with Posterior > 0
length(out$Parent.1[!is.na(out$Parent.1)])
length(unique(out$Offspring[!is.na(out$Parent.1)]))

# Offspring assigned >=1 parent with Posterior > cutoff
which(!is.na(out$Parent.1)&out$Posterior>cutoff)
length(which(!is.na(out$Parent.1)&out$Posterior>cutoff))
length(unique(which(!is.na(out$Parent.1)&out$Posterior>cutoff)))

# Mismatches
table(out$Mismatches[which(!is.na(out$Parent.1))])
table(out$Mismatches[which(!is.na(out$Parent.1)&out$Posterior>cutoff)])

# Create data.frame for offspring assigned >=1 parent
out_assign1=out[which(!is.na(out$Parent.1)&out$Posterior>cutoff),]
str(out_assign1)

# Create data.frame for offspring assigned 2 parents
out_assign2=out_assign1[which(!is.na(out_assign1$Parent.2)),]
str(out_assign2)




# List of parents
out_assign1$Parent.1

sum(table(table(out_assign1$Parent.1)))

# Histogram of family size
par(bg=colors()[c(356)])
par(mar=c(5.1,5.1,1.1,1.1))
plot(table(table(out_assign1$Parent.1)),type="h",lwd=15,axes=FALSE,ylab="Frequency",xlab="Family size",cex.lab=2,bty="n")
axis(side=1,1:26,cex.axis=2);axis(side=2,cex.axis=2)

plot(table(c(rep(0, 959-21), table(out_assign1$Parent.1))),type="h",lwd=15,axes=FALSE,ylab="Frequency",xlab="Family size",cex.lab=2,bty="n")
axis(side=1,1:26,cex.axis=2);axis(side=2,cex.axis=2)


#### Pool all parents into 1 .gcl object for ease of computations ####
PoolCollections.GCL(collections = c("CMFISHCR13", "CMFISHCRT13"), loci = loci188, IDs = NULL, newname = "CMFISHCR13ALL")
x <- FreqPop.GCL(sillyvec = "CMFISHCR13ALL", loci = loci188)
p <- x[, , "Allele 1"] / (x[, , "Allele 1"] + x[, , "Allele 2"])
maf <- pmin(p, 1-p)
dput(x = maf, file = "V:/WORK/Chum/AHRP/Objects/Phase 2/FishCreekAllParents_MAF.txt")



ParentAttributes=rbind(CMFISHCR13_Dam.CMFISHCRT13_Dam.gcl$attributes,CMFISHCR13_Sire.CMFISHCRT13_Sire.gcl$attributes)
str(ParentAttributes)
dimnames(ParentAttributes)[[2]]

# ID Parents
sapply(as.character(ParentAttributes$FISH_BARCODE),function(indv){unlist(strsplit(indv,split="130000"))[2]})

# Match Parent ID's in data.frame "out" with IDs in ParentAttributes object
par1_matches=match(out_assign1$Parent.1,sapply(as.character(ParentAttributes$FISH_BARCODE),function(indv){unlist(strsplit(indv,split="130000"))[2]}))
par2_matches=match(out_assign1$Parent.2,sapply(as.character(ParentAttributes$FISH_BARCODE),function(indv){unlist(strsplit(indv,split="130000"))[2]}))

str(out_assign1)
### Pair data for...

# FISH_BARCODE
out_assign1$Parent.1.FISH_BARCODE=as.character(ParentAttributes$FISH_BARCODE)[par1_matches]
out_assign1$Parent.2.FISH_BARCODE=as.character(ParentAttributes$FISH_BARCODE)[par2_matches]

# SEX
out_assign1$Parent.1.SEX=as.character(ParentAttributes$SEX)[par1_matches]
out_assign1$Parent.2.SEX=as.character(ParentAttributes$SEX)[par2_matches]

# OTOLITH_ORIGIN
out_assign1$Parent.1.OTOLITH_ORIGIN=as.character(ParentAttributes$OTOLITH_ORIGIN)[par1_matches]
out_assign1$Parent.2.OTOLITH_ORIGIN=as.character(ParentAttributes$OTOLITH_ORIGIN)[par2_matches]

# OTOLITH_AGE
out_assign1$Parent.1.OTOLITH_AGE=as.character(ParentAttributes$OTOLITH_AGE)[par1_matches]
out_assign1$Parent.2.OTOLITH_AGE=as.character(ParentAttributes$OTOLITH_AGE)[par2_matches]

out_assign1[which(out_assign1$Parent.2.FISH_BARCODE=="1300001390_29"),]
out_assign1[, c(1,3,5,20:25)]

length(unique(out_assign1$Parent.1))

# BOTTLE_ID
out_assign1$BOTTLE_ID <- CMFISHCR14a.gcl$attributes$BOTTLE_ID[match(out_assign1$Offspring, CMFISHCR14a.gcl$attributes$FK_FISH_ID)]
sort(table(out_assign1$BOTTLE_ID))
sort(table(out_assign1$Parent.1))

# How many bottles per parents?
table(out_assign1$Parent.1, out_assign1$BOTTLE_ID)
table(out_assign1$Parent.1, out_assign1$BOTTLE_ID)


out_assign1$Parent.1.SEX[match(dimnames(table(out_assign1$Parent.1, out_assign1$BOTTLE_ID))[[1]], out_assign1$Parent.1)]

#### Family Size for Dams
HDam = c(table(out_assign1$Parent.1[which(out_assign1$Parent.1.OTOLITH_ORIGIN=="H" & out_assign1$Parent.1.SEX=="F")]),
         table(out_assign1$Parent.2[which(out_assign1$Parent.2.OTOLITH_ORIGIN=="H" & out_assign1$Parent.2.SEX=="F")]))

NDam = c(table(out_assign1$Parent.1[which(out_assign1$Parent.1.OTOLITH_ORIGIN=="N" & out_assign1$Parent.1.SEX=="F")]),
         table(out_assign1$Parent.2[which(out_assign1$Parent.2.OTOLITH_ORIGIN=="N" & out_assign1$Parent.2.SEX=="F")]))

UDam = c(table(out_assign1$Parent.1[which(out_assign1$Parent.1.OTOLITH_ORIGIN=="U" & out_assign1$Parent.1.SEX=="F")]),
         table(out_assign1$Parent.2[which(out_assign1$Parent.2.OTOLITH_ORIGIN=="U" & out_assign1$Parent.2.SEX=="F")]))

#### Family Size for Sires
HSire = c(table(out_assign1$Parent.1[which(out_assign1$Parent.1.OTOLITH_ORIGIN=="H" & out_assign1$Parent.1.SEX=="M")]),
          table(out_assign1$Parent.2[which(out_assign1$Parent.2.OTOLITH_ORIGIN=="H" & out_assign1$Parent.2.SEX=="M")]))

NSire = c(table(out_assign1$Parent.1[which(out_assign1$Parent.1.OTOLITH_ORIGIN=="N" & out_assign1$Parent.1.SEX=="M")]),
          table(out_assign1$Parent.2[which(out_assign1$Parent.2.OTOLITH_ORIGIN=="N" & out_assign1$Parent.2.SEX=="M")]))

USire = c(table(out_assign1$Parent.1[which(out_assign1$Parent.1.OTOLITH_ORIGIN=="U" & out_assign1$Parent.1.SEX=="M")]),
          table(out_assign1$Parent.2[which(out_assign1$Parent.2.OTOLITH_ORIGIN=="U" & out_assign1$Parent.2.SEX=="M")]))

# Histogram Female ####
par(mar=c(5.1,5.1,2.1,1.1))
hist(HDam, col = "blue", breaks = 0:26, ylim = c(0,7), main = "Female", xlab = "Family size", ylab = "Frequency", cex.lab = 2, cex.axis = 2, cex.main = 3)
hist(NDam, col = "red", add = TRUE)
hist(UDam, col = "red", add = TRUE)
legend("topright", legend = c("Natural", "Hatchery"), fill = c("red", "blue"), cex = 2, bty = "n")

# Histogram Male ####
hist(HSire, col = "blue", breaks = 0:26, ylim = c(0,7), main = "Male", xlab = "Family size", ylab = "Frequency", cex.lab = 2, cex.axis = 2, cex.main = 3)
hist(NSire, col = "red", add = TRUE, breaks = 0:26)
hist(USire, col = "grey", add = TRUE, breaks = 0:26)
hist(1, col = "black", add = TRUE, breaks = 0:26)
hist(2, col = "purple", add = TRUE, breaks = 0:26)
legend("topright", legend = c("Natural", "Hatchery", "Unknown"), fill = c("red", "blue", "grey"), cex = 2, bty = "n")


##################### Including fish that produce 0 offspring!!!! ###########################################################
#### RS of Dams
table(ParentAttributes$OTOLITH_ORIGIN,ParentAttributes$SEX)

RS_H_Dam=sum(c(table(out_assign1$Parent.1[which(out_assign1$Parent.1.OTOLITH_ORIGIN=="H" & out_assign1$Parent.1.SEX=="F")]),
               table(out_assign1$Parent.2[which(out_assign1$Parent.2.OTOLITH_ORIGIN=="H" & out_assign1$Parent.2.SEX=="F")])))/
  table(ParentAttributes$OTOLITH_ORIGIN,ParentAttributes$SEX)["H","F"]

RS_N_Dam=sum(c(table(out_assign1$Parent.1[which(out_assign1$Parent.1.OTOLITH_ORIGIN=="N" & out_assign1$Parent.1.SEX=="F")]),
               table(out_assign1$Parent.2[which(out_assign1$Parent.2.OTOLITH_ORIGIN=="N" & out_assign1$Parent.2.SEX=="F")])))/
  table(ParentAttributes$OTOLITH_ORIGIN,ParentAttributes$SEX)["N","F"]

RSS_H_Dam=RS_H_Dam/RS_N_Dam


# Look at geometric means?
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

gm_mean(HDam);gm_mean(NDam)
gm_mean(HSire);gm_mean(NSire)


# Permutation test
n.H.Dams=table(ParentAttributes$OTOLITH_ORIGIN,ParentAttributes$SEX)["H","F"]
H_Dam_fams=c(table(out_assign1$Parent.1[which(out_assign1$Parent.1.OTOLITH_ORIGIN=="H" & out_assign1$Parent.1.SEX=="F")]),
             table(out_assign1$Parent.2[which(out_assign1$Parent.2.OTOLITH_ORIGIN=="H" & out_assign1$Parent.2.SEX=="F")]))
H_Dam_fams=as.numeric(c(rep(0,(n.H.Dams-length(H_Dam_fams))),H_Dam_fams))


n.N.Dams=table(ParentAttributes$OTOLITH_ORIGIN,ParentAttributes$SEX)["N","F"]
N_Dam_fams=c(table(out_assign1$Parent.1[which(out_assign1$Parent.1.OTOLITH_ORIGIN=="N" & out_assign1$Parent.1.SEX=="F")]),
             table(out_assign1$Parent.2[which(out_assign1$Parent.2.OTOLITH_ORIGIN=="N" & out_assign1$Parent.2.SEX=="F")]))
N_Dam_fams=as.numeric(c(rep(0,(n.N.Dams-length(N_Dam_fams))),N_Dam_fams))


reps=100000
true_diff=mean(N_Dam_fams)-mean(H_Dam_fams)
Noff=replicate(reps,sum(sample(c(N_Dam_fams,H_Dam_fams),size=n.N.Dams,replace=FALSE)))
Hoff=sum(N_Dam_fams,H_Dam_fams)-Noff
rsN=Noff/n.N.Dams
rsH=Hoff/n.H.Dams
diffs=rsN-rsH
perm_1tail_pvalue=round(sum(diffs>=true_diff)/reps,4)
hist(diffs)
abline(v=true_diff)


#### RS of Sires
table(ParentAttributes$OTOLITH_ORIGIN,ParentAttributes$SEX)

RS_H_Sire=sum(c(table(out_assign1$Parent.1[which(out_assign1$Parent.1.OTOLITH_ORIGIN=="H" & out_assign1$Parent.1.SEX=="M")]),
                table(out_assign1$Parent.2[which(out_assign1$Parent.2.OTOLITH_ORIGIN=="H" & out_assign1$Parent.2.SEX=="M")])))/
  table(ParentAttributes$OTOLITH_ORIGIN,ParentAttributes$SEX)["H","M"]

RS_N_Sire=sum(c(table(out_assign1$Parent.1[which(out_assign1$Parent.1.OTOLITH_ORIGIN=="N" & out_assign1$Parent.1.SEX=="M")]),
                table(out_assign1$Parent.2[which(out_assign1$Parent.2.OTOLITH_ORIGIN=="N" & out_assign1$Parent.2.SEX=="M")])))/
  table(ParentAttributes$OTOLITH_ORIGIN,ParentAttributes$SEX)["N","M"]

RSS_H_Sire=RS_H_Sire/RS_N_Sire



# Permutation test
n.H.Sires=table(ParentAttributes$OTOLITH_ORIGIN,ParentAttributes$SEX)["H","M"]
H_Sire_fams=c(table(out_assign1$Parent.1[which(out_assign1$Parent.1.OTOLITH_ORIGIN=="H" & out_assign1$Parent.1.SEX=="M")]),
             table(out_assign1$Parent.2[which(out_assign1$Parent.2.OTOLITH_ORIGIN=="H" & out_assign1$Parent.2.SEX=="M")]))
H_Sire_fams=as.numeric(c(rep(0,(n.H.Sires-length(H_Sire_fams))),H_Sire_fams))


n.N.Sires=table(ParentAttributes$OTOLITH_ORIGIN,ParentAttributes$SEX)["N","M"]
N_Sire_fams=c(table(out_assign1$Parent.1[which(out_assign1$Parent.1.OTOLITH_ORIGIN=="N" & out_assign1$Parent.1.SEX=="M")]),
             table(out_assign1$Parent.2[which(out_assign1$Parent.2.OTOLITH_ORIGIN=="N" & out_assign1$Parent.2.SEX=="M")]))
N_Sire_fams=as.numeric(c(rep(0,(n.N.Sires-length(N_Sire_fams))),N_Sire_fams))


reps=100000
true_diff=mean(N_Sire_fams)-mean(H_Sire_fams)
Noff=replicate(reps,sum(sample(c(N_Sire_fams,H_Sire_fams),size=n.N.Sires,replace=FALSE)))
Hoff=sum(N_Sire_fams,H_Sire_fams)-Noff
rsN=Noff/n.N.Sires
rsH=Hoff/n.H.Sires
diffs=rsN-rsH
perm_1tail_pvalue=round(sum(diffs>=true_diff)/reps,4)
hist(diffs)
abline(v=true_diff)



##################### Including fish that produce >= 1 offspring!!!! ###########################################################
#### RS of Dams
table(ParentAttributes$OTOLITH_ORIGIN,ParentAttributes$SEX)

RS_H_Dam=mean(c(table(out_assign1$Parent.1[which(out_assign1$Parent.1.OTOLITH_ORIGIN=="H" & out_assign1$Parent.1.SEX=="F")]),
                table(out_assign1$Parent.2[which(out_assign1$Parent.2.OTOLITH_ORIGIN=="H" & out_assign1$Parent.2.SEX=="F")])))

RS_N_Dam=mean(c(table(out_assign1$Parent.1[which(out_assign1$Parent.1.OTOLITH_ORIGIN=="N" & out_assign1$Parent.1.SEX=="F")]),
               table(out_assign1$Parent.2[which(out_assign1$Parent.2.OTOLITH_ORIGIN=="N" & out_assign1$Parent.2.SEX=="F")])))


RSS_H_Dam=RS_H_Dam/RS_N_Dam



# Permutation test

H_Dam_fams=c(table(out_assign1$Parent.1[which(out_assign1$Parent.1.OTOLITH_ORIGIN=="H" & out_assign1$Parent.1.SEX=="F")]),
             table(out_assign1$Parent.2[which(out_assign1$Parent.2.OTOLITH_ORIGIN=="H" & out_assign1$Parent.2.SEX=="F")]))
n.H.Dams=length(H_Dam_fams)


N_Dam_fams=c(table(out_assign1$Parent.1[which(out_assign1$Parent.1.OTOLITH_ORIGIN=="N" & out_assign1$Parent.1.SEX=="F")]),
             table(out_assign1$Parent.2[which(out_assign1$Parent.2.OTOLITH_ORIGIN=="N" & out_assign1$Parent.2.SEX=="F")]))
n.N.Dams=length(N_Dam_fams)


reps=100000
true_diff=mean(N_Dam_fams)-mean(H_Dam_fams)
Noff=replicate(reps,sum(sample(c(N_Dam_fams,H_Dam_fams),size=n.N.Dams,replace=FALSE)))
Hoff=sum(N_Dam_fams,H_Dam_fams)-Noff
rsN=Noff/n.N.Dams
rsH=Hoff/n.H.Dams
diffs=rsN-rsH
perm_1tail_pvalue=round(sum(diffs>=true_diff)/reps,4)
hist(diffs)
abline(v=true_diff)


#### RS of Sires
table(ParentAttributes$OTOLITH_ORIGIN,ParentAttributes$SEX)

RS_H_Sire=mean(c(table(out_assign1$Parent.1[which(out_assign1$Parent.1.OTOLITH_ORIGIN=="H" & out_assign1$Parent.1.SEX=="M")]),
                table(out_assign1$Parent.2[which(out_assign1$Parent.2.OTOLITH_ORIGIN=="H" & out_assign1$Parent.2.SEX=="M")])))

RS_N_Sire=mean(c(table(out_assign1$Parent.1[which(out_assign1$Parent.1.OTOLITH_ORIGIN=="N" & out_assign1$Parent.1.SEX=="M")]),
                table(out_assign1$Parent.2[which(out_assign1$Parent.2.OTOLITH_ORIGIN=="N" & out_assign1$Parent.2.SEX=="M")])))


RSS_H_Sire=RS_H_Sire/RS_N_Sire



# Permutation test

H_Sire_fams=c(table(out_assign1$Parent.1[which(out_assign1$Parent.1.OTOLITH_ORIGIN=="H" & out_assign1$Parent.1.SEX=="M")]),
             table(out_assign1$Parent.2[which(out_assign1$Parent.2.OTOLITH_ORIGIN=="H" & out_assign1$Parent.2.SEX=="M")]))
n.H.Sires=length(H_Sire_fams)


N_Sire_fams=c(table(out_assign1$Parent.1[which(out_assign1$Parent.1.OTOLITH_ORIGIN=="N" & out_assign1$Parent.1.SEX=="M")]),
             table(out_assign1$Parent.2[which(out_assign1$Parent.2.OTOLITH_ORIGIN=="N" & out_assign1$Parent.2.SEX=="M")]))
n.N.Sires=length(N_Sire_fams)


reps=100000
true_diff=mean(N_Sire_fams)-mean(H_Sire_fams)
Noff=replicate(reps,sum(sample(c(N_Sire_fams,H_Sire_fams),size=n.N.Sires,replace=FALSE)))
Hoff=sum(N_Sire_fams,H_Sire_fams)-Noff
rsN=Noff/n.N.Sires
rsH=Hoff/n.H.Sires
diffs=rsN-rsH
perm_1tail_pvalue=round(sum(diffs>=true_diff)/reps,4)
hist(diffs)
abline(v=true_diff)








### "Fri Feb 20 14:38:46 2015" Kyle Shedd
### Examining the distribution of posterior values by number of mismatches
date()

chum=read.table("clipboard",header=TRUE,sep="\t",stringsAsFactors=FALSE,na.strings="")
str(chum)
by(chum[,"Posterior"], chum[,"Mismatches"], summary)
by(chum[,"Posterior"], chum[,"Mismatches"], function(mismatch) {hist(mismatch, xlim=c(0,1), ylim=c(0,length(mismatch)), breaks=20, col=8)})
by(chum[,"Posterior"], chum[,"Mismatches"], function(mismatch) {sum(mismatch<0.8)/length(mismatch)*100})
## Percent of offspring with posterior < 0.8 for a given number of mismatches
# 0 mismatches - 17.2
# 1 mismatches - 66.7
# 2 mismatches - 87.5
# 3 mismatches - 100


# Same as above, but only for alevin that have parents w/ + LOD!
by(chum[!is.na(chum$Parent.1), "Posterior"], chum[!is.na(chum$Parent.1), "Mismatches"], summary)
by(chum[!is.na(chum$Parent.1), "Posterior"], chum[!is.na(chum$Parent.1), "Mismatches"], function(mismatch) {hist(mismatch, xlim=c(0,1), ylim=c(0,length(mismatch)), breaks=20, col=8)})
by(chum[!is.na(chum$Parent.1), "Posterior"], chum[!is.na(chum$Parent.1), "Mismatches"], function(mismatch) {sum(mismatch<0.8)/length(mismatch)*100})

sum(!is.na(CMFISHCR14a.gcl$counts[232,,1]))

table(!is.na(CMFISHCR14a.gcl$counts[232,c(loci176),1]))

# Distribution of # loci genotyped
hist(apply(CMFISHCR14a.gcl$counts[,c(loci176),1], 1, function(row) {sum(!is.na(row))})/length(loci176), col=8, main="% loci genotyped", xlab="% loci genotyed", breaks=20, xlim=c(0.8,1))
hist(apply(CMFISHCR13_Dam.CMFISHCRT13_Dam.gcl$counts[,c(loci176),1], 1, function(row) {sum(!is.na(row))})/length(loci176), col=8, main="% loci genotyped", xlab="% loci genotyed", breaks=20, xlim=c(0.8,1))
hist(apply(CMFISHCR13_Sire.CMFISHCRT13_Sire.gcl$counts[,c(loci176),1], 1, function(row) {sum(!is.na(row))})/length(loci176), col=8, main="% loci genotyped", xlab="% loci genotyed", breaks=20, xlim=c(0.8,1))





##### Try with fewer SNPs ####

topscoreloci24 <- dget(file="V:/WORK/Chum/AHRP/Objects/Phase 1/topscoreloci24.txt")
topscoreloci48 <- dget(file="V:/WORK/Chum/AHRP/Objects/Phase 1/topscoreloci48.txt")
topscoreloci72 <- dget(file="V:/WORK/Chum/AHRP/Objects/Phase 1/topscoreloci72.txt")
topscoreloci96 <- dget(file="V:/WORK/Chum/AHRP/Objects/Phase 1/topscoreloci96.txt")
topscoreloci120 <- dget(file="V:/WORK/Chum/AHRP/Objects/Phase 1/topscoreloci120.txt")
topscoreloci144 <- dget(file="V:/WORK/Chum/AHRP/Objects/Phase 1/topscoreloci144.txt")

##### 144 best QCed markers ####
##### Creat FRANz input file
# Input parameters
dir.create("FRANz/Run8_loci144_QCed")
path="FRANz/Run8_loci144_QCed/Run8_input.dat"
sillyvec=c("CMFISHCR13_Dam.CMFISHCRT13_Dam","CMFISHCR13_Sire.CMFISHCRT13_Sire","CMFISHCR14a")
loci=topscoreloci144
parents=1:2
offspring=3

nlocation=1
nloci=length(loci)
delim="/"
title="Alevin"

ngenotypes=sum(sapply(sillyvec,function(silly){get(paste(silly,".gcl",sep=''))$n}))
location="FishCreek"

# Make file!
file=paste(nlocation,nloci,delim,title,collapse=' ')
file=rbind(file,paste(ngenotypes,location,collapse=' '))
alleles=LocusControl$alleles[loci]
ploidy=LocusControl$ploidy[loci]
nalleles=LocusControl$nalleles[loci]
maxchar=nchar(nalleles)+1
maxchar_indv_name=10

for(silly in sillyvec){    
  my.gcl=get(paste(silly,".gcl",sep=""),pos=1)
  n=my.gcl$n
  # individual name
  invisible(ifelse(silly%in%sillyvec[parents],assign("fishIDs",sapply(as.character(my.gcl$attributes$FISH_BARCODE),function(ID){unlist(strsplit(ID,split="130000"))[2]})),assign("fishIDs",as.character(my.gcl$attributes$FK_FISH_ID))))
  # make it identifiable, but only 10 characters, take up extra with " "
  indv_IDs=sapply(fishIDs,function(ID){c(ID,rep(" ",times=maxchar_indv_name-nchar(ID)))})
  indv_IDs=sapply(indv_IDs,function(i){paste(i,collapse="")})
  # assign 2000 to F0 generation and 2001 to F1 generation
  invisible(ifelse(silly%in%sillyvec[parents],assign("generation",rep(2000,times=n)),assign("generation",rep(2001,n))))
  # assign sex to F0 generation and "?" to F1 generation
  invisible(ifelse(length(as.character(my.gcl$attributes$SEX))==0,assign("sex",rep("?",times=n)),assign("sex",as.character(my.gcl$attributes$SEX))))
  # IDs to look up genotypes
  IDs=dimnames(my.gcl$scores)[[1]]
  # raw scores
  scores=my.gcl$scores[,loci,]
  # genotypes in the format 00/00 where ? = NA
  counts=sapply(IDs,function(ID){paste(sapply(loci,function(locus){paste(sapply(1:ploidy[locus],function(ploid){ifelse(is.na(scores[ID,locus,ploid]),paste("?",collapse=""),paste(c(rep(0,maxchar[locus]-nchar(match(scores[ID,locus,ploid],alleles[[locus]]))),match(scores[ID,locus,ploid],alleles[[locus]])),collapse=""))}),collapse="/")}),collapse=" ")})
  # create indivdual genotype data by row  
  indvdata=as.character(sapply(1:n,function(ID){paste(indv_IDs[ID], "1", generation[ID], "?", sex[ID] ,counts[ID],collapse=" ")}))
  
  file=rbind(file,matrix(indvdata))
}

write.table(file,path,quote=FALSE,row.names=FALSE,col.names=FALSE)

# Run 8 <FRANz --Nmmax 1000 --Nfmax 1000 --femrepro 1:1 --malerepro 1:1 --typingerror 0.005 --updatefreqs --poutformat 2 --cervusgenotypeout Cervusgenotye.txt 
# --cervusoffspringout Cervusoffspring.txt V:\WORK\Chum\AHRP\FRANz\Run8_loci144_QCed\Run8_input.dat>

# Use system2 function to pass commands to command prompt from R!!!
system2("FRANz --Nmmax 2000 --Nfmax 2000 --femrepro 1:1 --malerepro 1:1 --typingerror 0.005 --updatefreqs --poutformat 2 --cervusgenotypeout Cervusgenotye.txt --cervusoffspringout Cervusoffspring.txt V:\\WORK\\Chum\\AHRP\\FRANz\\Run8_loci144_QCed\\Run8_input.dat")





##### 120 best QCed markers ####
##### Creat FRANz input file
# Input parameters
dir.create("FRANz/Run9_loci120_QCed")
path="FRANz/Run9_loci120_QCed/Run9_input.dat"
sillyvec=c("CMFISHCR13_Dam.CMFISHCRT13_Dam","CMFISHCR13_Sire.CMFISHCRT13_Sire","CMFISHCR14a")
loci=topscoreloci120
parents=1:2
offspring=3

nlocation=1
nloci=length(loci)
delim="/"
title="Alevin"

ngenotypes=sum(sapply(sillyvec,function(silly){get(paste(silly,".gcl",sep=''))$n}))
location="FishCreek"

# Make file!
file=paste(nlocation,nloci,delim,title,collapse=' ')
file=rbind(file,paste(ngenotypes,location,collapse=' '))
alleles=LocusControl$alleles[loci]
ploidy=LocusControl$ploidy[loci]
nalleles=LocusControl$nalleles[loci]
maxchar=nchar(nalleles)+1
maxchar_indv_name=10

for(silly in sillyvec){    
  my.gcl=get(paste(silly,".gcl",sep=""),pos=1)
  n=my.gcl$n
  # individual name
  invisible(ifelse(silly%in%sillyvec[parents],assign("fishIDs",sapply(as.character(my.gcl$attributes$FISH_BARCODE),function(ID){unlist(strsplit(ID,split="130000"))[2]})),assign("fishIDs",as.character(my.gcl$attributes$FK_FISH_ID))))
  # make it identifiable, but only 10 characters, take up extra with " "
  indv_IDs=sapply(fishIDs,function(ID){c(ID,rep(" ",times=maxchar_indv_name-nchar(ID)))})
  indv_IDs=sapply(indv_IDs,function(i){paste(i,collapse="")})
  # assign 2000 to F0 generation and 2001 to F1 generation
  invisible(ifelse(silly%in%sillyvec[parents],assign("generation",rep(2000,times=n)),assign("generation",rep(2001,n))))
  # assign sex to F0 generation and "?" to F1 generation
  invisible(ifelse(length(as.character(my.gcl$attributes$SEX))==0,assign("sex",rep("?",times=n)),assign("sex",as.character(my.gcl$attributes$SEX))))
  # IDs to look up genotypes
  IDs=dimnames(my.gcl$scores)[[1]]
  # raw scores
  scores=my.gcl$scores[,loci,]
  # genotypes in the format 00/00 where ? = NA
  counts=sapply(IDs,function(ID){paste(sapply(loci,function(locus){paste(sapply(1:ploidy[locus],function(ploid){ifelse(is.na(scores[ID,locus,ploid]),paste("?",collapse=""),paste(c(rep(0,maxchar[locus]-nchar(match(scores[ID,locus,ploid],alleles[[locus]]))),match(scores[ID,locus,ploid],alleles[[locus]])),collapse=""))}),collapse="/")}),collapse=" ")})
  # create indivdual genotype data by row  
  indvdata=as.character(sapply(1:n,function(ID){paste(indv_IDs[ID], "1", generation[ID], "?", sex[ID] ,counts[ID],collapse=" ")}))
  
  file=rbind(file,matrix(indvdata))
}

write.table(file,path,quote=FALSE,row.names=FALSE,col.names=FALSE)

# Run 9 <FRANz --Nmmax 1000 --Nfmax 1000 --femrepro 1:1 --malerepro 1:1 --typingerror 0.005 --updatefreqs --poutformat 2 --cervusgenotypeout Cervusgenotye.txt 
# --cervusoffspringout Cervusoffspring.txt V:\WORK\Chum\AHRP\FRANz\Run9_loci120_QCed\Run9_input.dat>

# Use system2 function to pass commands to command prompt from R!!!
system2("FRANz --Nmmax 2000 --Nfmax 2000 --femrepro 1:1 --malerepro 1:1 --typingerror 0.005 --updatefreqs --poutformat 2 --cervusgenotypeout Cervusgenotye.txt --cervusoffspringout Cervusoffspring.txt V:\\WORK\\Chum\\AHRP\\FRANz\\Run9_loci120_QCed\\Run9_input.dat")






##### 96 best QCed markers ####
##### Creat FRANz input file
# Input parameters
dir.create("FRANz/Run10_loci96_QCed")
path="FRANz/Run10_loci96_QCed/Run10_input.dat"
sillyvec=c("CMFISHCR13_Dam.CMFISHCRT13_Dam","CMFISHCR13_Sire.CMFISHCRT13_Sire","CMFISHCR14a")
loci=topscoreloci96
parents=1:2
offspring=3

nlocation=1
nloci=length(loci)
delim="/"
title="Alevin"

ngenotypes=sum(sapply(sillyvec,function(silly){get(paste(silly,".gcl",sep=''))$n}))
location="FishCreek"

# Make file!
file=paste(nlocation,nloci,delim,title,collapse=' ')
file=rbind(file,paste(ngenotypes,location,collapse=' '))
alleles=LocusControl$alleles[loci]
ploidy=LocusControl$ploidy[loci]
nalleles=LocusControl$nalleles[loci]
maxchar=nchar(nalleles)+1
maxchar_indv_name=10

for(silly in sillyvec){    
  my.gcl=get(paste(silly,".gcl",sep=""),pos=1)
  n=my.gcl$n
  # individual name
  invisible(ifelse(silly%in%sillyvec[parents],assign("fishIDs",sapply(as.character(my.gcl$attributes$FISH_BARCODE),function(ID){unlist(strsplit(ID,split="130000"))[2]})),assign("fishIDs",as.character(my.gcl$attributes$FK_FISH_ID))))
  # make it identifiable, but only 10 characters, take up extra with " "
  indv_IDs=sapply(fishIDs,function(ID){c(ID,rep(" ",times=maxchar_indv_name-nchar(ID)))})
  indv_IDs=sapply(indv_IDs,function(i){paste(i,collapse="")})
  # assign 2000 to F0 generation and 2001 to F1 generation
  invisible(ifelse(silly%in%sillyvec[parents],assign("generation",rep(2000,times=n)),assign("generation",rep(2001,n))))
  # assign sex to F0 generation and "?" to F1 generation
  invisible(ifelse(length(as.character(my.gcl$attributes$SEX))==0,assign("sex",rep("?",times=n)),assign("sex",as.character(my.gcl$attributes$SEX))))
  # IDs to look up genotypes
  IDs=dimnames(my.gcl$scores)[[1]]
  # raw scores
  scores=my.gcl$scores[,loci,]
  # genotypes in the format 00/00 where ? = NA
  counts=sapply(IDs,function(ID){paste(sapply(loci,function(locus){paste(sapply(1:ploidy[locus],function(ploid){ifelse(is.na(scores[ID,locus,ploid]),paste("?",collapse=""),paste(c(rep(0,maxchar[locus]-nchar(match(scores[ID,locus,ploid],alleles[[locus]]))),match(scores[ID,locus,ploid],alleles[[locus]])),collapse=""))}),collapse="/")}),collapse=" ")})
  # create indivdual genotype data by row  
  indvdata=as.character(sapply(1:n,function(ID){paste(indv_IDs[ID], "1", generation[ID], "?", sex[ID] ,counts[ID],collapse=" ")}))
  
  file=rbind(file,matrix(indvdata))
}

write.table(file,path,quote=FALSE,row.names=FALSE,col.names=FALSE)

# Run 10 <FRANz --Nmmax 1000 --Nfmax 1000 --femrepro 1:1 --malerepro 1:1 --typingerror 0.005 --updatefreqs --poutformat 2 --cervusgenotypeout Cervusgenotye.txt 
# --cervusoffspringout Cervusoffspring.txt V:\WORK\Chum\AHRP\FRANz\Run10_loci96_QCed\Run10_input.dat>

# Use system2 function to pass commands to command prompt from R!!!
system2("FRANz --Nmmax 2000 --Nfmax 2000 --femrepro 1:1 --malerepro 1:1 --typingerror 0.005 --updatefreqs --poutformat 2 --cervusgenotypeout Cervusgenotye.txt --cervusoffspringout Cervusoffspring.txt V:\\WORK\\Chum\\AHRP\\FRANz\\Run10_loci96_QCed\\Run10_input.dat")





##### 72 best QCed markers ####
##### Creat FRANz input file
# Input parameters
dir.create("FRANz/Run11_loci72_QCed")
path="FRANz/Run11_loci72_QCed/Run11_input.dat"
sillyvec=c("CMFISHCR13_Dam.CMFISHCRT13_Dam","CMFISHCR13_Sire.CMFISHCRT13_Sire","CMFISHCR14a")
loci=topscoreloci72
parents=1:2
offspring=3

nlocation=1
nloci=length(loci)
delim="/"
title="Alevin"

ngenotypes=sum(sapply(sillyvec,function(silly){get(paste(silly,".gcl",sep=''))$n}))
location="FishCreek"

# Make file!
file=paste(nlocation,nloci,delim,title,collapse=' ')
file=rbind(file,paste(ngenotypes,location,collapse=' '))
alleles=LocusControl$alleles[loci]
ploidy=LocusControl$ploidy[loci]
nalleles=LocusControl$nalleles[loci]
maxchar=nchar(nalleles)+1
maxchar_indv_name=10

for(silly in sillyvec){    
  my.gcl=get(paste(silly,".gcl",sep=""),pos=1)
  n=my.gcl$n
  # individual name
  invisible(ifelse(silly%in%sillyvec[parents],assign("fishIDs",sapply(as.character(my.gcl$attributes$FISH_BARCODE),function(ID){unlist(strsplit(ID,split="130000"))[2]})),assign("fishIDs",as.character(my.gcl$attributes$FK_FISH_ID))))
  # make it identifiable, but only 10 characters, take up extra with " "
  indv_IDs=sapply(fishIDs,function(ID){c(ID,rep(" ",times=maxchar_indv_name-nchar(ID)))})
  indv_IDs=sapply(indv_IDs,function(i){paste(i,collapse="")})
  # assign 2000 to F0 generation and 2001 to F1 generation
  invisible(ifelse(silly%in%sillyvec[parents],assign("generation",rep(2000,times=n)),assign("generation",rep(2001,n))))
  # assign sex to F0 generation and "?" to F1 generation
  invisible(ifelse(length(as.character(my.gcl$attributes$SEX))==0,assign("sex",rep("?",times=n)),assign("sex",as.character(my.gcl$attributes$SEX))))
  # IDs to look up genotypes
  IDs=dimnames(my.gcl$scores)[[1]]
  # raw scores
  scores=my.gcl$scores[,loci,]
  # genotypes in the format 00/00 where ? = NA
  counts=sapply(IDs,function(ID){paste(sapply(loci,function(locus){paste(sapply(1:ploidy[locus],function(ploid){ifelse(is.na(scores[ID,locus,ploid]),paste("?",collapse=""),paste(c(rep(0,maxchar[locus]-nchar(match(scores[ID,locus,ploid],alleles[[locus]]))),match(scores[ID,locus,ploid],alleles[[locus]])),collapse=""))}),collapse="/")}),collapse=" ")})
  # create indivdual genotype data by row  
  indvdata=as.character(sapply(1:n,function(ID){paste(indv_IDs[ID], "1", generation[ID], "?", sex[ID] ,counts[ID],collapse=" ")}))
  
  file=rbind(file,matrix(indvdata))
}

write.table(file,path,quote=FALSE,row.names=FALSE,col.names=FALSE)

# Run 11 <FRANz --Nmmax 1000 --Nfmax 1000 --femrepro 1:1 --malerepro 1:1 --typingerror 0.005 --updatefreqs --poutformat 2 --cervusgenotypeout Cervusgenotye.txt 
# --cervusoffspringout Cervusoffspring.txt V:\WORK\Chum\AHRP\FRANz\Run11_loci72_QCed\Run11_input.dat>

# Use system2 function to pass commands to command prompt from R!!!
system2("FRANz --Nmmax 2000 --Nfmax 2000 --femrepro 1:1 --malerepro 1:1 --typingerror 0.005 --updatefreqs --poutformat 2 --cervusgenotypeout Cervusgenotye.txt --cervusoffspringout Cervusoffspring.txt V:\\WORK\\Chum\\AHRP\\FRANz\\Run11_loci72_QCed\\Run11_input.dat")





##### 48 best QCed markers ####
##### Creat FRANz input file
# Input parameters
dir.create("FRANz/Run12_loci48_QCed")
path="FRANz/Run12_loci48_QCed/Run12_input.dat"
sillyvec=c("CMFISHCR13_Dam.CMFISHCRT13_Dam","CMFISHCR13_Sire.CMFISHCRT13_Sire","CMFISHCR14a")
loci=topscoreloci48
parents=1:2
offspring=3

nlocation=1
nloci=length(loci)
delim="/"
title="Alevin"

ngenotypes=sum(sapply(sillyvec,function(silly){get(paste(silly,".gcl",sep=''))$n}))
location="FishCreek"

# Make file!
file=paste(nlocation,nloci,delim,title,collapse=' ')
file=rbind(file,paste(ngenotypes,location,collapse=' '))
alleles=LocusControl$alleles[loci]
ploidy=LocusControl$ploidy[loci]
nalleles=LocusControl$nalleles[loci]
maxchar=nchar(nalleles)+1
maxchar_indv_name=10

for(silly in sillyvec){    
  my.gcl=get(paste(silly,".gcl",sep=""),pos=1)
  n=my.gcl$n
  # individual name
  invisible(ifelse(silly%in%sillyvec[parents],assign("fishIDs",sapply(as.character(my.gcl$attributes$FISH_BARCODE),function(ID){unlist(strsplit(ID,split="130000"))[2]})),assign("fishIDs",as.character(my.gcl$attributes$FK_FISH_ID))))
  # make it identifiable, but only 10 characters, take up extra with " "
  indv_IDs=sapply(fishIDs,function(ID){c(ID,rep(" ",times=maxchar_indv_name-nchar(ID)))})
  indv_IDs=sapply(indv_IDs,function(i){paste(i,collapse="")})
  # assign 2000 to F0 generation and 2001 to F1 generation
  invisible(ifelse(silly%in%sillyvec[parents],assign("generation",rep(2000,times=n)),assign("generation",rep(2001,n))))
  # assign sex to F0 generation and "?" to F1 generation
  invisible(ifelse(length(as.character(my.gcl$attributes$SEX))==0,assign("sex",rep("?",times=n)),assign("sex",as.character(my.gcl$attributes$SEX))))
  # IDs to look up genotypes
  IDs=dimnames(my.gcl$scores)[[1]]
  # raw scores
  scores=my.gcl$scores[,loci,]
  # genotypes in the format 00/00 where ? = NA
  counts=sapply(IDs,function(ID){paste(sapply(loci,function(locus){paste(sapply(1:ploidy[locus],function(ploid){ifelse(is.na(scores[ID,locus,ploid]),paste("?",collapse=""),paste(c(rep(0,maxchar[locus]-nchar(match(scores[ID,locus,ploid],alleles[[locus]]))),match(scores[ID,locus,ploid],alleles[[locus]])),collapse=""))}),collapse="/")}),collapse=" ")})
  # create indivdual genotype data by row  
  indvdata=as.character(sapply(1:n,function(ID){paste(indv_IDs[ID], "1", generation[ID], "?", sex[ID] ,counts[ID],collapse=" ")}))
  
  file=rbind(file,matrix(indvdata))
}

write.table(file,path,quote=FALSE,row.names=FALSE,col.names=FALSE)

# Run 12 <FRANz --Nmmax 1000 --Nfmax 1000 --femrepro 1:1 --malerepro 1:1 --typingerror 0.005 --updatefreqs --poutformat 2 --cervusgenotypeout Cervusgenotye.txt 
# --cervusoffspringout Cervusoffspring.txt V:\WORK\Chum\AHRP\FRANz\Run12_loci48_QCed\Run12_input.dat>

# Use system2 function to pass commands to command prompt from R!!!
system2("FRANz --Nmmax 2000 --Nfmax 2000 --femrepro 1:1 --malerepro 1:1 --typingerror 0.005 --updatefreqs --poutformat 2 --cervusgenotypeout Cervusgenotye.txt --cervusoffspringout Cervusoffspring.txt V:\\WORK\\Chum\\AHRP\\FRANz\\Run12_loci48_QCed\\Run12_input.dat")





##### 24 best QCed markers ####
##### Creat FRANz input file
# Input parameters
dir.create("FRANz/Run13_loci24_QCed")
path="FRANz/Run13_loci24_QCed/Run13_input.dat"
sillyvec=c("CMFISHCR13_Dam.CMFISHCRT13_Dam","CMFISHCR13_Sire.CMFISHCRT13_Sire","CMFISHCR14a")
loci=topscoreloci24
parents=1:2
offspring=3

nlocation=1
nloci=length(loci)
delim="/"
title="Alevin"

ngenotypes=sum(sapply(sillyvec,function(silly){get(paste(silly,".gcl",sep=''))$n}))
location="FishCreek"

# Make file!
file=paste(nlocation,nloci,delim,title,collapse=' ')
file=rbind(file,paste(ngenotypes,location,collapse=' '))
alleles=LocusControl$alleles[loci]
ploidy=LocusControl$ploidy[loci]
nalleles=LocusControl$nalleles[loci]
maxchar=nchar(nalleles)+1
maxchar_indv_name=10

for(silly in sillyvec){    
  my.gcl=get(paste(silly,".gcl",sep=""),pos=1)
  n=my.gcl$n
  # individual name
  invisible(ifelse(silly%in%sillyvec[parents],assign("fishIDs",sapply(as.character(my.gcl$attributes$FISH_BARCODE),function(ID){unlist(strsplit(ID,split="130000"))[2]})),assign("fishIDs",as.character(my.gcl$attributes$FK_FISH_ID))))
  # make it identifiable, but only 10 characters, take up extra with " "
  indv_IDs=sapply(fishIDs,function(ID){c(ID,rep(" ",times=maxchar_indv_name-nchar(ID)))})
  indv_IDs=sapply(indv_IDs,function(i){paste(i,collapse="")})
  # assign 2000 to F0 generation and 2001 to F1 generation
  invisible(ifelse(silly%in%sillyvec[parents],assign("generation",rep(2000,times=n)),assign("generation",rep(2001,n))))
  # assign sex to F0 generation and "?" to F1 generation
  invisible(ifelse(length(as.character(my.gcl$attributes$SEX))==0,assign("sex",rep("?",times=n)),assign("sex",as.character(my.gcl$attributes$SEX))))
  # IDs to look up genotypes
  IDs=dimnames(my.gcl$scores)[[1]]
  # raw scores
  scores=my.gcl$scores[,loci,]
  # genotypes in the format 00/00 where ? = NA
  counts=sapply(IDs,function(ID){paste(sapply(loci,function(locus){paste(sapply(1:ploidy[locus],function(ploid){ifelse(is.na(scores[ID,locus,ploid]),paste("?",collapse=""),paste(c(rep(0,maxchar[locus]-nchar(match(scores[ID,locus,ploid],alleles[[locus]]))),match(scores[ID,locus,ploid],alleles[[locus]])),collapse=""))}),collapse="/")}),collapse=" ")})
  # create indivdual genotype data by row  
  indvdata=as.character(sapply(1:n,function(ID){paste(indv_IDs[ID], "1", generation[ID], "?", sex[ID] ,counts[ID],collapse=" ")}))
  
  file=rbind(file,matrix(indvdata))
}

write.table(file,path,quote=FALSE,row.names=FALSE,col.names=FALSE)

# Run 13 <FRANz --Nmmax 1000 --Nfmax 1000 --femrepro 1:1 --malerepro 1:1 --typingerror 0.005 --updatefreqs --poutformat 2 --cervusgenotypeout Cervusgenotye.txt 
# --cervusoffspringout Cervusoffspring.txt V:\WORK\Chum\AHRP\FRANz\Run13_loci24_QCed\Run13_input.dat>

# Use system2 function to pass commands to command prompt from R!!!
system2("FRANz --Nmmax 2000 --Nfmax 2000 --femrepro 1:1 --malerepro 1:1 --typingerror 0.005 --updatefreqs --poutformat 2 --cervusgenotypeout Cervusgenotye.txt --cervusoffspringout Cervusoffspring.txt V:\\WORK\\Chum\\AHRP\\FRANz\\Run13_loci24_QCed\\Run13_input.dat")






#### Pipeline for comparing the number of assignments for different marker sets ####
# Establish Posterior p-value cutoff for parentage assingment
cutoff=0.80

## 176 loci
setwd("F:/2015 AHRP April Seattle")
setwd("FRANz/Run3_loci176_QCed/Run 7 postLOD typingerror_Nmax")
# Read in "parentage.csv" from FRANz
# Since the last column in the "parentage.csv" output file isn't named, I skip the first row and add the col.names manually, otherwise it coerced "Offspring" as the row.names
out_loci176 <- read.csv(file="parentage.csv",skip=1,header=FALSE,col.names=c("Offspring","Loci.Typed","Parent.1","Loci.Typed.1","Parent.2","Loci.Typed.2","LOD","Posterior","Common.Loci.Typed","Mismatches","n_f","n_m","Pair.LOD.Parent.1","Pair.LOD.Parent.2","Posterior.Parent.1","Posterior.Parent.2","ML"), as.is=TRUE)
str(out_loci176);head(out_loci176)
out_loci176$Parent.1=as.character(out_loci176$Parent.1);out_loci176$Parent.1[which(out_loci176$Parent.1=="")]=NA
out_loci176$Parent.2=as.character(out_loci176$Parent.2);out_loci176$Parent.2[which(out_loci176$Parent.2=="")]=NA

## 144 loci
setwd("F:/2015 AHRP April Seattle")
setwd("FRANz/Run8_loci144_QCed")
# Read in "parentage.csv" from FRANz
# Since the last column in the "parentage.csv" output file isn't named, I skip the first row and add the col.names manually, otherwise it coerced "Offspring" as the row.names
out_loci144 <- read.csv(file="parentage.csv",skip=1,header=FALSE,col.names=c("Offspring","Loci.Typed","Parent.1","Loci.Typed.1","Parent.2","Loci.Typed.2","LOD","Posterior","Common.Loci.Typed","Mismatches","n_f","n_m","Pair.LOD.Parent.1","Pair.LOD.Parent.2","Posterior.Parent.1","Posterior.Parent.2","ML"))
str(out_loci144);head(out_loci144)
out_loci144$Parent.1=as.character(out_loci144$Parent.1);out_loci144$Parent.1[which(out_loci144$Parent.1=="")]=NA
out_loci144$Parent.2=as.character(out_loci144$Parent.2);out_loci144$Parent.2[which(out_loci144$Parent.2=="")]=NA

## 120 loci
setwd("F:/2015 AHRP April Seattle")
setwd("FRANz/Run9_loci120_QCed")
# Read in "parentage.csv" from FRANz
# Since the last column in the "parentage.csv" output file isn't named, I skip the first row and add the col.names manually, otherwise it coerced "Offspring" as the row.names
out_loci120 <- read.csv(file="parentage.csv",skip=1,header=FALSE,col.names=c("Offspring","Loci.Typed","Parent.1","Loci.Typed.1","Parent.2","Loci.Typed.2","LOD","Posterior","Common.Loci.Typed","Mismatches","n_f","n_m","Pair.LOD.Parent.1","Pair.LOD.Parent.2","Posterior.Parent.1","Posterior.Parent.2","ML"))
str(out_loci120);head(out_loci120)
out_loci120$Parent.1=as.character(out_loci120$Parent.1);out_loci120$Parent.1[which(out_loci120$Parent.1=="")]=NA
out_loci120$Parent.2=as.character(out_loci120$Parent.2);out_loci120$Parent.2[which(out_loci120$Parent.2=="")]=NA

## 96 loci
setwd("F:/2015 AHRP April Seattle")
setwd("FRANz/Run10_loci96_QCed")
# Read in "parentage.csv" from FRANz
# Since the last column in the "parentage.csv" output file isn't named, I skip the first row and add the col.names manually, otherwise it coerced "Offspring" as the row.names
out_loci96 <- read.csv(file="parentage.csv",skip=1,header=FALSE,col.names=c("Offspring","Loci.Typed","Parent.1","Loci.Typed.1","Parent.2","Loci.Typed.2","LOD","Posterior","Common.Loci.Typed","Mismatches","n_f","n_m","Pair.LOD.Parent.1","Pair.LOD.Parent.2","Posterior.Parent.1","Posterior.Parent.2","ML"))
str(out_loci96);head(out_loci96)
out_loci96$Parent.1=as.character(out_loci96$Parent.1);out_loci96$Parent.1[which(out_loci96$Parent.1=="")]=NA
out_loci96$Parent.2=as.character(out_loci96$Parent.2);out_loci96$Parent.2[which(out_loci96$Parent.2=="")]=NA

## 72 loci
setwd("F:/2015 AHRP April Seattle")
setwd("FRANz/Run11_loci72_QCed")
# Read in "parentage.csv" from FRANz
# Since the last column in the "parentage.csv" output file isn't named, I skip the first row and add the col.names manually, otherwise it coerced "Offspring" as the row.names
out_loci72 <- read.csv(file="parentage.csv",skip=1,header=FALSE,col.names=c("Offspring","Loci.Typed","Parent.1","Loci.Typed.1","Parent.2","Loci.Typed.2","LOD","Posterior","Common.Loci.Typed","Mismatches","n_f","n_m","Pair.LOD.Parent.1","Pair.LOD.Parent.2","Posterior.Parent.1","Posterior.Parent.2","ML"))
str(out_loci72);head(out_loci72)
out_loci72$Parent.1=as.character(out_loci72$Parent.1);out_loci72$Parent.1[which(out_loci72$Parent.1=="")]=NA
out_loci72$Parent.2=as.character(out_loci72$Parent.2);out_loci72$Parent.2[which(out_loci72$Parent.2=="")]=NA

## 48 loci
setwd("F:/2015 AHRP April Seattle")
setwd("FRANz/Run12_loci48_QCed")
# Read in "parentage.csv" from FRANz
# Since the last column in the "parentage.csv" output file isn't named, I skip the first row and add the col.names manually, otherwise it coerced "Offspring" as the row.names
out_loci48 <- read.csv(file="parentage.csv",skip=1,header=FALSE,col.names=c("Offspring","Loci.Typed","Parent.1","Loci.Typed.1","Parent.2","Loci.Typed.2","LOD","Posterior","Common.Loci.Typed","Mismatches","n_f","n_m","Pair.LOD.Parent.1","Pair.LOD.Parent.2","Posterior.Parent.1","Posterior.Parent.2","ML"))
str(out_loci48);head(out_loci48)
out_loci48$Parent.1=as.character(out_loci48$Parent.1);out_loci48$Parent.1[which(out_loci48$Parent.1=="")]=NA
out_loci48$Parent.2=as.character(out_loci48$Parent.2);out_loci48$Parent.2[which(out_loci48$Parent.2=="")]=NA


# Number of unique parents (recaptures)
# sum(!is.na(unique(c(out$Parent.1,out$Parent.2)))) # old

length(unique(out$Parent.1[!is.na(out$Parent.1)]))
length(unique(out$Parent.1[which(!is.na(out$Parent.1)&out$Posterior>cutoff)]))

# Number of unique offspring assigned >=1 parent with Posterior > 0
length(unique(out$Offspring[!is.na(out$Parent.1)]))
length(unique(out$Offspring[which(!is.na(out$Parent.1)&out$Posterior>cutoff)]))



out_loci <- c("out_loci176","out_loci144","out_loci120","out_loci96","out_loci72","out_loci48")

sapply(out_loci, function(lst) {length(unique(get(lst)$Parent.1[which(!is.na(get(lst)$Parent.1)&get(lst)$Posterior>cutoff)]))}, simplify=FALSE, USE.NAMES=TRUE)
sapply(out_loci, function(lst) {length(unique(get(lst)$Parent.2[which(!is.na(get(lst)$Parent.1)&get(lst)$Posterior>cutoff)]))}, simplify=FALSE, USE.NAMES=TRUE)
sapply(out_loci, function(lst) {length(unique(get(lst)$Offspring[which(!is.na(get(lst)$Parent.1)&get(lst)$Posterior>cutoff)]))}, simplify=FALSE, USE.NAMES=TRUE)

sapply(out_loci, function(lst) {get(lst)[which(!is.na(get(lst)$Parent.1)&get(lst)$Posterior>cutoff),]}, simplify=FALSE, USE.NAMES=TRUE)


#### Compare Marker Sets ####
# Assume that the 176 markers are golden and then determine the % assigned correctly, % unassigned, % Type I, and % Type II error

# All offspring
unique(out_loci176$Offspring)

# Which Offspring have Parent.1 with Posterior > cutoff
assignments_gold <- out_loci176[which(!is.na(out_loci176$Parent.1)&out_loci176$Posterior>cutoff),c("Offspring", "Parent.1", "Parent.2", "Posterior.Parent.1", "Posterior.Parent.2")]
assignments_loci144 <- out_loci144[which(!is.na(out_loci144$Parent.1)&out_loci144$Posterior>cutoff),c("Offspring", "Parent.1", "Parent.2", "Posterior.Parent.1", "Posterior.Parent.2")]
assignments_loci120 <- out_loci120[which(!is.na(out_loci120$Parent.1)&out_loci120$Posterior>cutoff),c("Offspring", "Parent.1", "Parent.2", "Posterior.Parent.1", "Posterior.Parent.2")]
assignments_loci96 <- out_loci96[which(!is.na(out_loci96$Parent.1)&out_loci96$Posterior>cutoff),c("Offspring", "Parent.1", "Parent.2", "Posterior.Parent.1", "Posterior.Parent.2")]
assignments_loci72 <- out_loci72[which(!is.na(out_loci72$Parent.1)&out_loci72$Posterior>cutoff),c("Offspring", "Parent.1", "Parent.2", "Posterior.Parent.1", "Posterior.Parent.2")]
assignments_loci48 <- out_loci48[which(!is.na(out_loci48$Parent.1)&out_loci48$Posterior>cutoff),c("Offspring", "Parent.1", "Parent.2", "Posterior.Parent.1", "Posterior.Parent.2")]

# Offspring in subset that occur in assignments_gold
comp144 <- cbind(assignments_gold,assignments_loci144[match(assignments_gold$Offspring,assignments_loci144$Offspring),])
comp120 <- cbind(assignments_gold,assignments_loci120[match(assignments_gold$Offspring,assignments_loci120$Offspring),])
comp96 <- cbind(assignments_gold,assignments_loci96[match(assignments_gold$Offspring,assignments_loci96$Offspring),])
comp72 <- cbind(assignments_gold,assignments_loci72[match(assignments_gold$Offspring,assignments_loci72$Offspring),])
comp48 <- cbind(assignments_gold,assignments_loci48[match(assignments_gold$Offspring,assignments_loci48$Offspring),])


## Total number of assigments possible
offspring_ID <- CMFISHCR14a.gcl$attributes$FK_FISH_ID
NTot <- length(offspring_ID)


#### Determine Error Rates in subset loci assignments compared to "gold" assignments with all 176 loci ####
comp_mat <- matrix(data=NA, nrow=NTot, ncol=5, dimnames=list(offspring_ID, c(paste("loci",c(48,72,96,120,144),sep=""))))

for (markerset in c(48,72,96,120,144)){
  
  assignments_loci <- get(paste("assignments_loci",markerset,sep=""))
  comp <- get(paste("comp",markerset,sep=""))
  
  # Dealing with extra offspring that occur in subset and NOT in assignments_gold
  add_comp <- apply(assignments_loci[!assignments_loci$Offspring%in%assignments_gold$Offspring,], 1, function(row) {
    ifelse(is.na(row[3]), "Added_Excluded", "Added_Added")
    }) # 2
  names(add_comp) <- assignments_loci[!assignments_loci$Offspring%in%assignments_gold$Offspring, "Offspring"]
  
  # Dealing with Offspring that occur in assignments_gold
  gold_comp <- apply(comp, 1, function(row) {
    ifelse(is.na(row[3])&is.na(row[8]),
           ifelse(is.na(row[7]), "Type II_Excluded",
                  ifelse(row[2]==row[7], "Correct_Excluded", "Type I_Excluded")),
           ifelse(is.na(row[3])&!is.na(row[8]),
                  ifelse(row[10] >= cutoff,
                         ifelse(sum(row[c(7,8)] %in% row[2])==1, "Correct_Added", "Type I_Added"), "Correct_Excluded"),
                  ifelse(!is.na(row[3])&is.na(row[8]),
                         ifelse(sum(row[7] %in% row[c(2,3)])==1, "Correct_Type II", "Type I_Type II"),
                         ifelse(!is.na(row[3])&!is.na(row[8]),
                                ifelse(sum(row[c(7,8)] %in% row[c(2,3)])==0, "Type I_Type I",
                                       ifelse(sum(row[c(7,8)] %in% row[c(2,3)])==1, "Correct_Type I", "Correct_Correct"))))))
  })
  
  names(gold_comp) <- comp$Offspring
  
  relate_comp <- c(add_comp,gold_comp)
  
  # Deeling with the offspring that occur in neither subset or assignments_gold (these are all "Excluded_Excluded")
  exclude_comp <- rep(x="Excluded_Excluded", times=length(offspring_ID[!offspring_ID%in%names(relate_comp)]))
  names(exclude_comp) <- offspring_ID[!offspring_ID%in%names(relate_comp)]
  
  out_comp <- c(relate_comp,exclude_comp)
  
  out_comp <- out_comp[order(as.numeric(names(out_comp)))]

  comp_mat[,paste("loci",markerset,sep="")] <- out_comp
}

str(comp_mat)
apply(comp_mat, 2, table)


# Incorporate "gold" assignments
relate_gold <- ifelse(!is.na(assignments_gold[3]), "Correct_Correct", "Correct_Excluded")
names(relate_gold) <- assignments_gold$Offspring

exclude_gold <- rep(x="Excluded_Excluded", times=length(offspring_ID[!offspring_ID%in%names(relate_gold)]))
names(exclude_gold) <- offspring_ID[!offspring_ID%in%names(relate_gold)]

out_gold <- c(relate_gold, exclude_gold)

out_gold <- out_gold[order(as.numeric(names(out_gold)))]
table(out_gold)

### FINAL Error Rate Results (for this cutoff and the FRANz input parameters modeled)
OUT <- cbind(comp_mat, "loci176"=out_gold)

apply(OUT, 2, table)

OUT

#### Plot Error Rate Results by markerset ####
names(apply(OUT, 2, table)[[1]])
classes <- unique(unlist(lapply(apply(OUT, 2, table),names)))[c(4,3,1,7,5,2,8,6)]

# Create matrix with 
OUT_mat <- matrix(data=0, nrow=length(classes), ncol=ncol(OUT), dimnames=list(classes, colnames(OUT)))
for(markerset in paste("loci",c(48,72,96,120,144,176),sep="")){
  OUT_mat[match(names(apply(OUT, 2, table)[[markerset]]), classes), markerset] = apply(OUT, 2, table)[[markerset]]
}

# With all Offspring
require(gplots)

par(mar=c(4.1,4.1,1.1,9.1))
par(bg=colors()[c(356)])

barplot2(height=OUT_mat, beside=FALSE, names.arg=c(48,72,96,120,144,176), col=c("darkblue", "darkcyan", "cyan", "lightskyblue1", "green", "palegreen3", "lightgreen", "grey"), cex.axis=1.6, cex.names=1.6)
abline(h=0); mtext(text="Number of Loci", side=1, line=3, cex=2); mtext(text="Offspring", side=2, line=2.5, cex=2)
legend(x=7.3, y=500, legend=rev(rownames(OUT_mat)), col=rev(c("darkblue", "darkcyan", "cyan", "lightskyblue1", "green", "palegreen3", "lightgreen", "grey")), pch=15, xpd=TRUE, bty="n")

barplot2(height=OUT_mat[-8,], beside=FALSE, names.arg=c(48,72,96,120,144,176), col=c("darkblue", "darkcyan", "cyan", "lightskyblue1", "green", "palegreen3", "lightgreen"), cex.axis=1.6, cex.names=1.6)
abline(h=0); mtext(text="Number of Loci", side=1, line=3, cex=2); mtext(text="Offspring", side=2, line=2.5, cex=2)
legend(x=7.3, y=99.4, legend=rev(rownames(OUT_mat)), col=rev(c("darkblue", "darkcyan", "cyan", "lightskyblue1", "green", "palegreen3", "lightgreen", "grey")), pch=15, xpd=TRUE, bty="n")
legend(x=7.3, y=100, legend=rev(rownames(OUT_mat))[-1], col=rev(c("darkblue", "darkcyan", "cyan", "lightskyblue1", "green", "palegreen3", "lightgreen")), pch=15, xpd=TRUE, bty="n")

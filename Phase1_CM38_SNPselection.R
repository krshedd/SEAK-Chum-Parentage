## 2014 CM 38, Project 37 Chum SNP selection for the AHRP
## 380 fish (95 from each of the 4 pedigree streams), all 188 available chum markers
## Kyle Shedd Tuesday August 5 2014 08:23:50

ls()
rm(list=ls(all=TRUE))
search()
getwd()
#setwd("V:/WORK/Chum/AHRP")
setwd("F:/2015 AHRP April Seattle")
## load("Phase1_CM38_SNPselection.RData")


## save.image("V:/WORK/Chum/AHRP/Phase1_CM38_SNPselection.RData")
## load("V:/WORK/Chum/AHRP/Phase1_CM38_SNPselection.RData")

#This sources all of the new GCL functions to this workspace
source("V:/DATA/R_GEN/GCL Source Scripts/Functions.GCL.R")



## Pull all data for each silly code and create .gcl objects for each
AHRPsamples <- c("CMADMCR13","CMFISHCR13","CMPROSCR13","CMSAWCR13")

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

## Get LocusControl and unaltered .gcl's to start from scratch...
LocusControl <- dget(file = "Objects/Phase 1/OriginalLocusControl.txt")

for(silly in AHRPsamples) {
  assign(x = paste(silly, ".gcl", sep = ""), 
         value = dget(file = paste("Raw genotypes/Phase 1/", silly, ".txt", sep="")))
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




#######################################################################################################################################################################################################################
PHASE 1: unranked, veto of markers not in HWE, in LD, or <80% fish genotyped (relative to best marker) #############################################################################################################
#######################################################################################################################################################################################################################

phase1=array(dim=c(188,6,6),dimnames=list(loci188,c(AHRPsamples,"Overall","Keep"),c("HWE","LD","Success Rate (%)","MAF","sdMAF","Score")))

AHRPsamples <- unlist(strsplit(ls(pattern='.gcl')[c(1:4)],split='.gcl'))

#Testing GCL functions
HoFisFstTable.GCL(sillyvec=AHRPsamples,loci=loci188,dir="V:/WORK/Chum/AHRP/Data")

###########################################################################################################################################################################################################################################
# Need to pull stray out to calculate HWE and LD...so that we are only testing the "1 population" of Fish Creek, per se
# Give Eric L. the plate and well ID info so he can dig in the database to get otolith data on stray rate
###########################################################################################################################################################################################################################################

str(CMADMCR13.gcl)


ChumSamples=data.frame()
for (collection in AHRPsamples){
  ChumSamples=rbind(ChumSamples,get(paste(collection,".gcl",sep=""))$attributes[c(4,20:22)])
}
write.xlsx(ChumSamples,"Data/CM38_SamplesforSNPselection.xlsx")

# Read in SQL query to get otolith designations, received on 8/13/14
SQL=read.xlsx("V:/WORK/Chum/AHRP/Data/CM38_oto_data_20140813_SQLpull.xlsx",1)

# Create "KEY" in attributes table for DNATRAY_WELLNO
CMADMCR13.gcl$attributes$KEY=unlist(apply(cbind(CMADMCR13.gcl[[4]][[20]],CMADMCR13.gcl[[4]][[21]]),1,function(row){paste(row[1:2],collapse="_")}))
CMFISHCR13.gcl$attributes$KEY=unlist(apply(cbind(CMFISHCR13.gcl[[4]][[20]],CMFISHCR13.gcl[[4]][[21]]),1,function(row){paste(row[1:2],collapse="_")}))
CMPROSCR13.gcl$attributes$KEY=unlist(apply(cbind(CMPROSCR13.gcl[[4]][[20]],CMPROSCR13.gcl[[4]][[21]]),1,function(row){paste(row[1:2],collapse="_")}))
CMSAWCR13.gcl$attributes$KEY=unlist(apply(cbind(CMSAWCR13.gcl[[4]][[20]],CMSAWCR13.gcl[[4]][[21]]),1,function(row){paste(row[1:2],collapse="_")}))

str(SQL)
SQL$KEY=unlist(apply(cbind(SQL[[47]],SQL[[49]]),1,function(row){paste(row[1:2],collapse="_")}))

# Add otolith present
CMADMCR13.gcl$attributes$OTOLITH_MARK_PRESENT=ifelse(SQL$OTOLITH_MARK_PRESENT[match(CMADMCR13.gcl[[4]][[25]],SQL$KEY)]=="YES",TRUE,
                                                     ifelse(SQL$OTOLITH_MARK_PRESENT[match(CMADMCR13.gcl[[4]][[25]],SQL$KEY)]=="NO",FALSE,NA))
CMFISHCR13.gcl$attributes$OTOLITH_MARK_PRESENT=ifelse(SQL$OTOLITH_MARK_PRESENT[match(CMFISHCR13.gcl[[4]][[25]],SQL$KEY)]=="YES",TRUE,
                                                      ifelse(SQL$OTOLITH_MARK_PRESENT[match(CMFISHCR13.gcl[[4]][[25]],SQL$KEY)]=="NO",FALSE,NA))
CMPROSCR13.gcl$attributes$OTOLITH_MARK_PRESENT=ifelse(SQL$OTOLITH_MARK_PRESENT[match(CMPROSCR13.gcl[[4]][[25]],SQL$KEY)]=="YES",TRUE,
                                                      ifelse(SQL$OTOLITH_MARK_PRESENT[match(CMPROSCR13.gcl[[4]][[25]],SQL$KEY)]=="NO",FALSE,NA))
CMSAWCR13.gcl$attributes$OTOLITH_MARK_PRESENT=ifelse(SQL$OTOLITH_MARK_PRESENT[match(CMSAWCR13.gcl[[4]][[25]],SQL$KEY)]=="YES",TRUE,
                                                     ifelse(SQL$OTOLITH_MARK_PRESENT[match(CMSAWCR13.gcl[[4]][[25]],SQL$KEY)]=="NO",FALSE,NA))

# Add otolith call
CMADMCR13.gcl$attributes$OTOLITH_MARK_ID=SQL$OTOLITH_MARK_ID[match(CMADMCR13.gcl[[4]][[25]],SQL$KEY)]
CMFISHCR13.gcl$attributes$OTOLITH_MARK_ID=SQL$OTOLITH_MARK_ID[match(CMFISHCR13.gcl[[4]][[25]],SQL$KEY)]
CMPROSCR13.gcl$attributes$OTOLITH_MARK_ID=SQL$OTOLITH_MARK_ID[match(CMPROSCR13.gcl[[4]][[25]],SQL$KEY)]
CMSAWCR13.gcl$attributes$OTOLITH_MARK_ID=SQL$OTOLITH_MARK_ID[match(CMSAWCR13.gcl[[4]][[25]],SQL$KEY)]

# Add scale age
CMADMCR13.gcl$attributes$SW_AGE=SQL$SW_AGE[match(CMADMCR13.gcl[[4]][[25]],SQL$KEY)]
CMFISHCR13.gcl$attributes$SW_AGE=SQL$SW_AGE[match(CMFISHCR13.gcl[[4]][[25]],SQL$KEY)]
CMPROSCR13.gcl$attributes$SW_AGE=SQL$SW_AGE[match(CMPROSCR13.gcl[[4]][[25]],SQL$KEY)]
CMSAWCR13.gcl$attributes$SW_AGE=SQL$SW_AGE[match(CMSAWCR13.gcl[[4]][[25]],SQL$KEY)]

# Add otolith age for hatchery fish
x=unlist(strsplit(as.character(CMADMCR13.gcl$attributes$OTOLITH_MARK_ID),split="DIPAC0"))
CMADMCR13.gcl$attributes$OTO_AGE=2013-(as.numeric(x[nzchar(x)])+2000)-1;rm(x)
x=unlist(strsplit(as.character(CMFISHCR13.gcl$attributes$OTOLITH_MARK_ID),split="DIPAC0"))
CMFISHCR13.gcl$attributes$OTO_AGE=2013-(as.numeric(x[nzchar(x)])+2000)-1;rm(x)
x=unlist(strsplit(as.character(CMPROSCR13.gcl$attributes$OTOLITH_MARK_ID),split="DIPAC0"))
CMPROSCR13.gcl$attributes$OTO_AGE=2013-(as.numeric(x[nzchar(x)])+2000)-1;rm(x)
x=unlist(strsplit(as.character(CMSAWCR13.gcl$attributes$OTOLITH_MARK_ID),split="DIPAC0"))
CMSAWCR13.gcl$attributes$OTO_AGE=2013-(as.numeric(x[nzchar(x)])+2000)-1;rm(x)

# Determine stray rate per stream
tab=list()
for (i in 1:length(AHRPsamples)){
  tab[[paste(AHRPsamples[i])]]=(table(get(paste(AHRPsamples[i],".gcl",sep=""))[[4]][26]))
}

# More than one way to skin a cat
#for (collection in AHRPsamples){
#  print(table(get(paste(collection,".gcl",sep=""))[[4]][26]))
#}

stryrate=lapply(tab,function(row){round(row[2]/sum(row)*100,2)})

## Create new .gcl objects w/o Hatchery Fish

# Get ID's for all confirmed non-hatchery fish (NA's for Otolith excluded)
for(collection in AHRPsamples){
  assign(paste(collection,"natIDs",sep=""),AttributesToIDs.GCL(silly=collection,attribute="OTOLITH_MARK_PRESENT",matching="FALSE"))
}

# Make into list object
for(collection in AHRPsamples){
  assign(paste(collection,"natIDs",sep=""),list(as.numeric(na.omit(get(paste(collection,"natIDs",sep=""))))))
}

# Name list object
names(CMADMCR13natIDs)="CMADMCR13"
names(CMFISHCR13natIDs)="CMFISHCR13"
names(CMPROSCR13natIDs)="CMPROSCR13"
names(CMSAWCR13natIDs)="CMSAWCR13"

# Create new .gcl with just non-hatchery (natural-origin) fish
for(collection in AHRPsamples){
  PoolCollections.GCL(collection,loci=loci188,IDs=get(paste(collection,"natIDs",sep="")),newname=paste(collection,"nat",sep=""))
}

NATsamples=unlist(strsplit(ls(pattern='.gcl')[c(2,4,6,8)],split='.gcl'))

# Check sample sizes
for(collection in NATsamples){
  print(get(paste(collection,".gcl",sep=""))$n)
}
# Compare with attributes
lapply(tab,function(row){row[1]}) # Good ADM=90, FISH=35, PROS=66, SAW=58




###########################################################################################################################################################
Data QC/Massage
###########################################################################################################################################################

## Get sample size by locus
OriginalNATSampleSizebyLocus <- SampSizeByLocus.GCL(NATsamples,loci188)
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

write.xlsx(NATSampleSizes,file="Data/NATSampleSizes.xlsx")


###########################################################################################################################################################################################################################################

## HWE ########################################################################################
## Calculate HWE only with Natural-origin fish ################################################

gcl2Genepop.GCL(sillyvec=NATsamples,loci=loci188,path="Data/NATsamples.gen")

#library(adegenet)
#genind=read.genepop(file="V:/WORK/Chum/AHRP/Data/AHRPsamples.gen")
#HWEmatrix <- HWE.test.genind(genind,pop=NULL,permut=T,nsim=100000,hide.NA=F,res.type="matrix")
#t(HWEmatrix)


#phase1[,1:4,1]=t(HWEmatrix)
### Fisher's method for combined probabilities, note 2*4 (4 for 4 populations)
#phase1[,5,1]=t(t(round(pchisq(-2*apply(log(phase1[,1:4,1]),1,function(x) sum(x,na.rm=TRUE)),2*4,lower.tail=FALSE),3)))

#sum(t(t(sort(phase1[,5,1])))<0.05) #13
#sum(t(t(sort(phase1[,1,1])))<0.05) #11
#sum(t(t(sort(phase1[,2,1])))<0.05) #10
#sum(t(t(sort(phase1[,3,1])))<0.05) #11
#sum(t(t(sort(phase1[,4,1])))<0.05) #19

#pchisq(20.46398,2*4,lower.tail=FALSE)

###########################################################################################################################################################################################################################################
# NOTE: DO NOT USE ADEGENET FOR HWE,  USE GENEPOP
# Calculate HWE only with Natural-origin fish
# From Tech Doc 2, a marker fails if 1) overall p<0.01, or 2) any pop p<0.05)
###########################################################################################################################################################################################################################################



source("V:/DATA/R_GEN/tempGCL Source Scripts/ReadGenepopHWEKS.GCL.R")
HWE=ReadGenepopHWE.GCL(file="Data/NATsamples.txt.P")
str(HWE)
HWE[[2]]
table(HWE[[2]][,5]<0.01)
HWE[[2]][which(HWE[[2]][,5]<0.01), ]
rownames(HWE[[2]])[which(HWE[[2]][,5]<0.01)]

table(HWE[[2]][,1:4]<0.05)
rownames(HWE[[2]])[which(HWE[[2]][,1]<0.05)]
rownames(HWE[[2]])[which(HWE[[2]][,2]<0.05)]
rownames(HWE[[2]])[which(HWE[[2]][,3]<0.05)]
rownames(HWE[[2]])[which(HWE[[2]][,4]<0.05)]

t(t(apply(HWE[[2]][,1:4]<0.05,1,function(row){sum(row)})))
table(apply(HWE[[2]][,1:4]<0.05,1,function(row){sum(row)}))

failHWE=list(Overall=rownames(HWE[[2]])[which(HWE[[2]][,5]<0.01)],
             ADM=rownames(HWE[[2]])[which(HWE[[2]][,1]<0.05)],
             FISH=rownames(HWE[[2]])[which(HWE[[2]][,2]<0.05)],
             PROS=rownames(HWE[[2]])[which(HWE[[2]][,3]<0.05)],
             SAW=rownames(HWE[[2]])[which(HWE[[2]][,4]<0.05)])
failedHWE=unique(unlist(failHWE))

str(phase1)
phase1[,1:5,"HWE"]=HWE[[2]]
phase1[,6,"HWE"]=rownames(phase1)%in%failedHWE # 0 = keep, 1 = failed (toss)
phase1[,,1]
sum(phase1[,6,1]) # lost 21 markers due to HWE

# Any loci have 3 or more pops < 0.05? ########################################
names(which(apply(HWE[[2]][, c(1:4)], 1, function(row) {sum(row < 0.05)}) >= 3, useNames = TRUE))
names(which(HWE[[2]][, 5] < 0.1))
# USE THIS ####

keeploci=rownames(phase1)[which(phase1[,6,1]==0)]

## Check HWE with all fish (natural and hatchery)
gcl2Genepop.GCL(sillyvec=AHRPsamples,loci=loci188,path="Data/AHRPsamples.gen")
HWEall=ReadGenepopHWE.GCL(file="Data/AHRPsamples.txt.P")
failHWEall=list(Overall=rownames(HWEall[[2]])[which(HWEall[[2]][,5]<0.01)],
             ADM=rownames(HWEall[[2]])[which(HWEall[[2]][,1]<0.05)],
             FISH=rownames(HWEall[[2]])[which(HWEall[[2]][,2]<0.05)],
             PROS=rownames(HWEall[[2]])[which(HWEall[[2]][,3]<0.05)],
             SAW=rownames(HWEall[[2]])[which(HWEall[[2]][,4]<0.05)])

## LD #########################################################################################
## Calculate LD only with Natural-origin fish #################################################

#LinkageDisequilibriumForSNPs.GCL(sillyvec=AHRPsamples,loci=loci188) # incorrect dimensions
#AllPossiblePhenotypes.GCL(loci188) # also fails
LocusControl$alleles[loci188]
table(LocusControl$ploidy[loci188])

# Do this in Genepop and use Andy's function to read in results
source("V:/DATA/R_GEN/tempGCL Source Scripts/ReadGenepopDis.GCL.R")
LD=ReadGenepopDis.GCL(file="V:/EXEC/GenepopV4/AHRPsamples_original.txt.DIS") # only reads pops, not Fisher's combined
str(LD)
head(LD);tail(LD)

# My adaptation to Andy's function, gives a data.frame of p-values for Pop and Overall
source("V:/DATA/R_GEN/tempGCL Source Scripts/ReadGenepopDisKS.GCL.R")
LD2=ReadGenepopDis.GCL(file="Data/NATsamples.txt.DIS") 
str(LD2)
head(LD2);tail(LD2)

# How many pops in LD for each pair
LD2$fail=apply(LD2[,3:6]<0.05,1,function(row){sum(row)})
table(LD2$fail)

# Bonferroni correct
table(LD2[,"Overall"]<(0.05/dim(LD2)[1]))

# Which markers were involved in any LD for each pop
phase1[,1,"LD"]=rownames(phase1)%in%unique(c(as.character(LD2[which(LD2$CMADMCR13<0.05),1]),as.character(LD2[which(LD2$CMADMCR13<0.05),2]))) # 0 = not in any LD, 1 = involved in some LD
phase1[,2,"LD"]=rownames(phase1)%in%unique(c(as.character(LD2[which(LD2$CMFISHCR13<0.05),1]),as.character(LD2[which(LD2$CMFISHCR13<0.05),2]))) # 0 = not in any LD, 1 = involved in some LD
phase1[,3,"LD"]=rownames(phase1)%in%unique(c(as.character(LD2[which(LD2$CMPROSCR13<0.05),1]),as.character(LD2[which(LD2$CMPROSCR13<0.05),2]))) # 0 = not in any LD, 1 = involved in some LD
phase1[,4,"LD"]=rownames(phase1)%in%unique(c(as.character(LD2[which(LD2$CMSAWCR13<0.05),1]),as.character(LD2[which(LD2$CMSAWCR13<0.05),2]))) # 0 = not in any LD, 1 = involved in some LD
phase1[,5,"LD"]=rownames(phase1)%in%unique(c(as.character(LD2[which(LD2$Overall<0.05),1]),as.character(LD2[which(LD2$Overall<0.05),2]))) # 0 = not in any LD, 1 = involved in some LD
phase1[,,2]

# How many pairs in LD for 3 pops
failat3=LD2[which(LD2$fail>=3),]
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

x=FreqPop.GCL(sillyvec=AHRPsamples,loci=loci188)
str(x)
p=x[,,1]/(x[,,1]+x[,,2])

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
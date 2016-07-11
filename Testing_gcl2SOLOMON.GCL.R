setwd("V:/WORK/Chum/AHRP")
#CMFISHCR13.gcl=dget(file="Raw genotypes/Phase 2/CMFISHCR13.txt")
#CMFISHCRT13.gcl=dget(file="Raw genotypes/Phase 2/CMFISHCRT13.txt")
## load("V:/WORK/Chum/AHRP/Phase2_CM40_SNPselection_AlevinParentage.RData")


LocusControl=dget(file="Objects/Phase 2/OriginalLocusControl.txt")
loci188=LocusControl$locusnames





sillyvec=c("CMFISHCR13","CMFISHCRT13")
loci=loci188
output="SOLOMON/Parents.txt"

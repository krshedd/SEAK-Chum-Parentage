#### Looking at RAD SNPs from Seeb lab for lower 48 chum

setwd("V:/WORK/Chum/AHRP")

## save.image("V:/WORK/Chum/AHRP/RAD_SNPs.RData")
## load("V:/WORK/Chum/AHRP/RAD_SNPs.RData")
## rm(list=setdiff(ls(), c(objects(pattern="RAD"),"fun","hets","ps","phase1","allSNPs","SNPs192")))

RAD_part1=readClipboard()
RAD_part2=readClipboard()
allRAD=c(RAD_part1,RAD_part2)
RAD=matrix(allRAD,nrow=length(allRAD)/4,ncol=4,byrow=TRUE)
RAD2=RAD[,2:4]
rownames(RAD2)=RAD[,1]
colnames(RAD2)=c("Fst","Ho","He")
RAD2
RAD3=data.matrix(data.frame(RAD2,stringsAsFactors=FALSE))


fun=function(p) {2*p-2*p^2}
fun(0.2)

sum(RAD3[,"He"]>fun(0.2))
dim(RAD3)

writeClipboard(rownames(RAD3)[which(RAD3[,"He"]>fun(0.2))])

RAD3[which(RAD3[,"He"]>fun(0.2)),]


sum(RAD3[,"He"]>fun(0.2))
sum(RAD3[,"He"]>fun(0.25))
sum(RAD3[,"He"]>fun(0.3))
sum(RAD3[,"He"]>fun(0.35))
sum(RAD3[,"He"]>fun(0.4))
sum(RAD3[,"He"]>fun(0.45))

hets <- seq(0,0.5,by=0.001)
ps <- fun(hets)
names(ps) <- hets

RAD3 <- cbind(RAD3,"MAF"=sapply(RAD3[,"He"], function(x) {as.numeric(names(ps))[which.min(abs(x-ps))]}))

hist(RAD3[,"MAF"],col=8,xlim=c(0,0.5))

sum(RAD3[,"MAF"]>0.2)

names(which(RAD3[,"MAF"]>0.2))


### Compare 188 w/ these new RAD SNPs
t(t(sort(phase1[,"Overall","MAF"])))
sum(phase1[,"Overall","MAF"]>0.2)
## All Chum SNPs, assuming that the MAFs observed down south are similar to in SEAK
allSNPs <- c(phase1[,"Overall","MAF"],RAD3[,5])

length(allSNPs)

hist(allSNPs)
sum(allSNPs>0.2)
t(t(sort(allSNPs,decreasing=TRUE)))
## Best 192 SNPs ranked by MAF with WASSIP 188 + RAD 72
SNPs192 <- t(t(t(t(sort(allSNPs,decreasing=TRUE)))[1:192,]))

unlist(dimnames(SNPs192))
# Show MAF of all "best" 192 SNP panel
hist(SNPs192, col=8, xlim=c(0,0.5))
# Blue ones are the "new" RAD SNPs
hist(SNPs192[grep(pattern="RAD",unlist(dimnames(SNPs192))),],col=4,add=TRUE)
legend("topleft",legend="RAD",col=4,pch=15,bty="n",cex=2)




par(mar=c(5.1,5.1,4.1,2.1))
par(bg=colors()[c(356)])
## All 260
hist(allSNPs, col="gray70", main="Available Chum SNPs", cex.lab=2, cex.axis=2, cex.main=2)
# Which are RAD
hist(allSNPs[grep(pattern="RAD",names(allSNPs))],col="lightblue",add=TRUE)
# Best 192 SNPs
hist(SNPs192, col="gray30", add=TRUE)
# Which are RAD
hist(SNPs192[grep(pattern="RAD",unlist(dimnames(SNPs192))),],col="blue",add=TRUE)
abline(v=min(SNPs192), lwd=10)
text(x=min(SNPs192), y=32, labels="Use These", cex=3, pos=4)
legend("topright",legend=c("WASSIP 188", "RAD"),col=c("gray70","blue"),pch=15,bty="n",cex=2, pt.cex=4)

abline(v=sort(phase1[,"Overall","MAF"], decreasing=TRUE)[96], lwd=10)
text(x=min(SNPs192), y=32, labels="Use These", cex=3, pos=4)

legend("topright",legend=c("WASSIP 188", "RAD"),col=c("gray70","lightblue"),pch=15,bty="n",cex=2, pt.cex=4)




## All, but this will show RAD
hist(allSNPs, col="lightblue", main="Available Chum SNPs", cex.lab=2, cex.axis=2, cex.main=2, ylim=c(0,45), xlab="Minor Allele Frequency (MAF)")
legend("topright", legend=c("WASSIP 188", "Best 96", "RAD"), col=c("gray70", "gray30", "lightblue"), pch=15, bty="n", cex=2, pt.cex=4)


## WASSIP 188
hist(phase1[,"Overall","MAF"], col="gray70", main="Available Chum SNPs", cex.lab=2, cex.axis=2, cex.main=2, ylim=c(0,45), xlab="Minor Allele Frequency (MAF)", add=TRUE)
#hist(SNPs192[-grep(pattern="RAD",unlist(dimnames(SNPs192))),], col="gray30", add=TRUE)
abline(v=min(SNPs192), lwd=10)
text(x=min(SNPs192), y=30.5, labels="Use These", cex=3, pos=4)

legend("topright", legend=c("WASSIP 188"), col=c("gray70"), pch=15, bty="n", cex=2, pt.cex=4)
# Best 96
hist(sort(phase1[,"Overall","MAF"], decreasing=TRUE)[1:96], col="gray30", main="Available Chum SNPs", cex.lab=2, cex.axis=2, cex.main=2, ylim=c(0,45), xlab="Minor Allele Frequency (MAF)", add=TRUE)
legend("topright", legend=c("WASSIP 188", "Best 96"), col=c("gray70", "gray30"), pch=15, bty="n", cex=2, pt.cex=4)





# How many of these 192 are RAD?
names(SNPs192[grep(pattern="RAD",unlist(dimnames(SNPs192))),])
length(names(SNPs192[grep(pattern="RAD",unlist(dimnames(SNPs192))),])) # 51

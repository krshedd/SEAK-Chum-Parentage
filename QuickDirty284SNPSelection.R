rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures")
source("H:/R Source Scripts/Functions.GCL_KS.R")
source("~/../R/Functions.GCL.R")
CreateLocusControl.GCL(markersuite = "ChumParentage2015_284SNPs", username = "krshedd", password = "kokanee")
LOKI2R.GCL(sillyvec = "CMFISHCR13", username = "krshedd", password = "kokanee")
FreqPop.GCL(sillyvec = "CMFISHCR13", loci = LocusControl$locusnames)
x <- FreqPop.GCL(sillyvec = "CMFISHCR13", loci = LocusControl$locusnames)
str(x)
p <- x[, , 1] / (x[, , 1], x[, , 2])
p <- x[, , 1] / (x[, , 1] + x[, , 2])
p
hist(p)
q <- sapply(LocusControl$locusnames, function(locus) {pmin(x[1, locus, ]) / sum(x[1, locus, ])})
q
q <- sapply(LocusControl$locusnames, function(locus) {pmin(as.vector(x[1, locus, ])) / sum(x[1, locus, ])})
q
q <- sapply(LocusControl$locusnames, function(locus) {pmin(as.vector(x[1, locus, ]))})
q
sapply(LocusControl$locusnames, function(locus) {as.vector(x[1, locus, ])})
sapply(LocusControl$locusnames, function(locus) {pmin(x[1, locus, 1], x[1, locus, 2])})
sapply(LocusControl$locusnames, function(locus) {pmin(x[1, locus, 1], x[1, locus, 2]) / sum(x[1, locus, ])})
q <- sapply(LocusControl$locusnames, function(locus) {pmin(x[1, locus, 1], x[1, locus, 2]) / sum(x[1, locus, ])})
hist(q)
length(q)
length(p)
q["Oke_13594_MT"]
LocusControlOLD <- LocusControl
rm(LocusControl)
CreateLocusControl.GCL(markersuite = "Chum WASC96", username = "krshedd", password = "Ski10er!")
loci96WASC <- LocusControl$locusnames
q[loci96WASC]
length(q[loci96WASC])
hist(q[loci96WASC])
q[loci96WASC] > 0.1
names(q[loci96WASC] > 0.1)
names(q)[q[loci96WASC] > 0.1]
names(q[loci96WASC])[q[loci96WASC] > 0.1]
q[loci96WASC]
is.na(q[loci96WASC])
names(q[loci96WASC])is.na(q[loci96WASC])
names(q[loci96WASC])[is.na(q[loci96WASC])]
loci2remove <- names(q[loci96WASC])[is.na(q[loci96WASC])]
loci96WASC[loci2remove]
loci2remove <- which(q[loci96WASC])[is.na(q[loci96WASC])]
is.na(q[loci96WASC])
which(is.na(q[loci96WASC]))
loci2remove <- which(is.na(q[loci96WASC]))
length(loci2remove)
loci89WASC <- loci96WASC[loci2remove]
loci89WASC
loci89WASC <- loci96WASC[-loci2remove]
loci89WASC
q[loci89WASC]
hist(q[loci89WASC])
hist(q[loci89WASC], col = 8)
which(q[loci89WASC > 0.1])
which(q[loci89WASC] > 0.1)
loci89WASC[which(q[loci89WASC] > 0.1)]
loci89WASC[which(q[loci89WASC] > 0.05)]
loci284 <- LocusControlOLD$locusnames
loci284[which(q[loci284] > 0.05)]
sort(q)
sort(q, decreasing = TRUE)
names(sort(q, decreasing = TRUE)[1])
names(sort(q, decreasing = TRUE)[5])
names(sort(q, decreasing = TRUE)[1:192])
sort(q, decreasing = TRUE)[1:192]
summary(sort(q, decreasing = TRUE)[1:192])

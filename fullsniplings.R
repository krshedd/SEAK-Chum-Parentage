# Testing Eric Anderson's "fullsniplings" program; 2014 Fish Creek alevin data
# Kyle Shedd Wed May 06 11:03:47 2015
# =============================================================================

ls()
rm(list=ls(all=TRUE))
search()
getwd()
setwd("V:/WORK/Chum/AHRP")

## save.image("fullsniplings.RData")
## load("fullsniplings.RData")


require(Rcpp)
require(devtools)
require(doParallel)

# download, build and install the package
devtools::install_github(repo="fullsniplings", username="eriqande")

require(fullsniplings)

# example Chinook dataset
str(fs_dev_test_data$chinook_full_sibs_genos)
head(fs_dev_test_data$chinook_full_sibs_genos)
rownames(fs_dev_test_data$chinook_full_sibs_genos)

# main function
?run_mcmc

# test run with example Chinook dataset
quick_run  <- run_mcmc(fs_dev_test_data$chinook_full_sibs_genos, burn_in = 50, num_sweeps = 200)
str(quick_run$Partition)
str(quick_run$GP_result$VSL)
length(quick_run$VS)
quick_run$SibGroupsByName

# Plot results from sibship
plot(quick_run$Partition$Visits, type = "p", pch = 16, cex = quick_run$Partition$NumSibs)

# length of run used by Eric A. in his TBP paper
full_run  <- run_mcmc(fs_dev_test_data$chinook_full_sibs_genos, burn_in = 100, num_sweeps = 500)

# Check out test dataset to see how Eric A. did the analyses with fewer SNPs (wer they ranked?)
odd <- seq(from = 1, to = 190, by = 2)
even <- seq(from = 2, to = 190, by = 2)

p <- sapply(1:95, function(loci) {table(c(fs_dev_test_data$chinook_full_sibs_genos[, odd[loci]], fs_dev_test_data$chinook_full_sibs_genos[, even[loci]]))[1] /
                             length(!is.na(c(fs_dev_test_data$chinook_full_sibs_genos[, odd[loci]], fs_dev_test_data$chinook_full_sibs_genos[, even[loci]])))})
MAF <- pmin(p, 1-p)
hist(MAF)
plot(MAF, type = "o", pch = 16, xlab = "SNP")

## Create an input file from 2014 Fish Creek Chum alevin (i.e. dataframe) =====
# nrow = number of alevin genotyped
# ncol = 2 * number of SNPs genotyped
# values are integer values of genotypes (A = 1, C = 2, G = 3, T = 4)

setwd("V:/WORK/Chum/AHRP/Objects")
dput(x = CMFISHCR14a.gcl, file = "CMFISHCR14a.gcl.txt")
dput(x = loci176, file = "loci176.txt")

# Load Fish Creek 2014 alevin genotypes and loci176
CMFISHCR14a.gcl <- dget(file = "CMFISHCR14a.gcl.txt")
loci176 <- dget(file = "loci176.txt")

str(CMFISH14a.gcl)
str(CMFISH14a.gcl$scores)

# Subset genotype scores into loci176
alevin_genos_array <- CMFISH14a.gcl$scores[, loci176, ]
str(alevin_genos_array)

# Create a "two column" format dataframe of genotypes
alevin_genos_mat <- matrix(data = NA, nrow = dim(alevin_genos_array)[1], ncol = dim(alevin_genos_array)[2] * dim(alevin_genos_array)[3])
for(locus in 1:dim(alevin_genos_array)[2]){
  alevin_genos_mat[, (locus * 2 - 1):(locus * 2)] <- alevin_genos_array[, locus, ]
}

head(alevin_genos_mat)
rownames(alevin_genos_mat) <- CMFISH14a.gcl$attributes$FK_FISH_ID
colnames(alevin_genos_mat) <- paste(rep(loci176, each = 2), 1:2, sep = "_")

# Dataframe
alevin_genos_mat_num <- gsub(pattern = "A", replacement = 1, x = alevin_genos_mat)
alevin_genos_mat_num <- gsub(pattern = "C", replacement = 2, x = alevin_genos_mat_num)
alevin_genos_mat_num <- gsub(pattern = "G", replacement = 3, x = alevin_genos_mat_num)
alevin_genos_mat_num <- gsub(pattern = "T", replacement = 4, x = alevin_genos_mat_num)
str(alevin_genos_mat_num)
class(alevin_genos_mat_num) <- "integer"

alevin_genos_df <- as.data.frame(alevin_genos_mat_num)
str(alevin_genos_df)

# Run fullsniplings ===========================================================
full_run  <- run_mcmc(alevin_genos_df, burn_in = 100, num_sweeps = 500)

# Base plot
plot(full_run$Partition$Visits, type = "p", pch = 16, cex = full_run$Partition$NumSibs)

# Eric Anderson's ggplot from TPB publication
lscale <- 0.02
out_df <- data.frame(Index = seq(dim(full_run$Partition)[1]), Posterior = full_run$Partition$Visits / max(full_run$Partition$Visits), NumSibs = full_run$Partition$NumSibs)
out_df$ylo <- out_df$Posterior - lscale * out_df$NumSibs
out_df$yhi <- out_df$Posterior + lscale * out_df$NumSibs

require(ggplot2)
require(plyr)

bugle_theme <- function (base_size = 12, base_family = "")
{
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(axis.text = element_text(size = rel(0.8)),
          axis.ticks = element_line(colour = "black"),
          legend.key = element_rect(colour = "grey50"),
          legend.position = "bottom",
          panel.background = element_rect(fill = "white",
                                          colour = NA),
          panel.border = element_rect(fill = NA,
                                      colour = "grey50"),
          panel.grid.major = element_line(colour = "grey98",
                                          size = 0.2),
          panel.grid.minor = element_line(colour = "grey98",
                                          size = 0.5),
          strip.background = element_rect(fill = "grey80",
                                          colour = "grey50", size = 0.2)
    )
}

ggplot(data = out_df, aes(x = Index, y = Posterior)) +
  geom_linerange(aes(ymin = ylo,
                     ymax = yhi),
                 size = 0.1) +
  geom_line(size = 0.15) +
  bugle_theme()

# How does this line up to sampling? Are all alevin in a bottle siblings? =====
table(CMFISH14a.gcl$attributes$BOTTLE_ID)
sort(table(CMFISH14a.gcl$attributes$BOTTLE_ID), decreasing = TRUE)
out_df$NumSibs

# Note that the rownames are the index - 1, not the FK_FISH_ID's
full_run$Partition
bottle_data <- sapply(unique(CMFISH14a.gcl$attributes$BOTTLE_ID), function(bttle) {paste(which(CMFISH14a.gcl$attributes$BOTTLE_ID == bttle) - 1, collapse = "-")})
#sapply(unique(CMFISH14a.gcl$attributes$BOTTLE_ID), function(bttle) {paste(CMFISH14a.gcl$attributes$FK_FISH_ID[which(CMFISH14a.gcl$attributes$BOTTLE_ID == bttle)], collapse = "-")})

# Are all sibships from same bottle?
head(full_run$Partition)

## Note: considering all sibling relationships, regardless of posterior =======
# List by bottle of number of siblings belonging to that bottle index sibgroup
sibgroups_bottle_lst <- lapply(bottle_data, function(bttle) {sapply(strsplit(x = rownames(full_run$Partition), split = "-"), function(index) sum(index %in% unlist(strsplit(x = bttle, split = "-"))))})
# Number of fish per bottle
sapply(sibgroups_bottle_lst, sum)
sum(sapply(sibgroups_bottle_lst, sum)) # 553 alevin, this is the ML full-sib pedigree

# Number of sibling groups in a bottle
sapply(sibgroups_bottle_lst, function(bttle) {sum(bttle > 0)})
plot(table(sapply(sibgroups_bottle_lst, function(bttle) {sum(bttle > 0)})))


## Note: considering only sibling relationships with posterior > 0.8 ==========
# List by bottle of number of siblings belonging to that bottle index sibgroup
sibgroups_bottle_lst_posterior <- lapply(bottle_data, function(bttle) {sapply(strsplit(x = rownames(full_run$Partition), split = "-")[1:min(which(out_df$Posterior < 0.8)) - 1], function(index) sum(index %in% unlist(strsplit(x = bttle, split = "-"))))})
# Number of fish per bottle
sapply(sibgroups_bottle_lst_posterior, sum)
sum(sapply(sibgroups_bottle_lst_posterior, sum)) # 443 alevin, this is the ML full-sib pedigree for relationships with posterior > 0.8

# Number of sibling groups in a bottle
sapply(sibgroups_bottle_lst_posterior, function(bttle) {sum(bttle > 0)})
plot(table(sapply(sibgroups_bottle_lst_posterior, function(bttle) {sum(bttle > 0)})))

# Take-away ===================================================================
# Most bottles only contain 1 or two full-sib groupings, but some contain more

# Are all siblings in the same bottle? ========================================
bottle_sibgroups_lst <- setNames(object = lapply(rownames(full_run$Partition), function(sibgroup) {sapply(strsplit(x = bottle_data, split = "-"), function(bttle) sum(unlist(bttle) %in% unlist(strsplit(x = sibgroup, split = "-"))))}), 
                                 nm =  seq(length(bottle_sibgroups_lst)))

# Number of bottles for a sibling group ALL
sapply(bottle_sibgroups_lst, function(sibgroup) {sum(sibgroup > 0)})
table(sapply(bottle_sibgroups_lst, function(sibgroup) {sum(sibgroup > 0)}))

# Number of bottles for a sibling group with posterior > 0.8
sapply(bottle_sibgroups_lst[1:min(which(out_df$Posterior < 0.8)) - 1], function(sibgroup) {sum(sibgroup > 0)})
table(sapply(bottle_sibgroups_lst[1:min(which(out_df$Posterior < 0.8)) - 1], function(sibgroup) {sum(sibgroup > 0)}))

# Take-away ===================================================================
# Most full sibling groups only occur in 1 bottle, but 8 / 108 with posterior > 0.8 are in 2 bottles
# There are NO full sibling groups spread out over >2 bottles

bottle_sibgroups_lst[which(sapply(bottle_sibgroups_lst, function(sibgroup) {sum(sibgroup > 0)}) > 1)] # most are very close together

# Which sibgroups have posterior > 0.8 and are in more than 1 bottle? =========
# All
sibgroups_2bottle <- which(sapply(bottle_sibgroups_lst, function(sibgroup) {sum(sibgroup > 0)}) > 1) 
# Only with posterior > 0.8
sibgroups_2bottle[sibgroups_2bottle < min(which(out_df$Posterior < 0.8)) - 1]
# Subset bottle_sibgroups_lst
sibgroups_2bottle_lst <- setNames(object = bottle_sibgroups_lst[sibgroups_2bottle[sibgroups_2bottle < min(which(out_df$Posterior < 0.8)) - 1]], nm = sibgroups_2bottle[sibgroups_2bottle < min(which(out_df$Posterior < 0.8)) - 1])
# Just get bottle names
lapply(sibgroups_2bottle_lst, function(sibgroup) {names(sibgroup)[sibgroup > 0]})


# FURTHER WORK
# Check out "V:\DATA\Collection raw data sheets\Chum\SE Alaska\2014\Fish Creek Chum Alevin 2014.xlsx" for GPS coordinates

# All very close to each other, but not necessarily all from the same waypoint (all waypoints < 100' away)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Repeat with QCed markers - post SNP selection ==============================
## Create an input file from 2014 Fish Creek Chum alevin (i.e. dataframe) =====
# nrow = number of alevin genotyped
# ncol = 2 * number of SNPs genotyped
# values are integer values of genotypes (A = 1, C = 2, G = 3, T = 4)

setwd("V:/WORK/Chum/AHRP/Objects")

# Load Fish Creek 2014 alevin genotypes and loci176
CMFISHCR14a.gcl <- dget(file = "CMFISHCR14a.gcl.txt")
loci163 <- dget(file = "loci163.txt")
loci129 <- dget(file = "loci129.txt")
LocusControl <- dget(file = "Phase 1/OriginalLocusControl.txt")

str(CMFISHCR14a.gcl)
str(CMFISHCR14a.gcl$scores)

# Function to create input file for fullsniplings
fullsniplings.GCL <- function(silly, loci){
  
  temp.gcl = get(paste(silly, ".gcl", sep = ""))
  
  # Only include diploid loci in markerset
  loci = loci[which(LocusControl$ploidy[loci] == 2)]
  
  # Subset genotype scores for marker set
  alevin_genos_array = temp.gcl$scores[, loci, ]
  
  # Create a "two column" format dataframe of genotypes
  alevin_genos_mat <- matrix(data = NA, nrow = dim(alevin_genos_array)[1], ncol = dim(alevin_genos_array)[2] * dim(alevin_genos_array)[3], dimnames = list(temp.gcl$attributes$FK_FISH_ID, paste(rep(loci, each = 2), 1:2, sep = "_")))
  for(locus in 1:dim(alevin_genos_array)[2]){
    alevin_genos_mat[, (locus * 2 - 1):(locus * 2)] <- alevin_genos_array[, locus, ]
  }
  
  alevin_genos_mat_num <- gsub(pattern = "A", replacement = 1, x = alevin_genos_mat)
  alevin_genos_mat_num <- gsub(pattern = "C", replacement = 2, x = alevin_genos_mat_num)
  alevin_genos_mat_num <- gsub(pattern = "G", replacement = 3, x = alevin_genos_mat_num)
  alevin_genos_mat_num <- gsub(pattern = "T", replacement = 4, x = alevin_genos_mat_num)
  class(alevin_genos_mat_num) <- "integer"
  
  alevin_genos_df <- as.data.frame(alevin_genos_mat_num)
  
  assign(x = paste(silly, "_fullsniplings_input_loci_", length(loci), sep = ""), value = alevin_genos_df, pos = 1)
}

# All markers that passed gated measures (HWE and LD)
fullsniplings.GCL(silly = "CMFISHCR14a", loci = loci163)

# Subset of those that passed gated measures with (MAF > 0.05)
fullsniplings.GCL(silly = "CMFISHCR14a", loci = loci129)

# Run fullsniplings ===========================================================
full_run_loci160  <- run_mcmc(CMFISHCR14a_fullsniplings_input_loci_160, burn_in = 100, num_sweeps = 500)
dput(x = full_run_loci160, file = "fullsniplings_full_run_loci160.txt")

full_run_loci129  <- run_mcmc(CMFISHCR14a_fullsniplings_input_loci_129, burn_in = 100, num_sweeps = 500)
dput(x = full_run_loci129, file = "fullsniplings_full_run_loci129.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################### STOPPED HERE ############################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("V:/WORK/Chum/AHRP/Objects")
full_run_loci129 <- dget(file = "fullsniplings_full_run_loci129.txt")
full_run_loci160 <- dget(file = "fullsniplings_full_run_loci160.txt")


# Base plot
plot(full_run_loci160$Partition$Visits, type = "p", pch = 16, cex = full_run_loci160$Partition$NumSibs, ylim = c(0, 500))
abline(h = 400, lwd = 5, col = 2)

plot(full_run_loci129$Partition$Visits, type = "p", pch = 16, cex = full_run_loci129$Partition$NumSibs, ylim = c(0, 500))
abline(h = 400, lwd = 5, col = 2)

# Eric Anderson's ggplot from TPB publication
output <- full_run_loci160 # CHANGE THIS

lscale <- 0.02
out_df <- data.frame(Index = seq(dim(output$Partition)[1]), Posterior = output$Partition$Visits / max(output$Partition$Visits), NumSibs = output$Partition$NumSibs)
out_df$ylo <- out_df$Posterior - lscale * out_df$NumSibs
out_df$yhi <- out_df$Posterior + lscale * out_df$NumSibs

require(ggplot2)
require(plyr)

bugle_theme <- function (base_size = 12, base_family = "")
{
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(axis.text = element_text(size = rel(0.8)),
          axis.ticks = element_line(colour = "black"),
          legend.key = element_rect(colour = "grey50"),
          legend.position = "bottom",
          panel.background = element_rect(fill = "white",
                                          colour = NA),
          panel.border = element_rect(fill = NA,
                                      colour = "grey50"),
          panel.grid.major = element_line(colour = "grey98",
                                          size = 0.2),
          panel.grid.minor = element_line(colour = "grey98",
                                          size = 0.5),
          strip.background = element_rect(fill = "grey80",
                                          colour = "grey50", size = 0.2)
    )
}

ggplot(data = out_df, aes(x = Index, y = Posterior)) +
  geom_linerange(aes(ymin = ylo,
                     ymax = yhi),
                 size = 0.1) +
  geom_line(size = 0.15) +
  bugle_theme()

# How does this line up to sampling? Are all alevin in a bottle siblings? =====
table(CMFISHCR14a.gcl$attributes$BOTTLE_ID)
sort(table(CMFISHCR14a.gcl$attributes$BOTTLE_ID), decreasing = TRUE)
out_df$NumSibs

# Note that the rownames are the index - 1, not the FK_FISH_ID's
full_run_loci129$Partition
bottle_data <- sapply(unique(CMFISHCR14a.gcl$attributes$BOTTLE_ID), function(bttle) {paste(which(CMFISHCR14a.gcl$attributes$BOTTLE_ID == bttle) - 1, collapse = "-")})
#sapply(unique(CMFISHCR14a.gcl$attributes$BOTTLE_ID), function(bttle) {paste(CMFISHCR14a.gcl$attributes$FK_FISH_ID[which(CMFISHCR14a.gcl$attributes$BOTTLE_ID == bttle)], collapse = "-")})

# Are all sibships from same bottle?
head(full_run_loci129$Partition)

## Note: considering all sibling relationships, regardless of posterior =======
# List by bottle of number of siblings belonging to that bottle index sibgroup
sibgroups_bottle_lst <- lapply(bottle_data, function(bttle) {sapply(strsplit(x = rownames(full_run_loci129$Partition), split = "-"), function(index) sum(index %in% unlist(strsplit(x = bttle, split = "-"))))})
# Number of fish per bottle
sapply(sibgroups_bottle_lst, sum)
sum(sapply(sibgroups_bottle_lst, sum)) # 553 alevin, this is the ML full-sib pedigree

# Number of sibling groups in a bottle
sapply(sibgroups_bottle_lst, function(bttle) {sum(bttle > 0)})
plot(table(sapply(sibgroups_bottle_lst, function(bttle) {sum(bttle > 0)})))


## Note: considering only sibling relationships with posterior > 0.8 ==========
# List by bottle of number of siblings belonging to that bottle index sibgroup
sibgroups_bottle_lst_posterior <- lapply(bottle_data, function(bttle) {sapply(strsplit(x = rownames(full_run_loci129$Partition), split = "-")[1:min(which(out_df$Posterior < 0.8)) - 1], function(index) sum(index %in% unlist(strsplit(x = bttle, split = "-"))))})
# Number of fish per bottle
sapply(sibgroups_bottle_lst_posterior, sum)
sum(sapply(sibgroups_bottle_lst_posterior, sum)) # 443 alevin, this is the ML full-sib pedigree for relationships with posterior > 0.8

# Number of sibling groups in a bottle
sapply(sibgroups_bottle_lst_posterior, function(bttle) {sum(bttle > 0)})
plot(table(sapply(sibgroups_bottle_lst_posterior, function(bttle) {sum(bttle > 0)})))

# Take-away ===================================================================
# Most bottles only contain 1 or two full-sib groupings, but some contain more

# Are all siblings in the same bottle? ========================================
bottle_sibgroups_lst <- setNames(object = lapply(rownames(full_run_loci129$Partition), function(sibgroup) {sapply(strsplit(x = bottle_data, split = "-"), function(bttle) sum(unlist(bttle) %in% unlist(strsplit(x = sibgroup, split = "-"))))}), 
                                 nm =  seq(dim(full_run_loci129$Partition)[1]))

# Number of bottles for a sibling group ALL
sapply(bottle_sibgroups_lst, function(sibgroup) {sum(sibgroup > 0)})
table(sapply(bottle_sibgroups_lst, function(sibgroup) {sum(sibgroup > 0)}))

# Number of bottles for a sibling group with posterior > 0.8
bottle_sibgroups_lst_posterior <- bottle_sibgroups_lst[1:min(which(out_df$Posterior < 0.8)) - 1]

sapply(bottle_sibgroups_lst_posterior, function(sibgroup) {sum(sibgroup > 0)})
table(sapply(bottle_sibgroups_lst_posterior, function(sibgroup) {sum(sibgroup > 0)}))

# Take-away ===================================================================
# Most full sibling groups only occur in 1 bottle, but 7 / 105 with posterior > 0.8 are in 2 bottles
# There are NO full sibling groups spread out over >2 bottles

bottle_sibgroups_lst[which(sapply(bottle_sibgroups_lst, function(sibgroup) {sum(sibgroup > 0)}) > 1)] # most are very close together

# Which sibgroups have posterior > 0.8 and are in more than 1 bottle? =========
# All
sibgroups_2bottle <- which(sapply(bottle_sibgroups_lst, function(sibgroup) {sum(sibgroup > 0)}) > 1) 
# Only with posterior > 0.8
sibgroups_2bottle[sibgroups_2bottle < min(which(out_df$Posterior < 0.8)) - 1]
# Subset bottle_sibgroups_lst
sibgroups_2bottle_lst <- setNames(object = bottle_sibgroups_lst[sibgroups_2bottle[sibgroups_2bottle < min(which(out_df$Posterior < 0.8)) - 1]], nm = sibgroups_2bottle[sibgroups_2bottle < min(which(out_df$Posterior < 0.8)) - 1])
# Just get bottle names
lapply(sibgroups_2bottle_lst, function(sibgroup) {names(sibgroup)[sibgroup > 0]})


# Compare loci 160 to loci 129 ================================================
bottle_sibgroups_lst_loci160 <- bottle_sibgroups_lst
bottle_sibgroups_lst_posterior_loci160 <- bottle_sibgroups_lst_posterior
sibgroups_bottle_lst_loci160 <- sibgroups_bottle_lst
sibgroups_bottle_lst_posterior_loci160 <- sibgroups_bottle_lst_posterior
sibgroups_2bottle_lst_loci160 <- sibgroups_2bottle_lst
lapply(sibgroups_2bottle_lst_loci160, function(sibgroup) {names(sibgroup)[sibgroup > 0]})

bottle_sibgroups_lst_loci129 <- bottle_sibgroups_lst
bottle_sibgroups_lst_posterior_loci129 <- bottle_sibgroups_lst_posterior
sibgroups_bottle_lst_loci129 <- sibgroups_bottle_lst
sibgroups_bottle_lst_posterior_loci129 <- sibgroups_bottle_lst_posterior
sibgroups_2bottle_lst_loci129 <- sibgroups_2bottle_lst
lapply(sibgroups_2bottle_lst_loci129, function(sibgroup) {names(sibgroup)[sibgroup > 0]})

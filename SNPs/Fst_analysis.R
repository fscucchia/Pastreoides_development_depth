
###### Estimation of Fst statistic from transcriptome-derived SNPs ######

## Load libraries

library(gdsfmt)
library(SNPRelate)
library(hierfstat)
library(tidyr)
library(dplyr)

## Load the pruned SNPs

bed.fn <- "/data/home/mass/fscucchia/programs/plink_2/P_astreoides_vcf_pruned.bed"
fam.fn <- "/data/home/mass/fscucchia/programs/plink_2/P_astreoides_vcf_pruned.fam"
bim.fn <- "/data/home/mass/fscucchia/programs/plink_2/P_astreoides_vcf_pruned.bim"

## Create gds file
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "plink_pruned.gds")
snpgdsSummary("plink_pruned.gds")
# [1]
# The file name: /data/home/mass/fscucchia/Bermuda/output/SNPRelate/plink_pruned.gds
# The total number of samples: 12
# The total number of SNPs: 9060
# SNP genotypes are stored in SNP-major mode (Sample X SNP).

## Open the gds file
genofile <- snpgdsOpen("plink_pruned.gds")

## Get sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

## Get population information
pop_code <- scan("/data/home/mass/fscucchia/Bermuda/output/SNPRelate/pop.txt", what=character())
#see the file here https://github.com/fscucchia/Pastreoides_development_depth/blob/main/SNPs/pop.txt

head(cbind(sample.id, pop_code)) #assumes the order of sample IDs is as the same as population codes

## Make a PCA to look at the genetic clusters
pca <- snpgdsPCA(genofile, autosome.only=FALSE)
#[1] Working space: 12 samples, 9,060 SNPs

# Make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
    pop = factor(pop_code)[match(pca$sample.id, sample.id)],
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)

# Draw
pdf(paste0('PCA_prunedSNPs','.pdf'))
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop), xlab="eigenvector 2", ylab="eigenvector 1")
legend("bottomright", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))
dev.off()

# Look at the % variance explained by the components
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
#[1] 36.23 17.91 14.22 10.71  6.25  5.38


########## Fst  Estimation

x2<-2-snpgdsGetGeno(genofile) #snpgdsVCF2GDS stores the number of reference alleles, we want the the number of alternate alleles
#[1] Genotype matrix: 12 samples X 9060 SNPs

#add column with samples names 
vec <- c("D20","D21","D22","L7","L8","L9",
         "R13B","R15B","R33","S4","S5","S7") 

locus_transpose1 <- cbind(x2, sample = vec)  
locus_transpose2=as.data.frame(locus_transpose1)

#move last column with samples names to first
locus_transpose3 <- locus_transpose2 %>%
  select(sample, everything())

#add column with with groups as numbers, adult_meso is 1, planu_meso is 2 and so on
vec2 <- c("1","1","1","2","2","2","3","3","3",
          "4","4","4")

#add column with depths
vec3 <- c("45","45","45","45","45","45",
         "10","10","10","10","10","10") 

locus_transpose4 <- cbind(locus_transpose3, group = vec2) 
locus_transpose5 <- cbind(locus_transpose4, depth = vec3) 

#move last columns with depth and colony to 2nd and 3rd places

moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", 
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}

#move group column after sample column using the moveme function above
locus_transpose6=locus_transpose5[moveme(names(locus_transpose5), "group after sample")]

#move depth column after group column using the moveme function above
locus_transpose7=locus_transpose6[moveme(names(locus_transpose6), "depth after group")]

locus_transpose7$sample <- NULL #don't need this column

locus_transpose8 <- lapply(locus_transpose7,as.numeric) #pairwise.WCfst only works with numerics
locus_transpose9=as.data.frame(locus_transpose8)

Fst_between_groups = pairwise.WCfst(locus_transpose9[,-2],diploid=TRUE) #Fst between groups of samples

## Permute Fst to get p values
#see the original script here https://lists.r-forge.r-project.org/pipermail/adegenet-forum/2011-February/000214.html

library(adegenet)

NBPERM <- 999 # this is the number of permutations used for the p-values
mat.perm <- lapply(1:NBPERM, function(i) pairwise.WCfst(locus_transpose9[,-2],diploid=TRUE))

#Fst_between_groups contains original Fst values, mat.perm is a list with NPERM matrices of permuted Fst values. 
#To get e.g. right-tail p-values, you can just count the proportion of mat.obs >= mat.perm; e.g. for the first pair of populations:
mean(c(Fst_between_groups[1,2] < na.omit(sapply(1:NBPERM, function(i) mat.perm[[i]][1,2])), TRUE))

#In the above command, "na.omit(sapply(1:NBPERM, function(i) mat.perm[[i]][1,2]))" is a vector of permuted values for this pair of populations across 
#all replicates; c(..., TRUE) is added because the observed value is always added to the permuted values (it is one of the possible permutations of groups).
#In practice, it is easier to convert the results as objects of the class randtest (class for Monte Carlo test in ade4):
library(ade4)
Randtest <- as.randtest(na.omit(sapply(1:NBPERM, function(i) mat.perm[[i]][1,2])), Fst_between_groups[1,2], alter="greater")

#To have it done for all pairs of populations:
allTests <- list()
for(i in 1:(nrow(Fst_between_groups)-1)){
  for(j in 2:nrow(Fst_between_groups)){
    allTests[[paste(rownames(Fst_between_groups)[i],rownames(Fst_between_groups)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) mat.perm[[k]][i,j])), Fst_between_groups[i,j], alter="greater")
  }
}  

allTests
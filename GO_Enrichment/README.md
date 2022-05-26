## GO enrichment analysis for both the DE and WGCNA data of _P. astreoides_ planulae and adult samples from shallow and mesophotic reefs ######

### GO Enrichment - Differential Expression Data 
```
#This script is based on the work of Erin Chille (https://github.com/echille/Mcapitata_OA_Developmental_Gene_Expression_Timeseries/tree/main/5-Planula-GO-Enrichment-Analysis) with some modifications and adjustments.

### Set up workspace

setwd("C:/Federica/Clean_workspace/Bermuda_scripts/Enrichment")

#Load libraries

library(goseq)
library(tidyverse)
library(GSEABase)               #BiocManager::install("GSEABase")
library(data.table)
library(ggplot2)
library(cowplot)                #install.packages("cowplot")
library(patchwork)
library(dplyr)
library(tidyr)

#treatment information
treatmentinfo <- read.csv("RNAseq_data.csv", header = TRUE, sep = ";")
str(treatmentinfo)
head(treatmentinfo)

#gene count matrix
gcount <- as.data.frame(read.csv("gene_count_matrix.csv", row.names="gene_id"))
gcount$gene_id <- rownames(gcount)

## DESeq2 result - ADULTS
DEG.res <- read.csv("DEGs_adult_cluster.csv")[,-1] #this is the output of the DE script https://github.com/fscucchia/Pastreoides_development_depth/tree/main/DE

#Load Past annotations
annot_final <- read.csv("Past_annot.csv", header = TRUE, sep = ",")[,-1]

Go.ref <- subset(annot_final, select= c(SeqName, Length)) #Select only relevant information

names(Go.ref)[1] <- "gene_id" #rename column
names(Go.ref)[2] <- "length"

#Filter gcount by available annotations
Go.ref_gcount_merged <- merge(gcount, Go.ref, by = "gene_id")

#Make a dataframe containing the gene_ids and cluster for each cluster.
#Select only gene_id and cluster from DEseq2 res
DEGclust_adults <- subset(DEG.res, select=c(gene_id, cluster))
DEGclust_adults <- unique(DEGclust_adults)
clust1_adults <- filter(DEGclust_adults, cluster=="1")
nrow(clust1_adults) #nrow clust1
clust2_adults <- filter(DEGclust_adults, cluster=="2")
nrow(clust2_adults) #nrow clust2

#Set ID and gene length vectors, and make a binary matrix indicating which genes are differentially expressed. These are used as input to nullp, which for calculates a Probability Weighting Function for each set of DEGs.
#Make ID and length vectors
Go.ref_gcount_merged <- unique(Go.ref_gcount_merged)
dim(Go.ref_gcount_merged)
IDvector <- Go.ref_gcount_merged$gene_id
lengthVector <- Go.ref_gcount_merged$length

#Cluster 1
clust1.genes <- as.vector(clust1_adults$gene_id)
clust1.genes=as.integer(Go.ref_gcount_merged$gene_id%in%clust1.genes)
names(clust1.genes)=Go.ref_gcount_merged$gene_id
length(clust1.genes)
length(names(clust1.genes))
length(unique(names(clust1.genes)))

#Cluster 2
clust2.genes <- as.vector(clust2_adults$gene_id)
clust2.genes=as.integer(Go.ref_gcount_merged$gene_id%in%clust2.genes)
names(clust2.genes)=Go.ref_gcount_merged$gene_id
length(clust2.genes)

pwf.C1<-nullp(DEgenes=clust1.genes, id=IDvector, bias.data=lengthVector) #weight vector by length of gene
pwf.C2<-nullp(clust2.genes, IDvector, bias.data=lengthVector) #weight vector by length of gene

## Prepare GO term dataframe

GO.annot <- subset(annot_final, select= c(SeqName, GO_IDs)) #Select only relevant information
names(GO.annot)[1] <- "gene_id" #rename column

GO.annot.na <- filter(GO.annot, GO_IDs!="NA;NA") #Remove NAs
GO.annot.na_cleaned <- GO.annot.na[!grepl(";;", GO.annot.na$GO_IDs), ] #Remove ";;", that is genes with no GO term

splitted <- strsplit(as.character(GO.annot.na_cleaned$GO_IDs), ";") #split into multiple GO ids
GO.terms <- data.frame(v1 = rep.int(GO.annot.na_cleaned$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
colnames(GO.terms) <- c("gene_id", "GO.ID")

GO.terms$GO.ID <- as.factor(GO.terms$GO.ID)
GO.terms$gene_id <- as.factor(GO.terms$gene_id)
GO.terms$GO.ID <- gsub(" ", "", GO.terms$GO.ID)
GO.terms <- unique(GO.terms)

#Find enriched GO terms, "selection-unbiased testing for category enrichment amongst significantly expressed genes for RNA-seq data"
GOwall.C1_adults <- goseq(pwf.C1, GO.ref$gene_id, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
GOwall.C2_adults <- goseq(pwf.C2, GO.ref$gene_id, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)

#Find only enriched GO terms that are statistically significant at cutoff
C1.GO.sigp_adults<-GOwall.C1_adults$category[GOwall.C1_adults$over_represented_pvalue<.05]
C1.GO.sigp_adults<-data.frame(C1.GO.sigp_adults)
colnames(C1.GO.sigp_adults) <- c("category")
C1.GO.sigp_adults <- merge(C1.GO.sigp_adults, GOwall.C1_adults, by="category")
C1.GO.sigp_adults <- C1.GO.sigp_adults[order(C1.GO.sigp_adults$ontology, C1.GO.sigp_adults$over_represented_pvalue, -C1.GO.sigp_adults$numDEInCat),]
C1.GO.sigp_adults$term <- as.factor(C1.GO.sigp_adults$term)

C2.GO.sigp_adults<-GOwall.C2_adults$category[GOwall.C2_adults$over_represented_pvalue<.05]
C2.GO.sigp_adults<-data.frame(C2.GO.sigp_adults)
colnames(C2.GO.sigp_adults) <- c("category")
C2.GO.sigp_adults <- merge(C2.GO.sigp_adults, GOwall.C2_adults, by="category")
C2.GO.sigp_adults <- C2.GO.sigp_adults[order(C2.GO.sigp_adults$ontology, C2.GO.sigp_adults$over_represented_pvalue, -C2.GO.sigp_adults$numDEInCat),]
C2.GO.sigp_adults$term <- as.factor(C2.GO.sigp_adults$term)

###Save significant terms
write.csv(C1.GO.sigp_adults, "C1.GO.sigp_adults.csv", row.names = FALSE)
write.csv(C2.GO.sigp_adults, file = "C2.GO.sigp_adults.csv", row.names = FALSE)


## DESeq2 result - PLANULAE
DEG.res <- read.csv("DEGs_planulae_cluster.csv", header = TRUE, sep = "," )[,-1]

#filter DEGs for log2FoldChange>|2|
DEG.res <- filter(DEG.res, log2FoldChange > 2 | log2FoldChange < (-2))
nrow(DEG.res) #from 2898 to 561

#Make a dataframe containing the gene_ids and cluster for each cluster.
#Select only gene_id and cluster from DEseq2 res
DEGclust_planulae <- subset(DEG.res, select=c(gene_id, cluster))
DEGclust_planulae$cluster <- as.factor(DEGclust_planulae$cluster)
DEGclust_planulae <- unique(DEGclust_planulae)
clust1_planulae <- filter(DEGclust_planulae, cluster=="1")
clust2_planulae <- filter(DEGclust_planulae, cluster=="2")

#Cluster 1
clust1.genes <- as.vector(clust1_planulae$gene_id)
clust1.genes=as.integer(Go.ref_gcount_merged$gene_id%in%clust1.genes)
names(clust1.genes)=Go.ref_gcount_merged$gene_id
length(clust1.genes)
length(names(clust1.genes))
length(unique(names(clust1.genes)))

#Cluster 2
clust2.genes <- as.vector(clust2_planulae$gene_id)
clust2.genes=as.integer(Go.ref_gcount_merged$gene_id%in%clust2.genes)
names(clust2.genes)=Go.ref_gcount_merged$gene_id
length(clust2.genes)

pwf.C1<-nullp(DEgenes=clust1.genes, id=IDvector, bias.data=lengthVector) #weight vector by length of gene
pwf.C2<-nullp(clust2.genes, IDvector, bias.data=lengthVector) #weight vector by length of gene

#Find enriched GO terms, "selection-unbiased testing for category enrichment amongst significantly expressed genes for RNA-seq data"
GOwall.C1_planulae <- goseq(pwf.C1, GO.ref$gene_id, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
GOwall.C2_planulae <- goseq(pwf.C2, GO.ref$gene_id, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)

#Find only enriched GO terms that are statistically significant at cutoff
C1.GO.sigp_planulae<-GOwall.C1_planulae$category[GOwall.C1_planulae$over_represented_pvalue<.05]
C1.GO.sigp_planulae<-data.frame(C1.GO.sigp_planulae)
colnames(C1.GO.sigp_planulae) <- c("category")
C1.GO.sigp_planulae <- merge(C1.GO.sigp_planulae, GOwall.C1_planulae, by="category")
C1.GO.sigp_planulae <- C1.GO.sigp_planulae[order(C1.GO.sigp_planulae$ontology, C1.GO.sigp_planulae$over_represented_pvalue, -C1.GO.sigp_planulae$numDEInCat),]
C1.GO.sigp_planulae$term <- as.factor(C1.GO.sigp_planulae$term)

C2.GO.sigp_planulae<-GOwall.C2_planulae$category[GOwall.C2_planulae$over_represented_pvalue<.05]
C2.GO.sigp_planulae<-data.frame(C2.GO.sigp_planulae)
colnames(C2.GO.sigp_planulae) <- c("category")
C2.GO.sigp_planulae <- merge(C2.GO.sigp_planulae, GOwall.C2_planulae, by="category")
C2.GO.sigp_planulae <- C2.GO.sigp_planulae[order(C2.GO.sigp_planulae$ontology, C2.GO.sigp_planulae$over_represented_pvalue, -C2.GO.sigp_planulae$numDEInCat),]
C2.GO.sigp_planulae$term <- as.factor(C2.GO.sigp_planulae$term)

###Save significant terms
write.csv(C1.GO.sigp_planulae, "C1.GO.sigp_planulae.csv", row.names = FALSE)
write.csv(C2.GO.sigp_planulae, file = "C2.GO.sigp_planulae.csv", row.names = FALSE)
```

### GO Enrichment - WGCNA Data 
```






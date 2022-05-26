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
setwd("/data/home/mass/fscucchia/Bermuda/output/WGCNA")
load(".RData") #load all WCGNA data

#Prepare dataframe
genes.GO <- as.data.frame(t(datExpr))
genes.GO <- cbind(gene_id = rownames(genes.GO), genes.GO)
rownames(genes.GO) <- NULL

head(genes.GO)
# [1]       gene_id      D20      D21      D22        L7        L8        L9
# 1 Pastreoides01556 4.372019 4.372019 4.372019  4.372019  6.353233  6.146853
# 2 Pastreoides37429 4.372019 6.525170 5.977893  7.810618  8.158924  7.427095
# 3 Pastreoides47565 4.372019 6.897309 7.535716  9.190954  9.277749  8.973132
# 4 Pastreoides47564 4.372019 8.246516 8.384271  9.567475  9.538375  8.995897
# 5 Pastreoides48395 4.372019 5.958641 5.727800  4.990296  4.974196  5.337775
# 6 Pastreoides47566 7.903908 8.644152 8.889982 11.205382 11.347467 11.173282
#       R13B     R15B      R33        S4        S5        S7
# 1 4.372019 4.372019 4.372019  5.385981  5.232322  5.464851
# 2 6.324533 7.206471 7.095487  7.300074  7.339899  7.140677
# 3 6.857971 6.872340 7.886073  9.316904  9.383418  9.299059
# 4 8.168371 8.157733 8.090708  9.508662  9.662097  9.380658
# 5 4.372019 6.427276 6.887108  6.964790  6.639221  7.210621
# 6 8.938027 8.797427 8.465760 11.348584 11.344754 11.234663

### Select interesting Clusters (the ones with a pattern by depth)
#here I looked at the difference in mean module eigengene value (boxplot profiles) and decided to analyze cluster 5 and 7, that are the ones that 
#show the greater change in mean value in the meso vs shal adult comparison and in the meso vs shal planulae comparison.

## Cluster 5 modules
genes_clust5.GO <- genes.GO[,c(1:13)]
geneColor <- geneInfo %>% select(gene_id, moduleColor) #Make a dataframe containing just gene_id and moduleColor

module_cluster5 <- geneColor %>% filter(moduleColor == "thistle3" | moduleColor == "bisque4"| moduleColor == "coral3" | moduleColor == "coral1" | moduleColor == "turquoise")  #Make a dataframe containing just gene_id and moduleColor of cluster1

geneColor$gene_id <- as.factor(geneColor$gene_id) #Make factor for merge
genes_clust5.GO$gene_id <- as.factor(genes_clust5.GO$gene_id) #Make factor for merge
genes_clust5.GO <- merge(module_cluster5, genes_clust5.GO)


## Cluster 7 modules
genes_clust7.GO <- genes.GO[,c(1:13)]
geneColor <- geneInfo %>% select(gene_id, moduleColor) #Make a dataframe containing just gene_id and moduleColor

module_cluster7 <- geneColor %>% filter(moduleColor == "indianred3" | moduleColor == "brown2"| moduleColor == "paleturquoise")  #Make a dataframe containing just gene_id and moduleColor of cluster1

geneColor$gene_id <- as.factor(geneColor$gene_id) #Make factor for merge
genes_clust7.GO$gene_id <- as.factor(genes_clust7.GO$gene_id) #Make factor for merge
genes_clust7.GO <- merge(module_cluster7, genes_clust7.GO)

#Build a data frame that links the gene IDs, modules, and counts of expressed genes and the gene lengths.

GOref <- merge(genes.GO, GO.annot, by.x="gene_id")
GOref <- merge(genes.GO, GO.annot, by ="gene_id")
head(GOref)
dim(GOref)

##GOseq requires a vector of all genes and all differentially expressed genes. 

####### Cluster 5 ######
cluster5_gene.vector <- as.vector(genes_clust5.GO$gene_id)
cluster5_gene.vector=as.integer(GOref$gene_id %in% cluster5_gene.vector)
names(cluster5_gene.vector)=GOref$gene_id
head(cluster5_gene.vector)

cluster5_ID.vector <- as.character(GOref$gene_id) #Make ID vector
head(cluster5_ID.vector)
dim(cluster5_ID.vector)
cluster5_Length.vector <- GOref$length #Make length vector
head(cluster5_Length.vector)

#Calculate Probability Weighting Function
pwf.cluster5 <- nullp(cluster5_gene.vector, cluster5_ID.vector, bias.data=cluster5_Length.vector) #weight vector by length of gene

#Prepare GO term dataframe
GO.annot <- select(geneInfo, gene_id, Annotation.GO.ID)
splitted <- strsplit(as.character(GO.annot$Annotation.GO.ID), ";") #split into multiple GO ids
GO.terms <- data.frame(v1 = rep.int(GO.annot$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
colnames(GO.terms) <- c("gene_id", "GO.ID")
head(GO.terms)
tail(GO.terms)

GO.terms$GO.ID<- as.character(GO.terms$GO.ID)
GO.terms$GO.ID <- replace_na(GO.terms$GO.ID, "unknown")
GO.terms$GO.ID <- as.factor(GO.terms$GO.ID)
GO.terms$gene_id <- as.factor(GO.terms$gene_id)
GO.terms$GO.ID <- gsub(" ", "", GO.terms$GO.ID)
dim(GO.terms)

##Find enriched GO terms, "selection-unbiased testing for category enrichment amongst significantly expressed genes for RNA-seq data"
GOwall.cluster5<- goseq(pwf.cluster5, GOref$gene_id, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
dim(GOwall.cluster5)
head(GOwall.cluster5)

#Find only enriched GO terms that are statistically significant at cutoff
cluster5.GO.05<-GOwall.cluster5$category[GOwall.cluster5$over_represented_pvalue<.05]
cluster5.GO.05<-data.frame(cluster5.GO.05)
colnames(cluster5.GO.05) <- c("category")
cluster5.GO.05 <- merge(cluster5.GO.05, GOwall.cluster5, by="category")
cluster5.GO.05 <- cluster5.GO.05[order(cluster5.GO.05$ontology, cluster5.GO.05$over_represented_pvalue, -cluster5.GO.05$numDEInCat),]
cluster5.GO.05$term <- as.factor(cluster5.GO.05$term)
dim(cluster5.GO.05) #Number of sig GO terms

#Save significant terms
write.csv(cluster5.GO.05, file = "GO.05.cluster5.csv", row.names = FALSE)

####### Cluster 7 ######
cluster7_gene.vector <- as.vector(genes_clust7.GO$gene_id)
cluster7_gene.vector=as.integer(GOref$gene_id %in% cluster7_gene.vector)
names(cluster7_gene.vector)=GOref$gene_id
head(cluster7_gene.vector)

cluster7_ID.vector <- as.character(GOref$gene_id) #Make ID vector
head(cluster7_ID.vector)
dim(cluster7_ID.vector)
cluster7_Length.vector <- GOref$length #Make length vector
head(cluster7_Length.vector)

#Calculate Probability Weighting Function
pwf.cluster7 <- nullp(cluster7_gene.vector, cluster7_ID.vector, bias.data=cluster7_Length.vector) #weight vector by length of gene

#Prepare GO term dataframe
GO.annot <- select(geneInfo, gene_id, Annotation.GO.ID)
splitted <- strsplit(as.character(GO.annot$Annotation.GO.ID), ";") #split into multiple GO ids
GO.terms <- data.frame(v1 = rep.int(GO.annot$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
colnames(GO.terms) <- c("gene_id", "GO.ID")
head(GO.terms)
tail(GO.terms)

GO.terms$GO.ID<- as.character(GO.terms$GO.ID)
GO.terms$GO.ID <- replace_na(GO.terms$GO.ID, "unknown")
GO.terms$GO.ID <- as.factor(GO.terms$GO.ID)
GO.terms$gene_id <- as.factor(GO.terms$gene_id)
GO.terms$GO.ID <- gsub(" ", "", GO.terms$GO.ID)
dim(GO.terms)

##Find enriched GO terms, "selection-unbiased testing for category enrichment amongst significantly expressed genes for RNA-seq data"
GOwall.cluster7<- goseq(pwf.cluster7, GOref$gene_id, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
dim(GOwall.cluster7)
head(GOwall.cluster7)

#Find only enriched GO terms that are statistically significant at cutoff
cluster7.GO.05<-GOwall.cluster7$category[GOwall.cluster7$over_represented_pvalue<.05]
cluster7.GO.05<-data.frame(cluster7.GO.05)
colnames(cluster7.GO.05) <- c("category")
cluster7.GO.05 <- merge(cluster7.GO.05, GOwall.cluster7, by="category")
cluster7.GO.05 <- cluster7.GO.05[order(cluster7.GO.05$ontology, cluster7.GO.05$over_represented_pvalue, -cluster7.GO.05$numDEInCat),]
cluster7.GO.05$term <- as.factor(cluster7.GO.05$term)
dim(cluster7.GO.05) #Number of sig GO terms

#Save significant terms
write.csv(cluster7.GO.05, file = "GO.05.cluster7.csv", row.names = FALSE)
```

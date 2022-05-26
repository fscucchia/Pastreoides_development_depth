## RNAseq Weighted Correlation Network Aanalysis of _P. astreoides_ planulae and adult samples from shallow and mesophotic reefs ######

This script is based on the work of [Erin Chille](https://github.com/echille/Mcapitata_Developmental_Gene_Expression_Timeseries/tree/master/2a-WGCNA), with some modifications and adjustments.

### Set up workspace in R 

```
#Load libraries

library("WGCNA")              #BiocManager::install("WGCNA")
library("flashClust")         #install.packages("flashClust")
library("pheatmap")  
library("clusterProfiler")    #BiocManager::install("clusterProfiler")
library("simplifyEnrichment") #BiocManager::install("simplifyEnrichment")
library("genefilter")           #BiocManager::install("genefilter") 
library("DESeq2")               #BiocManager::install("DESeq2")
library("factoextra")           #install.packages("factoextra")
library("NbClust")              #install.packages("NbClust")
library("ComplexHeatmap")       #BiocManager::install("ComplexHeatmap")
library("tidyverse")            
library("RColorBrewer")
library("ggplot2")              
library("goseq")                #BiocManager::install("goseq")
library("gridExtra")            #install.packages("gridExtra")
library("VennDiagram")          #install.packages("VennDiagram")
library("patchwork")            #install.packages("patchwork")
library("dplyr")

#treatment information
treatmentinfo <- read.csv("RNAseq_data.csv", header = TRUE, sep = ";")

#gene count matrix
gcount <- as.data.frame(read.csv("gene_count_matrix.csv", row.names="gene_id"), colClasses = double)
```

### Quality-filter gene counts
```
#Set filter values for PoverA: smallest sample size per treat is 3, so 3/12 (12 samples) is 0.25
#This means that 3 out of 12 (0.25) samples need to have counts over 10.
#So P=25 percent of the samples have counts over A=10. 
filt <- filterfun(pOverA(0.25,10))

#create filter for the counts data
gfilt <- genefilter(gcount, filt)

#identify genes to keep by count filter
gkeep <- gcount[gfilt,]

#identify gene lists
gn.keep <- rownames(gkeep)

#gene count data filtered in PoverA, P percent of the samples have counts over A
gcount_filt <- as.data.frame(gcount[which(rownames(gcount) %in% gn.keep),])

#How many rows do we have before and after filtering?
nrow(gcount) #Before                
# [1] 64636

nrow(gcount_filt) #After
# [1] 35776
```

### Read normalization
```
# Normalize our read counts using VST-normalization in DESeq2
# Construct the DESeq2 dataset

#Merge the age and depth columns into a new column, group. Set group as a factor.
treatmentinfo$group <- factor(treatmentinfo$group, levels = c("adult_meso","adult_shal","planu_meso", "planu_shal"))

treatmentinfo$depth <- factor(treatmentinfo$depth)

#Create a DESeqDataSet design from gene count matrix and labels. Here we set the design to look at 
#any differences in gene expression across samples attributed to depth.

#Set DESeq2 design
gdds <- DESeqDataSetFromMatrix(countData = gcount_filt,
                                  colData = treatmentinfo,
                                  design = ~depth)

# Log-transform the count data using a variance stabilizing transforamtion (vst). 

SF.gdds <- estimateSizeFactors( gdds ) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than for to use vst
print(sizeFactors(SF.gdds)) #View size factors

gvst <- vst(gdds, blind=FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
```

### Compile WGCNA Dataset
```
# Transpose the filtered gene count matrix so that the gene IDs are rows and the sample IDs are columns.

datExpr <- as.data.frame(t(assay(gvst))) #transpose to output to a new data frame with the column names as row names. And make all data numeric

#Check for genes and samples with too many missing values with goodSamplesGenes. There shouldn't be any because we performed pre-filtering
gsg = goodSamplesGenes(datExpr, verbose = 3)
# [1] allOK is TRUE

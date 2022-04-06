## RNAseq Differential Expression Analysis of _P. astreoides_ planulae and adult samples from shallow and mesophotic reefs ######

### Set up workspace in R Studio

```{r}
# Load libraries
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

#load treatment information
treatmentinfo <- read.csv("RNAseq_data2.csv", header = TRUE, sep = ";")
str(treatmentinfo)
head(treatmentinfo)

#load gene count matrix
gcount <- as.data.frame(read.csv("gene_count_matrix.csv", row.names="gene_id")) #the gene_count_matrix.csv is the output of Stringtie
head(gcount)
```
### Construct DESeq2 dataset
#### Pre-filter gene counts
Pre-filtering our dataset to reduce the memory size dataframe, increase the speed of the transformation and improve sensitivity of statistical analysis by removing low-coverage counts. 

```
#Set filter values for PoverA: smallest sample size per treat. is 3, so 3/12 (12 samples) is 0.25. 
#This means that 3 out of 12 (0.25) samples need to have counts over 10. So P=25 percent of the samples have counts over A=10. 
filt <- filterfun(pOverA(0.25,10))

#create filter for the counts data
gfilt <- genefilter(gcount, filt)

#identify genes to keep by count filter
gkeep <- gcount[gfilt,]

#identify gene lists
gn.keep <- rownames(gkeep)

#gene count data filtered in PoverA, P percent of the samples have counts over A
gcount_filt <- as.data.frame(gcount[which(rownames(gcount) %in% gn.keep),])

#Merge the age and depth columns into a new column, group. Set group as a factor.
treatmentinfo$group <- factor(treatmentinfo$group, levels = c("adult_meso","adult_shal","planu_meso", "planu_shal"))

#Create a DESeqDataSet design from gene count matrix and labels. Here we set the design to look at the interaction of age and depth to test for any differences in gene #expression across samples attributed to depth and developmental stage.

#Set DESeq2 design
gdds <- DESeqDataSetFromMatrix(countData = gcount_filt,
                                  colData = treatmentinfo,
                                  design = ~group)
```

#### Visualize gene count data
```{r}
## Log-transform the count data
#Log-transform the data using a variance stabilizing transforamtion (vst) for visualization purposes. This transformation deals with the sampling variability of low counts by #calculating within-group variability.  

#To use vst we first need to calculate the size factors of the samples, that is an estimate of how many reads each sample has compared to the others. 
SF.gdds <- estimateSizeFactors(gdds) #size factors should be less than for to use vst
print(sizeFactors(SF.gdds)) #View size factors. In this case size factors are all less than 4, so vst can be used.

gvst <- vst(gdds, blind=FALSE) 

### Principal component plot of samples

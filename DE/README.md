## RNAseq Differential Expression Analysis of _P. astreoides_ planulae and adult samples from shallow and mesophotic reefs ######

This script is based on the work of [Erin Chille](https://github.com/echille/Mcapitata_OA_Developmental_Gene_Expression_Timeseries/blob/main/4-Differential-Gene-Expression-Analysis/pHTreatment_RNAseqDE.Rmd), with some modifications and adjustments.

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

gPCAdata <- plotPCA(gvst, intgroup = c("group"), returnData=TRUE)
percentVar <- round(100*attr(gPCAdata, "percentVar")) #plot PCA of samples with all data
PCA <- ggplot(gPCAdata, aes(PC1, PC2, color=group)) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(labels = c("adult_meso","adult_shal","planu_meso", "planu_shal"), values = c("adult_meso"="blue", "adult_shal"="indianred3", "planu_meso"="deepskyblue", "planu_shal"="orange")) +
  coord_fixed() +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) + #Set the plot background
  theme(legend.position = ("top")); PCA #set title attributes`

ggsave(file = "PCA_all_vst.png", PCA)
```

#### Keep only selected life stages (adults) from treatment info and count data
```{r}
treatmentinfo_adult <- filter(treatmentinfo, age=="adult")
gcount_adult <- gcount[, treatmentinfo_adult$sample_id]

#create filter for the counts data
filt_adult <- filterfun(pOverA(0.5,10))
gfilt_adult <- genefilter(gcount_adult, filt_adult)

#identify genes to keep by count filter
gkeep_adult <- gcount_adult[gfilt_adult,]

#identify gene lists
gn.keep_adult <- rownames(gkeep_adult)

#gene count data filtered in PoverA, P percent of the samples have counts over A
gcount_filt_adult <- as.data.frame(gcount_adult[which(rownames(gcount_adult) %in% gn.keep_adult),])

### Construct the DESeq dataset

treatmentinfo_adult$depth <- factor(treatmentinfo_adult$depth, levels = c("mesophotic","shallow"))

#Set DESeq2 design
gdds_adult <- DESeqDataSetFromMatrix(countData = gcount_filt_adult,
                                  colData = treatmentinfo_adult,
                                  design = ~depth)
```

#### Visualize gene count data
```{r}
## Log-transform the count data
SF.gdds_adult <- estimateSizeFactors( gdds_adult) 
print(sizeFactors(SF.gdds_adult)) #need to be less than 4

gvst_adult <- vst(gdds_adult, blind=FALSE) 

### Principal component plot of adult samples

gPCAdata_adult <- plotPCA(gvst_adult, intgroup = c("depth"), returnData=TRUE)
percentVar_adult <- round(100*attr(gPCAdata_adult, "percentVar")) #plot PCA of samples with all data
PCA_adult <- ggplot(gPCAdata_adult, aes(PC1, PC2, color=depth)) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar_adult[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_adult[2],"% variance")) +
  scale_color_manual(labels = c("mesophotic", "shallow"), values = c("mesophotic"="blue", "shallow"="indianred3")) +
  coord_fixed() +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) + #Set the plot background
  theme(legend.position = ("top")); PCA_adult #set title attributes

ggsave(file = "Adult_PCA_vst.png", PCA_adult)
```

#### Keep only selected life stages (planulae) from treatment info and count data
```{r}
treatmentinfo_planulae <- filter(treatmentinfo, age=="planulae")
gcount_planulae <- gcount[, treatmentinfo_planulae$sample_id]

#create filter for the counts data
filt_planulae <- filterfun(pOverA(0.5,10))
gfilt_planulae <- genefilter(gcount_planulae, filt_planulae)

#identify genes to keep by count filter
gkeep_planulae <- gcount_planulae[gfilt_planulae,]

#identify gene lists
gn.keep_planulae <- rownames(gkeep_planulae)

#gene count data filtered in PoverA, P percent of the samples have counts over A
gcount_filt_planulae <- as.data.frame(gcount_planulae[which(rownames(gcount_planulae) %in% gn.keep_planulae),])

### Construct the DESeq dataset

treatmentinfo_planulae$depth <- factor(treatmentinfo_planulae$depth, levels = c("mesophotic","shallow"))

#Set DESeq2 design
gdds_planulae <- DESeqDataSetFromMatrix(countData = gcount_filt_planulae,
                                     colData = treatmentinfo_planulae,
                                     design = ~depth)
```

#### Visualize gene count data
```{r}
## Log-transform the count data
SF.gdds_planulae <- estimateSizeFactors( gdds_planulae) 
print(sizeFactors(SF.gdds_planulae)) 

gvst_planulae <- vst(gdds_planulae, blind=FALSE) 

#### Principal component plot of adult samples

gPCAdata_planulae <- plotPCA(gvst_planulae, intgroup = c("depth"), returnData=TRUE)
percentVar_planulae <- round(100*attr(gPCAdata_planulae, "percentVar")) #plot PCA of samples with all data
PCA_planulae <- ggplot(gPCAdata_planulae, aes(PC1, PC2, color=depth)) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar_planulae[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_planulae[2],"% variance")) +
  scale_color_manual(labels = c("mesophotic", "shallow"), values = c("mesophotic"="deepskyblue", "shallow"="orange")) +
  coord_fixed() +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) + #Set the plot background
  theme(legend.position = ("top")); PCA_planulae #set title attributes

ggsave(file = "planulae_PCA_vst.png", PCA_planulae)
```

#### Keep only selected depth (shallow) from treatment info and count data
```{r}
treatmentinfo_shallow <- filter(treatmentinfo, depth=="shallow")
gcount_shallow <- gcount[, treatmentinfo_shallow$sample_id]

#create filter for the counts data
filt_shallow <- filterfun(pOverA(0.5,10))
gfilt_shallow <- genefilter(gcount_shallow, filt_shallow)

#identify genes to keep by count filter
gkeep_shallow <- gcount_shallow[gfilt_shallow,]

#identify gene lists
gn.keep_shallow <- rownames(gkeep_shallow)

#gene count data filtered in PoverA, P percent of the samples have counts over A
gcount_filt_shallow <- as.data.frame(gcount_shallow[which(rownames(gcount_shallow) %in% gn.keep_shallow),])

### Construct the DESeq dataset
treatmentinfo_shallow$age <- factor(treatmentinfo_shallow$age, levels = c("adult","planulae"))

#Set DESeq2 design
gdds_shallow <- DESeqDataSetFromMatrix(countData = gcount_filt_shallow,
                                     colData = treatmentinfo_shallow,
                                     design = ~age)
```

#### Visualize gene count data
```{r}
## Log-transform the count data
SF.gdds_shallow <- estimateSizeFactors( gdds_shallow) 
print(sizeFactors(SF.gdds_shallow)) 

gvst_shallow <- vst(gdds_shallow, blind=FALSE) 

#### Principal component plot of shallow samples

gPCAdata_shallow <- plotPCA(gvst_shallow, intgroup = c("age"), returnData=TRUE)
percentVar_shallow <- round(100*attr(gPCAdata_shallow, "percentVar")) #plot PCA of samples with all data
PCA_shallow <- ggplot(gPCAdata_shallow, aes(PC1, PC2, color=age)) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar_shallow[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_shallow[2],"% variance")) +
  scale_color_manual(labels = c("adult", "planulae"), values = c("adult"="indianred3", "planulae"="orange")) +
  coord_fixed() +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) + #Set the plot background
  theme(legend.position = ("top")); PCA_shallow #set title attributes

ggsave(file = "shallow_PCA_vst.png", PCA_shallow)
```

#### Keep only selected depth (mesophotic) from treatment info and count data
```{r}
treatmentinfo_meso <- filter(treatmentinfo, depth=="mesophotic")
gcount_meso <- gcount[, treatmentinfo_meso$sample_id]

#create filter for the counts data
filt_meso <- filterfun(pOverA(0.5,10))
gfilt_meso <- genefilter(gcount_meso, filt_meso)

#identify genes to keep by count filter
gkeep_meso <- gcount_meso[gfilt_meso,]

#identify gene lists
gn.keep_meso <- rownames(gkeep_meso)

#gene count data filtered in PoverA, P percent of the samples have counts over A
gcount_filt_meso <- as.data.frame(gcount_meso[which(rownames(gcount_meso) %in% gn.keep_meso),])

### Construct the DESeq dataset
treatmentinfo_meso$age <- factor(treatmentinfo_meso$age, levels = c("adult","planulae"))

#Set DESeq2 design
gdds_meso <- DESeqDataSetFromMatrix(countData = gcount_filt_meso,
                                       colData = treatmentinfo_meso,
                                       design = ~age)
```

#### Visualize gene count data
```{r}
## Log-transform the count data
SF.gdds_meso <- estimateSizeFactors( gdds_meso) 
print(sizeFactors(SF.gdds_meso)) 

gvst_meso <- vst(gdds_meso, blind=FALSE) 

#### Principal component plot of mesophotic samples

gPCAdata_meso<- plotPCA(gvst_meso, intgroup = c("age"), returnData=TRUE)
percentVar_meso <- round(100*attr(gPCAdata_meso, "percentVar")) #plot PCA of samples with all data
PCA_meso <- ggplot(gPCAdata_meso, aes(PC1, PC2, color=age)) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar_meso[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_meso[2],"% variance")) +
  scale_color_manual(labels = c("adult", "planulae"), values = c("adult"="blue", "planulae"="deepskyblue")) +
  coord_fixed() +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) + #Set the plot background
  theme(legend.position = ("top")); PCA_meso #set title attributes

ggsave(file = "meso_PCA_vst.png", PCA_meso)
```

### Differential Gene Expression Analysis 

#### Run DE analysis

#DESEq2 itrnally applies the median of ratios method for normalization.
DEG_adult <- DESeq(gdds_adult) #run differential expression test by group using the Wald model

#Explore significant p-values for meso and shallow adults
DEG_adult.results <- results(DEG_adult, contrast= c("depth","mesophotic","shallow"))
write.csv(DEG_adult.results, "DEGs_meso_vs_shal_adults_ALL.csv")
head(DEG_adult.results)
sum(DEG_adult.results$padj < 0.05, na.rm=TRUE)  #439 genes
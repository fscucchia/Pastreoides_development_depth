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

#### Run DE analysis - Adult samples
```{r}
#DESEq2 internally applies the median of ratios method for normalization.
DEG_adult <- DESeq(gdds_adult) #run differential expression test by group using the Wald model

#Explore significant p-values for meso and shallow adults
DEG_adult.results <- results(DEG_adult, contrast= c("depth","mesophotic","shallow"))
write.csv(DEG_adult.results, "DEGs_meso_vs_shal_adults_ALL.csv")
head(DEG_adult.results)
sum(DEG_adult.results$padj < 0.05, na.rm=TRUE)  #439 genes

# Volcano Plot (significance as a function of fold change)
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(DEG_adult.results)
#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value (pajd))))
#Color the significant points with log2 fold change >2 red ()
#with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
with(subset(topT, padj<0.05), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))

#save list of DEGs
DEGs_adult <- as.data.frame(subset(DEG_adult.results, padj<0.05))
DEGs_adult_ordered <- order(DEGs_adult$padj) #Order p-values by smallest value first
DEGs_adult$contrast <- as.factor(c("adult_mesophotic_vs_shallow"))
DEGs_adult$gene_id  <- rownames(DEGs_adult)
rownames(DEGs_adult) <- NULL
write.csv(DEGs_adult, "DEGs_meso_vs_shal_adults.csv")
```
#### Plot DEGs 
```{r}
#identify signficant pvalues with 5%FDR
DEG_adult.results$gene_id  <- rownames(gdds_adult)
sig_adults <- subset(DEG_adult.results, padj<0.05,)
rownames(sig_adults) <- sig_adults[,7] #rename rownames of sig_adults as column 7
#subset list of sig transcripts from original count data
sig.list.adults <- gdds_adult[which(rownames(gdds_adult) %in% rownames(sig_adults)),]

#apply a vst transformation to minimize effects of small counts and normalize by library size
gdds_adult <- estimateSizeFactors(gdds_adult) 
print(sizeFactors(gdds_adult))
rsig <- vst(sig.list.adults, blind=FALSE)

PCA.sig.adults <- plotPCA(rsig, intgroup="depth")#Plot PCA of all samples for DEG only
PCA.sig.adults #view plot

pdf(file="PCA_sigDEGs_adults.pdf")
PCA.sig.adults #view plot
dev.off()

## Heatmap 
#make an expression object
#difference in expression compared to average across all samples
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(MetBrewer)

colors = met.brewer("OKeeffe1",n=15,type="continuous", direction=-1)

Sym.mat <- assay(rsig)
Sym.mat <- Sym.mat - rowMeans(Sym.mat)
Sym.df <- data.frame(colData(rsig)[c("depth")])
pdf(file="SigDEGs_Heatmap_adults_changeColors.pdf")
pheatmap(Sym.mat, annotation_col = Sym.df, clustering_method = "average",
         clustering_distance_rows="euclidean", 
         #color = viridis(400),
         color = colors,
         show_rownames =FALSE, cluster_cols=TRUE,
         show_colnames =TRUE) #plot heatmap of all DEG by group

dev.off()
```
#### Run DE analysis - Planulae samples
```{r}
DEG_planulae <- DESeq(gdds_planulae) #run differential expression test by group using the Wald model

DEG_planulae.results <- results(DEG_planulae, contrast= c("depth","mesophotic","shallow"))
write.csv(DEG_planulae.results, "DEGs_planulae_mesophotic_vs_shallow_ALL.csv")
head(DEG_planulae.results)
sum(DEG_planulae.results$padj < 0.05, na.rm=TRUE)  #2898 genes

#save DEGs list
DEGs_planulae <- as.data.frame(subset(DEG_planulae.results, padj<0.05))
DEGs_planulae_ordered <- order(DEGs_planulae$padj) #Order p-values by smallest value first
DEGs_planulae$contrast <- as.factor(c("planulae_mesophotic_vs_shallow"))
DEGs_planulae$gene_id  <- rownames(DEGs_planulae)
rownames(DEGs_planulae) <- NULL
write.csv(DEGs_planulae, "DEGs_meso_vs_shal_planulae.csv")

###Volcano plot 
# Volcano Plot (significance as a function of fold change)
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(DEG_planulae.results)
#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value (pajd))))
#Color the significant points 
with(subset(topT, padj<0.05), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
```
#### Plot DEGs 
```{r}
#identify signficant pvalues with 5%FDR
DEG_planulae.results$gene_id  <- rownames(gdds_planulae)
sig_planulae <- subset(DEG_planulae.results, padj<0.05,)
rownames(sig_planulae) <- sig_planulae[,7] #rename rownames of sig_planulae as column 7
#subset list of sig transcripts from original count data
sig.list.planulae <- gdds_planulae[which(rownames(gdds_planulae) %in% rownames(sig_planulae)),]


#apply a vst transformation to minimize effects of small counts and normalize wrt library size
gdds_planulae <- estimateSizeFactors(gdds_planulae) 
print(sizeFactors(gdds_planulae)) 
rsig <- vst(sig.list.planulae, blind=FALSE)

PCA.sig.planulae <- plotPCA(rsig, intgroup="depth")#Plot PCA of all samples for DEG only
PCA.sig.planulae #view plot

pdf(file="PCA_sigDEGs_planulae.pdf")
PCA.sig.planulae #view plot
dev.off()

### Heatmap 
#make an expression object
#difference in expression compared to average across all samples
library(pheatmap)
Sym.mat <- assay(rsig)
Sym.mat <- Sym.mat - rowMeans(Sym.mat)
Sym.df <- data.frame(colData(rsig)[c("depth")])
colors = met.brewer("OKeeffe1",n=20,type="continuous", direction=-1)

pdf(file="SigDEGs_Heatmap_planulae_changeColors.pdf")
pheatmap(Sym.mat, annotation_col = Sym.df, clustering_method = "average",
         clustering_distance_rows="euclidean", 
         color=colors,
         show_rownames =FALSE, cluster_cols=TRUE,
         show_colnames =TRUE) #plot heatmap of all DEG by group

dev.off()
```

#### Run DE analysis - Shallow samples
```{r}
DEG_shallow <- DESeq(gdds_shallow) #run differential expression test by group using the Wald model

DEG_shallow.results <- results(DEG_shallow, contrast= c("age","planulae","adult"))  #adult is the reference
write.csv(DEG_shallow.results, "DEGs_shal_planulae_vs_adults_ALL.csv")
head(DEG_shallow.results)
sum(DEG_shallow.results$padj < 0.05, na.rm=TRUE)  #10664 genes

#save list of DEGs
DEGs_shallow <- as.data.frame(subset(DEG_shallow.results, padj<0.05))
DEGs_shallow_ordered <- order(DEGs_shallow$padj) #Order p-values by smallest value first
DEGs_shallow$contrast <- as.factor(c("shallow_planulae_vs_adults"))
DEGs_shallow$gene_id  <- rownames(DEGs_shallow)
rownames(DEGs_shallow) <- NULL
write.csv(DEGs_shallow, "DEGs_shal_planulae_vs_adults.csv")

### Volcano plot 
# Volcano Plot (significance as a function of fold change)
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(DEG_shallow.results)
#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value (pajd))))
#Color the significant points 
with(subset(topT, padj<0.05), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
```
#### Plot DEGs 
```{r}
#identify signficant pvalues with 5%FDR
DEG_shallow.results$gene_id  <- rownames(gdds_shallow)
sig_shallow <- subset(DEG_shallow.results, padj<0.05,)
rownames(sig_shallow) <- sig_shallow[,7] #rename rownames of sig_shallow as column 7
#subset list of sig transcripts from original count data
sig.list.shallow <- gdds_shallow[which(rownames(gdds_shallow) %in% rownames(sig_shallow)),]


#apply a vst transformation to minimize effects of small counts and normalize wrt library size
gdds_shallow <- estimateSizeFactors(gdds_shallow) 
print(sizeFactors(gdds_shallow)) 
rsig <- vst(sig.list.shallow, blind=FALSE)

PCA.sig.shallow <- plotPCA(rsig, intgroup="age")#Plot PCA of all samples for DEG only
PCA.sig.shallow #view plot

pdf(file="PCA_sigDEGs_shallow.pdf")
PCA.sig.shallow #view plot
dev.off()

### Heatmap 
#make an expression object
#difference in expression compared to average across all samples
library(pheatmap)
Sym.mat <- assay(rsig)
Sym.mat <- Sym.mat - rowMeans(Sym.mat)
Sym.df <- data.frame(colData(rsig)[c("age")])
colors = met.brewer("OKeeffe1",n=15,type="continuous", direction=-1)

pdf(file="SigDEGs_Heatmap_shallow_changeColors.pdf")
pheatmap(Sym.mat, annotation_col = Sym.df, clustering_method = "average",
         clustering_distance_rows="euclidean", 
         color=colors,
         show_rownames =FALSE, cluster_cols=TRUE,
         show_colnames =TRUE) #plot heatmap of all DEG by group

dev.off()
```

#### Run DE analysis - Mesophotic samples
```{r}
DEG_meso <- DESeq(gdds_meso) #run differential expression test by group using the Wald model

DEG_meso.results <- results(DEG_meso, contrast= c("age","planulae","adult"))  #adult is the reference
write.csv(DEG_meso.results, "DEGs_meso_planulae_vs_adults_ALL.csv")
head(DEG_meso.results)
sum(DEG_meso.results$padj < 0.05, na.rm=TRUE)  #8257 genes

#save list of DEGs
DEGs_meso <- as.data.frame(subset(DEG_meso.results, padj<0.05))
DEGs_meso_ordered <- order(DEGs_meso$padj) #Order p-values by smallest value first
DEGs_meso$contrast <- as.factor(c("meso_planulae_vs_adults"))
DEGs_meso$gene_id  <- rownames(DEGs_meso)
rownames(DEGs_meso) <- NULL
write.csv(DEGs_meso, "DEGs_meso_planulae_vs_adults.csv")

###Volcano plot 
# Volcano Plot (significance as a function of fold change)
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(DEG_meso.results)
#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value (pajd))))
#Color the significant points 
with(subset(topT, padj<0.05), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
```
#### Plot DEGs 
```{r}
###### Plot DEGs 
#identify signficant pvalues with 5%FDR
DEG_meso.results$gene_id  <- rownames(gdds_meso)
sig_meso <- subset(DEG_meso.results, padj<0.05,)
rownames(sig_meso) <- sig_meso[,7] #rename rownames of sig_meso as column 7
#subset list of sig transcripts from original count data
sig.list.meso <- gdds_meso[which(rownames(gdds_meso) %in% rownames(sig_meso)),]


#apply a vst transformation to minimize effects of small counts and normalize wrt library size
gdds_meso <- estimateSizeFactors(gdds_meso) 
print(sizeFactors(gdds_meso)) 
rsig <- vst(sig.list.meso, blind=FALSE)

PCA.sig.meso <- plotPCA(rsig, intgroup="age")#Plot PCA of all samples for DEG only
PCA.sig.meso #view plot

pdf(file="PCA_sigDEGs_mesophotic.pdf")
PCA.sig.meso #view plot
dev.off()

### Heatmap 
#make an expression object
#difference in expression compared to average across all samples
library(pheatmap)
Sym.mat <- assay(rsig)
Sym.mat <- Sym.mat - rowMeans(Sym.mat)
Sym.df <- data.frame(colData(rsig)[c("age")])
pdf(file="SigDEGs_Heatmap_mesophotic_changeColors.pdf")
pheatmap(Sym.mat, annotation_col = Sym.df, clustering_method = "average",
         clustering_distance_rows="euclidean", 
         color=colors,
         show_rownames =FALSE, cluster_cols=TRUE,
         show_colnames =TRUE) #plot heatmap of all DEG by group

dev.off()


#### merge all DEGs tables
DEGs_all <- bind_rows(DEGs_adult, DEGs_planulae, DEGs_shallow, DEGs_meso)
```

### Visualize differentially-expressed genes: ADULTS
```{r}
#### Subset and Log-transform the count data
#Subset the gene count matrix by the list of DEGs
DEG_adult.results$gene_id  <- rownames(gdds_adult)
sig_adults <- subset(DEG_adult.results, padj<0.05,)
rownames(sig_adults) <- sig_adults[,7] #rename rownames of sig_adults as column 7
#subset list of sig transcripts from original count data
sig.list.adults <- gdds_adult[which(rownames(gdds_adult) %in% rownames(sig_adults)),]

#apply a rlog transformation to minimize effects of small counts and normalize wrt library size
rsig <- rlog(sig.list.adults, blind=FALSE)


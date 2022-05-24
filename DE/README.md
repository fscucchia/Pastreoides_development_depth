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
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot"))
with(subset(topT , padj<.05), points(log2FoldChange, -log10(padj), pch=20,col="black"))
with(subset(topT, padj<0.05 & (log2FoldChange)>0), points(log2FoldChange, -log10(padj), pch=20, col="violetred2",cex=0.5))
with(subset(topT, padj<0.05 & (log2FoldChange)<(0*-0)), points(log2FoldChange, -log10(padj), pch=20, col="lightseagreen",cex=0.5))

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
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot"))
with(subset(topT , padj<.05), points(log2FoldChange, -log10(padj), pch=20,col="black"))
with(subset(topT, padj<0.05 & (log2FoldChange)>0), points(log2FoldChange, -log10(padj), pch=20, col="violetred2",cex=0.5))
with(subset(topT, padj<0.05 & (log2FoldChange)<(0*-0)), points(log2FoldChange, -log10(padj), pch=20, col="lightseagreen",cex=0.5))
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

#### merge all DEGs tables
```{r}
DEGs_all <- bind_rows(DEGs_adult, DEGs_planulae)
```

#### Visualize differentially-expressed genes: ADULTS
```{r}
##### Subset and Log-transform the count data
#Subset the gene count matrix by the list of DEGs
DEG_adult.results$gene_id  <- rownames(gdds_adult)
sig_adults <- subset(DEG_adult.results, padj<0.05,)
rownames(sig_adults) <- sig_adults[,7] #rename rownames of sig_adults as column 7
#subset list of sig transcripts from original count data
sig.list.adults <- gdds_adult[which(rownames(gdds_adult) %in% rownames(sig_adults)),]

#apply a vst transformation to minimize effects of small counts and normalize wrt library size
rsig <- vst(sig.list.adults, blind=FALSE)

#Make a matrix for computing similarity
mat.adults <- assay(rsig)#[pln_topVarGenes, ] #make an expression object
mat.adults <- mat.adults - rowMeans(mat.adults) #difference in expression compared to average across all samples

### Compute the optimal number of clusters for plotting
#Find the optimum number of clusters using 30 indexes with the NbClust() package. 
#This took about 35 minutes to run, so it is commented out.

nb <- NbClust(mat.adults, distance = "euclidean", min.nc = 2,
         max.nc = 10, method = "kmeans")
#According to the majority rule, the best number of clusters is  2 
 
fviz_nbclust(nb)

calc.kmeans <- kmeans(mat.adults, 2)
cluster_res <- data.frame(gene_id = names(calc.kmeans$cluster), cluster = calc.kmeans$cluster)
# 

DEGs_adult <- as.data.frame(subset(DEG_adult.results, padj<0.05))
write.csv(DEGs_adult, "DEGs_adult_withoutCluster.csv")
DEGs_adult_ordered <- order(DEGs_adult$padj) #Order p-values by smallest value first
DEGs_adult$contrast <- as.factor(c("adult_mesophotic_vs_shallow"))

DEGs_adult_cluster <- merge(DEGs_adult, cluster_res, by = "gene_id")
write.csv(DEGs_adult_cluster, "DEGs_adult_cluster.csv")


#### Plot a heatmap of differentially-expressed genes

DEGs_adult_cluster_heat <- read.csv("DEGs_adult_cluster.csv", header = TRUE, sep = ",")[,-c(1)]
DEGs_adult_cluster_heat <- subset(DEGs_adult_cluster, select = c(gene_id, cluster))

#Prepare annotations
hm_ann_row <- unique(DEGs_adult_cluster_heat)
rownames(hm_ann_row) <- hm_ann_row$gene_id
hm_ann_row <- subset(hm_ann_row, select=cluster)
hm_ann_row$cluster <- gsub(1,"Cluster1",hm_ann_row$cluster)
hm_ann_row$cluster <- gsub(2,"Cluster2",hm_ann_row$cluster)
hm_ann_row <- as.matrix(hm_ann_row[rownames(mat.adults),])
hmdepth_adults <- colData(grlog_adult)[c("depth")]

hmdepth_adults$depth <- factor(hmdepth_adults$depth, levels=c("shallow", "mesophotic"))
hm_ann_col <- HeatmapAnnotation(df=hmdepth_adults, col = list(depth=c("mesophotic" ="cadetblue", "shallow"  ="indianred3"))) #make dataframe for column naming

DEGs_adult_clusterHeatmap <-  Heatmap(mat.adults, column_title = "Depth", 
                           name = "expression",
                           show_row_names = FALSE, top_annotation = hm_ann_col, show_column_names = FALSE, row_dend_side = "left" ,
                           column_split = 2, column_dend_height = unit(0.5, "in"),
                           km = 2, row_km_repeats = 100, row_title = c("Cluster1", "Cluster2"),
                           row_gap = unit(2.5, "mm"), border = TRUE,
                           column_names_gp =  gpar(fontsize = 10)); DEGs_adult_clusterHeatmap
png("DEGs_adult_clusterHeatmap.png")
DEGs_adult_clusterHeatmap
dev.off()
```

#### Visualize differentially-expressed genes: PLANULAE
```{r}
##### Subset and Log-transform the count data
#Subset the gene count matrix by the list of DEGs
DEG_planulae.results$gene_id  <- rownames(gdds_planulae)
sig_planulae <- subset(DEG_planulae.results, padj<0.05,)
rownames(sig_planulae) <- sig_planulae[,8] #rename rownames of sig_planulae as column 8 (gene_ID)
#subset list of sig transcripts from original count data
sig.list.planulae <- gdds_planulae[which(rownames(gdds_planulae) %in% rownames(sig_planulae)),]

#apply a vst transformation to minimize effects of small counts and normalize wrt library size
planulae_sig <- vst(sig.list.planulae, blind=FALSE)

#Make a matrix for computing similarity
mat.planulae <- assay(planulae_sig)#[pln_topVarGenes, ] #make an expression object
mat.planulae <- mat.planulae - rowMeans(mat.planulae) #difference in expression compared to average across all samples

### Compute the optimal number of clusters for plotting
#Find the optimum number of clusters using 30 indexes with the NbClust() package. 
#This took about 35 minutes to run, so it is commented out.

nb_planulae <- NbClust(mat.planulae, distance = "euclidean", min.nc = 2,
              max.nc = 10, method = "kmeans")

#
lista.methods = c("kl", "ch", "hartigan","mcclain", "gamma", "gplus",
                  "tau", "dunn", "sdindex", "sdbw", "cindex", "silhouette",
                  "ball","ptbiserial", "gap","frey")
lista.distance = c("metodo","euclidean", "maximum", "manhattan", "canberra")

tabla = as.data.frame(matrix(ncol = length(lista.distance), nrow = length(lista.methods)))
names(tabla) = lista.distance

for (j in 2:length(lista.distance)){
  for(i in 1:length(lista.methods)){
    
    nb_planulae = NbClust(mat.planulae, distance = lista.distance[j],
                 min.nc = 2, max.nc = 10, 
                 method = "kmeans", index =lista.methods[i])
    tabla[i,j] = nb_planulae$Best.nc[1]
    tabla[i,1] = lista.methods[i]
    
  }}

tabla
#

#According to the majority rule, the best number of clusters is  2 

fviz_nbclust(nb_planulae)

calc.kmeans <- kmeans(mat.planulae, 2)
cluster_res <- data.frame(gene_id = names(calc.kmeans$cluster), cluster = calc.kmeans$cluster)
# 

DEGs_planulae <- as.data.frame(subset(DEG_planulae.results, padj<0.05))
write.csv(DEGs_planulae, "DEGs_planulae_withoutCluster.csv")
DEGs_planulae_ordered <- order(DEGs_planulae$padj) #Order p-values by smallest value first
DEGs_planulae$contrast <- as.factor(c("planulae_mesophotic_vs_shallow"))

DEGs_planulae_cluster <- merge(DEGs_planulae, cluster_res, by = "gene_id")
write.csv(DEGs_planulae_cluster, "DEGs_planulae_cluster.csv")


#### Plot a heatmap of differentially-expressed genes

DEGs_planulae_cluster_heat <- read.csv("DEGs_planulae_cluster.csv", header = TRUE, sep = ",")[,-c(1)]
DEGs_planulae_cluster_heat <- subset(DEGs_planulae_cluster, select = c(gene_id, cluster))

#Prepare annotations
hm_ann_row <- unique(DEGs_planulae_cluster_heat)
rownames(hm_ann_row) <- hm_ann_row$gene_id
hm_ann_row <- subset(hm_ann_row, select=cluster)
hm_ann_row$cluster <- gsub(1,"Cluster1",hm_ann_row$cluster)
hm_ann_row$cluster <- gsub(2,"Cluster2",hm_ann_row$cluster)
hm_ann_row <- as.matrix(hm_ann_row[rownames(mat.planulae),])
hmdepth_planulae <- colData(gvst_planulae)[c("depth")]

hmdepth_planulaedepth <- factor(hmdepth_planulae$depth, levels=c("shallow", "mesophotic"))
hm_ann_col <- HeatmapAnnotation(df=hmdepth_planulae, col = list(depth=c("mesophotic" ="cadetblue", "shallow"  ="indianred3"))) #make dataframe for column naming

DEGs_planulae_clusterHeatmap <-  Heatmap(mat.planulae, column_title = "Depth", 
                                      name = "expression",
                                      show_row_names = FALSE, top_annotation = hm_ann_col, show_column_names = FALSE, row_dend_side = "left" ,
                                      column_split = 2, column_dend_height = unit(0.5, "in"),
                                      km = 2, row_km_repeats = 100, row_title = c("Cluster1", "Cluster2"),
                                      row_gap = unit(2.5, "mm"), border = TRUE,
                                      column_names_gp =  gpar(fontsize = 10)); DEGs_planulae_clusterHeatmap
png("DEGs_planulae_clusterHeatmap.png")
DEGs_planulae_clusterHeatmap
dev.off()
```



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

# Cluster the samples to look for obvious outliers
#Look for outliers by examining the sample tree:

sampleTree = hclust(dist(datExpr), method = "average")

# Plot the sample tree
pdf(paste0('sampleTree','.pdf'))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

### Network construction and consensus module detection
# Choosing a soft-thresholding power: Analysis of network topology 
# The soft thresholding power is the number to which the co-expression similarity is raised to calculate adjacency. 

#Choose a set of soft-thresholding powers
powers <- c(seq(from = 1, to=19, by=2), c(21:30)) #Create a string of numbers from 1 through 10, and even numbers from 10 through 20
powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 2))

#Call the network topology analysis function
sft <-pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

#Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9;
#Scale-free topology fit index as a function of the soft-thresholding power
pdf(paste0('network','.pdf'))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# # # this line corresponds to using an R^2 cut-off
abline(h=0.8,col="red")
# # # Mean connectivity as a function of the soft-thresholding power

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#The lowest scale-free topology fit index R^2 recommended by Langfelder and Horvath is 0.8. 
#From the graph, it appears that our soft thresholding power is 17 because it is the lowest 
#power before the R^2=0.8 threshold that maximizes with model fit (s17 is the number right above the red line).

### Network construction and module detection:
# Co-expression adjacency and topological overlap matrix similarity
# Co-expression similarity and adjacency, using the soft thresholding power 17 and translate the adjacency into topological overlap matrix to calculate 
# the corresponding dissimilarity. I will use a signed network because we have a relatively high softPower, according 
# to >12 (https://peterlangfelder.com/2018/11/25/__trashed/). 
# Moreover, in expression data where you are interested in when expression on one gene increases or decreases with expression level of another you would use a signed network (when you are interested in the direction of change, correlation and anti-correlation, you use a signed network).

options(stringsAsFactors = FALSE)
enableWGCNAThreads() #Allow multi-threading within WGCNA

#Run analysis
softPower=17 #Set softPower to 17
adjacency=adjacency(datExpr, power=softPower,type="signed") #Calculate adjacency
TOM= TOMsimilarity(adjacency,TOMType = "signed") #Translate adjacency into topological overlap matrix
dissTOM= 1-TOM #Calculate dissimilarity in TOM
save(adjacency, TOM, dissTOM, file = "/data/home/mass/fscucchia/Bermuda/output/WGCNA/adjTOM.RData")
save(dissTOM, file = "/data/home/mass/fscucchia/Bermuda/output/WGCNA/dissTOM.RData") 

# Clustering using TOM
#Form distance matrix
geneTree= flashClust(as.dist(dissTOM), method="average")

#We will now plot a dendrogram of genes. Each leaf corresponds to a gene, branches grouping together densely are interconnected, highly co-expressed genes
pdf(file="/data/home/mass/fscucchia/Bermuda/output/WGCNA/dissTOMClustering.pdf", width=20, height=20)
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)
dev.off()
```

### Modules identification
```
# Module identification is cutting the branches off the tree in the dendrogram above. We want large modules, so we set the minimum module size 
# relatively high (minimum size = 30).

minModuleSize = 30 #default value used most often
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize)
table(dynamicMods) #list modules and respective sizes
save(dynamicMods, geneTree, file = "/data/home/mass/fscucchia/Bermuda/output/WGCNA/dyMod_geneTree.RData")

dyMod_geneTree <- load(file = "/data/home/mass/fscucchia/pH_exp_1/WGCNA/dyMod_geneTree.RData")

dyMod_geneTree
## [1] "dynamicMods" "geneTree"

# Plot the module assignment under the gene dendrogram
dynamicColors = labels2colors(dynamicMods) # Convert numeric labels into colors
table(dynamicColors)
# [1] dynamicColors
# antiquewhite2   antiquewhite4         bisque4           black            blue
#              58              97             121            1080            3576
#           blue2           brown          brown2          brown4           coral
#              74            3006              75             127              58
#          coral1          coral2          coral3            cyan       darkgreen
#             104              95              56             582             321
#        darkgrey     darkmagenta  darkolivegreen darkolivegreen2 darkolivegreen4
#             290             185             186              30              75
#      darkorange     darkorange2         darkred   darkseagreen3   darkseagreen4
#             264             134             347              59             105
#   darkslateblue   darkturquoise      darkviolet      firebrick3      firebrick4
#             117             317              74              31              78
#     floralwhite           green     greenyellow          grey60        honeydew
#             139            1506             662             507              59
#       honeydew1      indianred3      indianred4           ivory  lavenderblush2
#             105              35              80             140              60
#  lavenderblush3      lightblue4      lightcoral       lightcyan      lightcyan1
#             108              41              80             512             143
#      lightgreen      lightpink3      lightpink4  lightslateblue  lightsteelblue
#             470              61             110              43              87
# lightsteelblue1     lightyellow         magenta        magenta4          maroon
#             157             461             768              61             111
#    mediumorchid   mediumpurple1   mediumpurple2   mediumpurple3   mediumpurple4
#              95              47              92             163              55
#    midnightblue    navajowhite1    navajowhite2          orange      orangered1
#             562              63             112             267              49
#      orangered3      orangered4   paleturquoise  palevioletred2  palevioletred3
#              93             165             208              66             112
#            pink           pink4            plum           plum1           plum2
#             815              49              93             166             115
#           plum3          purple             red       royalblue     saddlebrown
#              73             693            1084             393             233
#          salmon         salmon2         salmon4         sienna3         sienna4
#             607              69             113             179              50
#         skyblue        skyblue1        skyblue2        skyblue3        skyblue4
#             246              93              95             170              53
#       steelblue             tan         thistle        thistle1        thistle2
#             223             653              70             113             114
#        thistle3       turquoise          violet           white          yellow
#              72            5144             190             251            2864
#         yellow3         yellow4     yellowgreen
#              52              94             170





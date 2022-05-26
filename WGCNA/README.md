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

pdf(file="/data/home/mass/fscucchia/Bermuda/output/WGCNA/dissTOMColorClustering.pdf", width=20, height=20)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

# Merge modules whose expression profiles are very similar or choose not to merge
# Plot module similarity based on eigengene value

#Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors, softPower = 17)
MEs = MEList$eigengenes

#Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

#Cluster again and plot the results
METree = flashClust(as.dist(MEDiss), method = "average")

pdf(file="/data/home/mass/fscucchia/Bermuda/output/WGCNA/eigengeneClustering1.pdf", width = 20)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()

#Merge modules with >85% eigengene similarity (most studies use80-90% similarity)

MEDissThres= 0.15 #merge modules that are 85% similar

pdf(file="/data/home/mass/fscucchia/Bermuda/output/WGCNA/eigengeneClustering2.pdf", width = 20)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h=MEDissThres, col="red")
dev.off()

merge= mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
# [1] mergeCloseModules: Merging modules whose distance is less than 0.15
#    multiSetMEs: Calculating module MEs.
#      Working on set 1 ...
#      moduleEigengenes: Calculating 103 module eigengenes in given set.
#    multiSetMEs: Calculating module MEs.
#      Working on set 1 ...
#      moduleEigengenes: Calculating 61 module eigengenes in given set.
#    multiSetMEs: Calculating module MEs.
#      Working on set 1 ...
#      moduleEigengenes: Calculating 56 module eigengenes in given set.
#    Calculating new MEs...
#    multiSetMEs: Calculating module MEs.
#      Working on set 1 ...
#      moduleEigengenes: Calculating 56 module eigengenes in given set.

mergedColors= merge$colors
mergedMEs= merge$newMEs

pdf(file="/data/home/mass/fscucchia/Bermuda/output/WGCNA/mergedClusters.pdf", width=20, height=20)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
dev.off()

#Save new colors

moduleColors = mergedColors # Rename to moduleColors
colorOrder = c("grey", standardColors(50)); # Construct numerical labels corresponding to the colors
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
ncol(MEs) 
# [1] 56
#Instead of 103 modules, we now have 56, a much more reasonable number.

# Plot new tree
#Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
#Cluster again and plot the results
pdf(file="/data/home/mass/fscucchia/Bermuda/output/WGCNA/eigengeneClustering3.pdf")
METree = flashClust(as.dist(MEDiss), method = "average")
MEtreePlot = plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()

# Relating modules to group, quantifying module–trait associations
#Prepare trait data. Data has to be numeric, so I replaced the group for numeric values

allTraits <- names(treatmentinfo$group)
allTraits$adult_meso <- c(1,1,1,0,0,0,0,0,0,0,0,0)
allTraits$planu_meso <- c(0,0,0,1,1,1,0,0,0,0,0,0)
allTraits$adult_shal <- c(0,0,0,0,0,0,1,1,1,0,0,0)
allTraits$planu_shal <- c(0,0,0,0,0,0,0,0,0,1,1,1)

datTraits <- as.data.frame(allTraits)
dim(datTraits)
## [1] 12  4

rownames(datTraits) <- treatmentinfo$sample_id
print(datTraits)
# [1]     adult_meso planu_meso adult_shal planu_shal
# D20           1          0          0          0
# D21           1          0          0          0
# D22           1          0          0          0
# L7            0          1          0          0
# L8            0          1          0          0
# L9            0          1          0          0
# R13B          0          0          1          0
# R15B          0          0          1          0
# R33           0          0          1          0
# S4            0          0          0          1
# S5            0          0          0          1
# S7            0          0          0          1

#Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors,softPower=17)$eigengenes
MEs = orderMEs(MEs0)
names(MEs)
# [1] "MEskyblue2"        "MElavenderblush2"  "MEviolet"
#  [4] "MEdarkseagreen4"   "MEmediumpurple4"   "MEdarkolivegreen2"
#  [7] "MEindianred4"      "MEthistle2"        "MEsienna4"
# [10] "MEthistle"         "MEskyblue4"        "MEcoral2"
# [13] "MEdarkmagenta"     "MEbrown4"          "MEdarkolivegreen"
# [16] "MEplum2"           "MEfirebrick3"      "MEmediumpurple3"
# [19] "MEantiquewhite2"   "MEmediumpurple1"   "MEdarkolivegreen4"
# [22] "MEyellowgreen"     "MElightcyan1"      "MEdarkgreen"
# [25] "MEdarkslateblue"   "MEthistle3"        "MEbisque4"
# [28] "MEcoral3"          "MEcoral1"          "MEturquoise"
# [31] "MEhoneydew"        "MEplum"            "MEtan"
# [34] "MEdarkgrey"        "MEgrey60"          "MElightsteelblue"
# [37] "MEorangered4"      "MEindianred3"      "MEbrown2"
# [40] "MEpaleturquoise"   "MEmaroon"          "MEpalevioletred3"
# [43] "MElightgreen"      "MEyellow3"         "MEblack"
# [46] "MEdarkorange"      "MEcyan"            "MElightcyan"
# [49] "MEfloralwhite"     "MEantiquewhite4"   "MEplum3"
# [52] "MEdarkviolet"      "MEyellow4"         "MEnavajowhite2"
# [55] "MEblue2"           "MEorangered3"
```

### Correlations of traits and genes with eigengenes
```
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
Colors=sub("ME","",names(MEs))

moduleTraitTree = hclust(dist(t(moduleTraitCor)), method = "average");
pdf(file="/data/home/mass/fscucchia/Bermuda/output/WGCNA/Life-stage-depth clustering based on module-trait correlation.pdf")
plot(moduleTraitTree, main = "Group clustering based on module-trait correlation", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

#Correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datExpr)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples);
```

### Plot module-trait associations
```
#Plot module trait correlations as a heatmap

textMatrix = paste(signif(moduleTraitCor, 2), "\n(" ,signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
head(textMatrix)

pdf(file="/data/home/mass/fscucchia/Bermuda/output/WGCNA/Module-trait-relationships.pdf")
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits),  yLabels = names(MEs), ySymbols = names(MEs), 
     cex.lab.y= 0.55, cex.lab.x= 0.55, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = TRUE, cex.text = 0.25, textAdj = , 
     zlim = c(-1,1), main = paste("Module-trait relationships"))

labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits),  yLabels = names(MEs), ySymbols = names(MEs), cex.lab.y= 0.55, 
     cex.lab.x= 0.55, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = TRUE, cex.text = 0.25, textAdj = , zlim = c(-1,1), 
     main = paste("Module-trait relationships"))
dev.off()


# Plot as clustered Heatmap
#add bold sigignificant p-values, dendrogram with WGCNA MEtree cut-off, module clusters

#Create list of pvalues for eigengene correlation with specific life stages
heatmappval <- signif(moduleTraitPvalue, 1)

#Make list of heatmap row colors
htmap.colors <- names(MEs)
htmap.colors <- gsub("ME", "", htmap.colors)

library(dendsort)
row_dend = dendsort(hclust(dist(moduleTraitCor)))
col_dend = dendsort(hclust(dist(t(moduleTraitCor))))

pdf(file = "/data/home/mass/fscucchia/Bermuda/output/WGCNA/Module-trait-relationship-heatmap3.pdf", height = 11.5, width = 8)
ht=Heatmap(moduleTraitCor, name = "Eigengene", column_title = "Module-Trait Eigengene Correlation", 
        col = blueWhiteRed(50), 
        row_names_side = "left", row_dend_side = "left",
        width = unit(4, "in"), height = unit(8.5, "in"), 
        column_order = 1:4, column_dend_reorder = FALSE, cluster_columns = hclust(dist(t(moduleTraitCor)), method = "average"), column_split = 3, column_dend_height = unit(0.5, "in"),
        cluster_rows = METree, row_split = 10, row_gap = unit(2.5, "mm"), border = TRUE,
        cell_fun = function(j, i, x, y, w, h, col) {
        if(heatmappval[i, j] <= 0.05) {
            grid.text(sprintf("%s", heatmappval[i, j]), x, y, gp = gpar(fontsize = 8, fontface = "bold"))
        }
        else {
            grid.text(sprintf("%s", heatmappval[i, j]), x, y, gp = gpar(fontsize = 8, fontface = "plain"))
        }},
        column_names_gp =  gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10, alpha = 0.75, border = TRUE, fill = htmap.colors))
draw(ht)
dev.off()

# Create dataframe that associates module colors with clusters based on the above heatmap and the MEtree dendrogram.
MEcluster1 <- data.frame(moduleColor = c("skyblue2","lavenderblush2", "violet", "darkseagreen4", "mediumpurple4"), moduleCluster = c(1))
MEcluster2 <- data.frame(moduleColor = c("darkolivegreen2", "indianred4", "thistle2", "sienna4", "thistle"), moduleCluster = c(2))
MEcluster3 <- data.frame(moduleColor = c("skyblue4","coral2", "darkmagenta", "brown4", "darkolivegreen", "plum2",  "firebrick3", "mediumpurple3"), moduleCluster = c(3))
MEcluster4 <- data.frame(moduleColor = c("antiquewhite2", "mediumpurple1", "darkolivegreen4","yellowgreen", "lightcyan1","darkgreen", "darkslateblue"), moduleCluster = c(4))
MEcluster5 <- data.frame(moduleColor = c("thistle3", "bisque4", "coral3","coral1", "turquoise"), moduleCluster = c(5))
MEcluster6 <- data.frame(moduleColor = c("honeydew", "plum" , "tan", "darkgrey", "grey60", "lightsteelblue", "orangered4" ), moduleCluster = c(6)) 
MEcluster7 <- data.frame(moduleColor = c("indianred3","brown2","paleturquoise" ), moduleCluster = c(7)) 
MEcluster8 <- data.frame(moduleColor = c("maroon", "palevioletred3", "lightgreen", "yellow3", "black", "darkorange", "cyan", "lightcyan"), moduleCluster = c(8)) 
MEcluster9 <- data.frame(moduleColor = c("floralwhite", "antiquewhite4", "plum3"), moduleCluster = c(9)) 
MEcluster10 <- data.frame(moduleColor = c("darkviolet", "yellow4", "navajowhite2", "blue2","orangered3"), moduleCluster = c(10)) 

moduleCluster = bind_rows(MEcluster1, MEcluster2, MEcluster3, MEcluster4, MEcluster5, MEcluster6, MEcluster7, MEcluster8, MEcluster9, MEcluster10)

head(moduleCluster)
# [1]     moduleColor moduleCluster
# 1        skyblue2             1
# 2  lavenderblush2             1
# 3          violet             1
# 4   darkseagreen4             1
# 5   mediumpurple4             1
# 6 darkolivegreen2             2
```

### Make dataframe for Strader plots
```
# View module eigengene data
head(MEs)
#[1] MEskyblue2 MElavenderblush2    MEviolet MEdarkseagreen4 MEmediumpurple4
# D20  0.4272220        0.3072634  0.77585753       0.8163653      0.60269797
# D21 -0.2428876        0.2638944  0.11605634      -0.1314675     -0.29729763
# D22 -0.2641934        0.2760190 -0.07048764       0.0401772      0.08470422
# L7  -0.1604286       -0.2992337 -0.14998878      -0.2713475     -0.16581250
# L8  -0.1368332       -0.3111931 -0.17572893      -0.3692494     -0.25930120
# L9  -0.2145480       -0.2977739 -0.27233310      -0.1800915     -0.22182692

names(MEs)
# [1] "MEskyblue2"        "MElavenderblush2"  "MEviolet"
#  [4] "MEdarkseagreen4"   "MEmediumpurple4"   "MEdarkolivegreen2"
#  [7] "MEindianred4"      "MEthistle2"        "MEsienna4"
# [10] "MEthistle"         "MEskyblue4"        "MEcoral2"
# [13] "MEdarkmagenta"     "MEbrown4"          "MEdarkolivegreen"    

Strader_MEs <- MEs
Strader_MEs$group <- treatmentinfo$group
Strader_MEs$sample_id <- rownames(Strader_MEs)
head(Strader_MEs)
# [1] MEskyblue2 MElavenderblush2    MEviolet MEdarkseagreen4 MEmediumpurple4
# D20  0.4272220        0.3072634  0.77585753       0.8163653      0.60269797
# D21 -0.2428876        0.2638944  0.11605634      -0.1314675     -0.29729763
# D22 -0.2641934        0.2760190 -0.07048764       0.0401772      0.08470422
# L7  -0.1604286       -0.2992337 -0.14998878      -0.2713475     -0.16581250
# L8  -0.1368332       -0.3111931 -0.17572893      -0.3692494     -0.25930120
# L9  -0.2145480       -0.2977739 -0.27233310      -0.1800915     -0.22182692

head(Strader_MEs)
# [1] adult_meso adult_meso adult_meso planu_meso planu_meso planu_meso
#  [7] adult_shal adult_shal adult_shal planu_shal planu_shal planu_shal
# Levels: adult_meso adult_shal planu_meso planu_shal

# Calculate 10 over-arching expression patterns using mean eigengene for each module in a cluster
# Create a column to the Strader_MEs data frame containing age-depth groups
Strader_MEs$group <- c("adult_meso", "adult_meso","adult_meso","planu_meso", "planu_meso","planu_meso","adult_shal","adult_shal","adult_shal", "planu_shal","planu_shal","planu_shal")

C1_Strader_MEs <- select(Strader_MEs, MEskyblue2:MEmediumpurple4)
C1_Strader_MEs$Mean <- rowMeans(C1_Strader_MEs)
C2_Strader_MEs <- select(Strader_MEs, MEdarkolivegreen2:MEthistle)
C2_Strader_MEs$Mean <- rowMeans(C2_Strader_MEs)
C3_Strader_MEs <- select(Strader_MEs, MEskyblue4:MEmediumpurple3)
C3_Strader_MEs$Mean <- rowMeans(C3_Strader_MEs)
C4_Strader_MEs <- select(Strader_MEs, MEantiquewhite2:MEdarkslateblue)
C4_Strader_MEs$Mean <- rowMeans(C4_Strader_MEs)
C5_Strader_MEs <- select(Strader_MEs, MEthistle3:MEturquoise)
C5_Strader_MEs$Mean <- rowMeans(C5_Strader_MEs)
C6_Strader_MEs <- select(Strader_MEs, MEhoneydew:MEorangered4)
C6_Strader_MEs$Mean <- rowMeans(C6_Strader_MEs)
C7_Strader_MEs <- select(Strader_MEs, MEindianred3:MEpaleturquoise)
C7_Strader_MEs$Mean <- rowMeans(C7_Strader_MEs)
C8_Strader_MEs <- select(Strader_MEs, MEmaroon:MElightcyan)
C8_Strader_MEs$Mean <- rowMeans(C8_Strader_MEs)
C9_Strader_MEs <- select(Strader_MEs, MEfloralwhite:MEplum3)
C9_Strader_MEs$Mean <- rowMeans(C9_Strader_MEs)
C10_Strader_MEs <- select(Strader_MEs, MEdarkviolet:MEorangered3)
C10_Strader_MEs$Mean <- rowMeans(C10_Strader_MEs)

expressionProfile_data <- as.data.frame(cbind(group = Strader_MEs$group, cluster1= C1_Strader_MEs$Mean, cluster2 = C2_Strader_MEs$Mean, 
                          cluster3 = C3_Strader_MEs$Mean, cluster4 = C4_Strader_MEs$Mean,cluster5 = C5_Strader_MEs$Mean,cluster6 = C6_Strader_MEs$Mean,
                          cluster7 = C7_Strader_MEs$Mean,cluster8 = C8_Strader_MEs$Mean,cluster9 = C9_Strader_MEs$Mean,cluster10 = C10_Strader_MEs$Mean))

# head(expressionProfile_data)
# [1] group            cluster1           cluster2            cluster3
# 1 adult_meso   0.585881229174129  0.630360718407027   0.539039190845934
# 2 adult_meso -0.0583403960509763 -0.382486696072669   0.241216052692213
# 3 adult_meso  0.0132438841407502 -0.262214789143298   0.196104735206657
# 4 adult_shal  -0.209362218494154 0.0401859629678783 -0.0781554139884685
# 5 adult_shal  -0.250461149130286 0.0351985753526948 -0.0950156522179035
# 6 adult_shal  -0.237314687600277 0.0485129182101299 -0.0695585278294926



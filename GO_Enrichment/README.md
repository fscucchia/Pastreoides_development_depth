## GO enrichment analysis for both the DE and WGCNA data of _P. astreoides_ planulae and adult samples from shallow and mesophotic reefs ######

### GO Enrichment - Differential Expression Data 
```
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

#DESeq2 result
DEG.res <- read.csv("DEGs_adult_cluster.csv")[,-1] #this is the output of the DE script https://github.com/fscucchia/Pastreoides_development_depth/tree/main/DE

#Load Past annotations
annot_final <- read.csv("Past_annot.csv", header = TRUE, sep = ",")[,-1]
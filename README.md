# Pastreoides_development_depth

This electronic notebook provides the scripts employed to analyze gene expression dynamics across developmental stages (planulae and adults) of the coral _Porites astreoides_ inhabiting shallow (10 m) and mesophotic (45 m) reefs in Bermuda.


### RNA-Seq reads quality filtering and mapping

**[Quality filtering and mapping](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/Filtering_and_Mapping)** - RNA-Seq reads processing included adapter trimming using Cutadapt v2.6 ([Martin, 2011](https://doi.org/10.14806/ej.17.1.200)) and quality filtering using Trimmomatic v0.39 ([Bolger et al., 2014](https://doi.org/10.1093/bioinformatics/btu170)). Reads were aligned to the [_P. astreoides_ genome assembly](https://osf.io/ed8xu/) using HISAT2 v2.2.1 ([Kim et al., 2019](https://www.nature.com/articles/s41587-019-0201-4)). Transcripts assembly and quantification were performed using Stringtie v2.2.5 ([Pertea et al., 2015](https://www.nature.com/articles/nbt.3122)).

### Algal symbiont species identification

**[Coral host and symbiont species identification](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/Species_Identification)** - High quality reads were blasted using Diamond v2.0.11 ([Buchfink et al., 2021](https://www.nature.com/articles/s41592-021-01101-x)) against the [NCBI](https://www.ncbi.nlm.nih.gov/), [Reefgenomics](http://reefgenomics.org/), [Marinegenomics](https://marinegenomics.oist.jp/gallery) and [UQ eSpace](https://espace.library.uq.edu.au/view/UQ:f1b3a11) proteomes databases of Symbiodiniaceae species _Symbiodinium microadraticum_, _Symbiodinium tridacnidorum_, _Symbiodinium necroappetens_, _Symbiodinium natans_, _Symbiodinium linuacheae_, _Cladocopium goreaui_, _Cladocopium_ C15, _Fugacium kawagutii_,  and _Durusdinium trenchii_  (formerly _Symbiodinium_ spp. clades A, C, F, and D ([LaJeunesse et al., 2018](https://doi.org/10.1016/j.cub.2018.07.008))).

### SNPs characterization

**[Coral host SNPs characterization](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/SNPs)** - Single nucleotide polymorphisms (SNPs) analysis was conducted using the Genome Analysis Toolkit framework (GATK, v4.2.0; ([McKenna et al., 2010](https://doi.org/10.1101/gr.107524.110))) following the recommended RNA-Seq SNPs practice of the Broad Institute ([(Auwera et al. 2013)](https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/0471250953.bi1110s43)), with necessary adjustments for genotype calling in non-model organisms where variants sites are not known beforehand. HISAT-aligned reads were sorted and marked for duplicates, variant calling was performed with the GATK HaplotypeCaller tool ([McKenna et al., 2010](https://doi.org/10.1101/gr.107524.110)) and genotypes were then jointly called using the GATK GenotypeGVCFs tool. The GATK SelectVariants and VariantFiltration tools were used to filter the joined variant-calling matrices for quality by depth. Filtering for linkage disequilibrium was carried out using PLINK (v2.0, ([Purcell et al. 2007](https://www.cell.com/ajhg/fulltext/S0002-9297(07)61352-4)).
To assess genetic differentiation among age-depth groups, the fixation index (_Fst_)([Weir & Cockerham, 1984](https://doi.org/10.1111/j.1558-5646.1984.tb05657.x)) was estimated using the R package [HIERFSTAT](https://cran.r-project.org/web/packages/hierfstat/index.html) v0.5.10. 

### Differential expression

**[Coral host differential expression](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/DE)** - DE analysis was conducted using Bioconductor DEseq2 v1.26.0 ([Love et al., 2014](https://doi.org/10.1186/s13059-014-0550-8)) in the R environment (v3.6.3) by analyzing planulae and adult samples considering a single factor (depth) with two levels (shallow, mesophotic).

### Weighted correlation network analysis

**[Coral host WGCNA](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/WGCNA)** - WGCNA analysis ([Langfelder and Horvath 2008](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559)) was performed using the R package WGCNA (v1.70.3), with soft thresholding power and adjacency of type “signed”.  

### Gene ontology enrichment

**[Coral host gene ontology enrichment analysis](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/GO_Enrichment)** - GO annotation of the _P. astreoides_ genome was retrieved from the [Past_Genome Project](https://osf.io/ed8xu/). GO enrichment analysis was performed for both the DE and WGCNA data using the package Goseq (v1.42.0; [Young et al. 2010](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-2-r14)) in the R environment.  

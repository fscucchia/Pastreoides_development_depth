# Pastreoides_development_depth

This electronic notebook provides the scripts employed to analyze gene expression dynamics across developmental stages (larvae and adults) of the coral _Porites astreoides_ inhabiting shallow (10 m) and mesophotic (45 m) reefs in Bermuda.

In progress... 😁👩🏻‍💻


### RNA-Seq reads quality filtering and mapping

**[Quality filtering and mapping](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/Filtering_and_Mapping)** - Details the processing and analyses of the _P. astreoides_ transcriptome sequencing data. RNA-Seq reads processing included adapter trimming using Cutadapt v1.15 ([Martin, 2011](https://doi.org/10.14806/ej.17.1.200)) and quality filtering using Trimmomatic v0.3 ([Bolger et al., 2014](https://doi.org/10.1093/bioinformatics/btu170)). Reads were aligned to the host genome assembly using HISAT2 ([Kim et al., 2019](https://www.nature.com/articles/s41587-019-0201-4)). Transcripts assembly and quantification were performed using Stringtie ([Pertea et al., 2015](https://www.nature.com/articles/nbt.3122)).

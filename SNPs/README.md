## RNAseq short variant (SNPs) analysis of _P. astreoides_ samples ######

This script is based on the pipeline for short variant analysis employed for the [Coral_Stress_Phenome](https://github.com/hputnam/Coral_Stress_Phenome/tree/main/Genotype_Analysis/Pocillopora_acuta_PacBio_Assembly/RNAseq_short_variant_analysis) project, with some modifications and adjustments.

Additionally, the pipeline by Dmytro Kryvokhyzha for [genotype calling in a non-model organism](https://evodify.com/gatk-in-non-model-organism/) was checked for comparison.

---

**Installing programs**

1) GATK 
Followed the directions of the [Broad Institute](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4):                                  
- Download [gatk v4.2.0.0](https://github.com/broadinstitute/gatk/releases/tag/4.2.0.0)
```
cd /data/home/mass/fscucchia/programs
unzip gatk-4.2.0.0.zip
```
Added the path to .bashrc
export PATH="/data/home/mass/fscucchia/programs/gatk-4.2.0.0/:$PATH"

---


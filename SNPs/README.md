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
```
export PATH="/data/home/mass/fscucchia/programs/gatk-4.2.0.0/:$PATH"
```
2) rgsam

- Download from the [github repo](https://github.com/djhshih/rgsam)
```
cd /data/home/mass/fscucchia/programs
unzip rgsam-master.zip
cd /data/home/mass/fscucchia/programs/rgsam-master
make 
make install
```
3) PLINK

- Download PLINK from [here](https://www.cog-genomics.org/plink2)

- Added to bashrc
```
export PATH="/data/home/mass/fscucchia/programs/plink2:$PATH"
```
---

### 01- FastqToSam + collect RG + MergeBamAlignment 

Convert paired-fastq to BAM file (sorted by read name), add read group info (RG) to aligned reads + run the MergeBamAlignment command, which also filters the alinged read (e.g. removes secondary alignments), 

#### FastqToSam
Run [`FastqToSam_RUN.sh`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/SNPs/FastqToSam_RUN.sh), which calls for the script [`FastqToSam.sh`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/SNPs/FastqToSam.sh). 
_Took 13 hours for 12 samples._

#### Collect RG
Run [`collectRG_rgsam_RUN.sh`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/SNPs/collectRG_rgsam_RUN.sh), which calls for the script [`collectRG_rgsam.sh`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/SNPs/collectRG_rgsam.sh). I tried to run all samples in array, but for some reason it did not work.
So I had to run each sample individually.
I ran it again changing R1 with R3 in the `collectRG_rgsam.sh` script. 
_Took 10 hours for 12 samples._

#### Prepare your reference file
The GATK uses two files to access and safety check access to the reference files: a .dict dictionary of the contig names and sizes, and a .fai fasta index file to allow efficient random access to the reference bases. You have to generate these files in order to be able to use a fasta file as reference.

- I used CreateSequenceDictionary.jar from Picard to create a .dict file from a fasta file. This produces a SAM-style header file describing the contents of the fasta file.
Run script [`gatk_CreateSequenceDictionary.sh`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/SNPs/gatk_CreateSequenceDictionary.sh) argument 1.

- I used the faidx command in samtools to prepare the fasta index file. This file describes byte offsets in the fasta file for each contig, allowing to compute exactly where a particular reference base at contig:pos is in the fasta file. This produces a text file with one record per line for each of the fasta contigs. Each record is of the: contig, size, location, basesPerLine, bytesPerLine.
Run script [`gatk_CreateSequenceDictionary.sh`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/SNPs/gatk_CreateSequenceDictionary.sh) argument 2.

#### Merge unalinged bam file 
Merge unalinged bam file (now with read group info) with aligned bam file (read group info from unalinged bam is transfered to aligned bam).
Run [`MergeBamAlignment_RUN.sh`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/SNPs/MergeBamAlignment_RUN.sh), which calls for the script [`MergeBamAlignment.sh`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/SNPs/MergeBamAlignment.sh). I tried to run all samples in array, but it did not work again, like the script above. So I run each sample individually.
_Took 13 hours for 12 samples._

### 02- MarkDuplicates
Potential PCR duplicates need to be marked with Picard Tools.

- Merge read groups belonging to the same sample into a single BAM file. Run script [`MarkDuplicates.sh`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/SNPs/MarkDuplicates.sh).
_Took 5 hours._

### 03- SplitNCigarReads
The ‘CIGAR’ (Compact Idiosyncratic Gapped Alignment Report) string is how the SAM/BAM format represents spliced alignments. Understanding the CIGAR string will help you understand how your query sequence aligns to the reference genome.  
Run script [`SplitNCigarReads.sh`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/SNPs/SplitNCigarReads.sh). _Took 2 days_.  
This will split reads that contain Ns in their cigar string (e.g. spanning splicing events in RNAseq data), it identifies all N cigar elements and creates k+1 new reads (where k is the number of N cigar elements). 
This is to distinguish between deletions in exons and large skips due to introns. For mRNA-to-genome alignment, an N operation represents an intron. 

---
**Steps 4 and 5 are part of the GATK "Best Practices" guide but can't really be undertaken with non-model genomes**. These steps in fact require "known sites", i.e. sites where we know beforehand that SNPs occure. This info is not avaliable for non-model systems. Since sites not in this list are considered putative errors that need to be corrected, these steps have to be skipped. 

---

### 06- HaplotypeCaller
Run script [`HaplotypeCaller.sh`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/SNPs/HaplotypeCaller.sh). _Took 21 hours_.

This assumes:  
--sample-ploidy 2 (default)  
--heterozygosity 0.001 (deafult; dont have prior info to update this with)

### 07- Combine *.g.vcf.gz files and call genotypes
Run scripts [`CombineGVCFs.sh`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/SNPs/CombineGVCFs.sh) and [`GenotypeGVCFs.sh`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/SNPs/GenotypeGVCFs.sh). _Together, they took 1.5 hours._

### 08- Select SNPs and Indels
Run script [`Variant_Filtration.sh`](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/SNPs/Variant_Filtration.sh) argument 1.

##### Make diagnostic plots for Variants Scores: 1st-pass filtering
- Extract Variant Quality Scores and Plot  
Run script [`Variant_Filtration.sh`](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/SNPs/Variant_Filtration.sh) argument 2.  
- Make diagnostic plots for Variants Scores  
Run the [`Diagnostic plots for Variants Scores.r`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/SNPs/Diagnostic%20plots%20for%20Variants%20Scores.r) in R.

### 09- Apply Variant filtering
Run script [`Variant_Filtration.sh`](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/SNPs/Variant_Filtration.sh) argument 4, which calls for the script [`Variant_filt.sh`](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/SNPs/Variant_filt.sh)(needs to be in the output Variant_Filtration directory). Here I put all the values of the filtering parameters (SNP_QD_MIN=2.00, SNP_MQ_MIN=50.00, etc.) that I could see in the plots from the previous step. 

<details>
<summary>Troubleshooting tips!</summary>
<br>
NOTE: --filter-expression "QD < $INDEL_QD_MIN" and the other parameters need to be run with the $ sign, just the number does not work.
</details>

##### Check number PASSED after the first filtering
```
cd /data/home/mass/fscucchia/Bermuda/output/Variant_Filtration
OUT="GVCFall_morefilter"
zcat "${OUT}_SNPs.vcf.gz" | grep -v '^#' | wc -l 
#949411
zcat "${OUT}_SNPs_VarScores_filter.vcf.gz" | grep 'PASS' | wc -l
#457216 (48% remain after first-pass filtering) 
zcat "${OUT}_INDELs.vcf.gz" | grep -v '^#' | wc -l
#230546
zcat "${OUT}_INDELs_VarScores_filter.vcf.gz" | grep 'PASS' | wc -l
#157252 (68% remain after first-pass filtering)
```
#### Extract only variants that PASSED filtering
```
zcat "${OUT}_SNPs_VarScores_filter.vcf.gz" | grep -E '^#|PASS' > "${OUT}_SNPs_VarScores_filterPASSED.vcf"
gatk IndexFeatureFile --input "${OUT}_SNPs_VarScores_filterPASSED.vcf" 1> "${OUT}_SNPs_VarScores_filterPASSED.vcf.log" 2>&1
zcat "${OUT}_INDELs_VarScores_filter.vcf.gz" | grep -E '^#|PASS' > "${OUT}_INDELs_VarScores_filterPASSED.vcf"
gatk IndexFeatureFile --input "${OUT}_INDELs_VarScores_filterPASSED.vcf" 1> "${OUT}_INDELs_VarScores_filterPASSED.vcf.log" 2>&1
```
#### Check filtering worked
Show that no variants are left below our threasholds.
Run script [`Variant_Filtration.sh`](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/SNPs/Variant_Filtration.sh) argument 5, which also calls for the R scripts [`check_filtering.R`](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/SNPs/check_filtering.R) and
[`plot_variants_scores2.R`](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/SNPs/plot_variants_scores2.R) (in the output Variant_Filtration directory).

##### 2nd-pass filtering, filter genotypes
When all low confidence variant sites are removed, filter VCF files for genotype quality.
Run script [`Variant_Filtration.sh`](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/SNPs/Variant_Filtration.sh) argument 6.

##### Extract and plot DP info for each sample 
From all samples before previous filtering. Run the following command:
```
for ((i=3; i<=25; i +=2)); do cut -f $i,$((i+1)) GVCFall.DP.table | awk '$1 != "./." {print $2}' > $i.DP; done
```
where 3-25 is odd numbers for 12 samples. This numbering is required because every sample is represented by two columns (GT, DP). 
Split it by samples and keep only positions that have been genotyped (!="./.").

##### Visualize the extracted information, plot DP distribution and define cut-off:
Run the [`plot_DP_distribution.R`](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/SNPs/plot_DP_distribution.R) in R.

View the content of GVCFall.DP.percentiles.txt and select the acceptable cut-off. Discard the genotypes below the 5th percentile and above the 99th percentile.

##### Apply DP filtering
Run script [`Variant_Filtration.sh`](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/SNPs/Variant_Filtration.sh) argument 7, which also calls the R script [`Variant_filt_DP.sh`](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/SNPs/Variant_filt_DP.sh) (in the Variant_Filtration directory).

##### Set filtered sites to no call
Run script [`Variant_Filtration.sh`](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/SNPs/Variant_Filtration.sh) argument 8.

##### 3nd-pass filtering
Run script [`Variant_Filtration.sh`](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/SNPs/Variant_Filtration.sh) argument 10.

### 10- VCF to Table 
Run script [`Variant_Filtration.sh`](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/SNPs/Variant_Filtration.sh) argument 11. 
Got my output files (genotype tables) and vcf files ready for subsequent analyses.

### 11- Filtering for linkage disequilibrium 
- Run script [`genotype_plink.sh`](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/SNPs/genotype_plink.sh) argument 1, which manually sets the variant ID in the final vcf file (it replaces the dot in the ID column with a manually set ID). This is needed for argument 4 to work.
- Run script [`genotype_plink.sh`](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/SNPs/genotype_plink.sh) argument 2, which converts the vcf file into a bed file.

<details>
<summary>Troubleshooting tips!</summary>
<br>
I got the error "Error: _____ cannot contain multiallelic variants" while trying to make the bed file, so I added `--max-alleles 2` in the command.
I also got the error "Error: Multiple instances of '_' in sample ID", so I added `--double-id` in the command.
</details>

- Run script [`genotype_plink.sh`](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/SNPs/genotype_plink.sh) argument 3, which performs the pruning of SNPs that are strongly genetically linked. 
- Run script [`genotype_plink.sh`](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/SNPs/genotype_plink.sh) argument 4, which extracts the pruned SNPs data from the vcf file. 

### 12- Fst analysis 
- Run the R script [`Fst_analysis.R`](https://github.com/fscucchia/Pastreoides_development_depth/tree/main/SNPs/Fst_analysis.R).

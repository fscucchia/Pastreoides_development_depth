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
Run script [`HaplotypeCaller`](). _Took almost 2 days_.
This assumes:  
--sample-ploidy 2 (default)  
--heterozygosity 0.001 (deafult; dont have prior info to update this with)

### 07- Combine *.g.vcf.gz files and call genotypes
Run scripts [`CombineGVCFs.sh`]() and [`GenotypeGVCFs.sh`](). _Together, they took 1.5 hours._

### 08- Select SNPs and Indels
Run script [`Variant_Filtration.sh`]() argument 1.

##### Make diagnostic plots for Variants Scores: 1st-pass filtering
- Extract Variant Quality Scores and Plot  
Run script [`Variant_Filtration.sh`]() argument 2.  
- Make diagnostic plots for Variants Scores  
Run the [`Diagnostic plots for Variants Scores.r`]() in R.

### 09- Apply Variant filtering
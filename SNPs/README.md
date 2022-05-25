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
Run [`FastqToSam_RUN.sh`](), which calls for the script [`FastqToSam.sh`](). 
Took 13 hours for 12 samples.

#### Collect RG
Run [`collectRG_rgsam_RUN.sh`](), which calls for the script [`collectRG_rgsam.sh`](). I tried to run all samples in array, but for some reason it did not work.
So I had to run each sample individually.
I ran it again changing R1 with R3 in the `collectRG_rgsam.sh` script. 
Took 10 hours for 12 samples.

#### Prepare your reference file
The GATK uses two files to access and safety check access to the reference files: a .dict dictionary of the contig names and sizes, and a .fai fasta index file to allow efficient random access to the reference bases. You have to generate these files in order to be able to use a fasta file as reference.

- I used CreateSequenceDictionary.jar from Picard to create a .dict file from a fasta file. This produces a SAM-style header file describing the contents of the fasta file.
Run script [`gatk_CreateSequenceDictionary.sh`]() argument 1.

- I used the faidx command in samtools to prepare the fasta index file. This file describes byte offsets in the fasta file for each contig, allowing to compute exactly where a particular reference base at contig:pos is in the fasta file. This produces a text file with one record per line for each of the fasta contigs. Each record is of the: contig, size, location, basesPerLine, bytesPerLine.
Run script [`gatk_CreateSequenceDictionary.sh`]() argument 2.



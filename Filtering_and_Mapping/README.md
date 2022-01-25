
The following document contains the bioinformatic pipeline used for cleaning, aligning and assembling _P. astreoides_ raw RNA sequences.

---

**Tools used**  

Quality check: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [MultiQC](https://multiqc.info/)  
Quality trimming: [Cutadapt](https://cutadapt.readthedocs.io/en/stable/), [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)  
Alignment to the reference genome: [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml)  
Preparation of alignment for the assembly: [SAMtools](http://www.htslib.org/doc/samtools.html)  
Transcripts assembly and quantification: [StringTie](https://ccb.jhu.edu/software/stringtie/) 

---

### Concatenate reads from different lines 

[Design table]() with samples list.

#### Concatenate forward R1

This is done with the ```cat``` linux command.

```
cat D21_S1_L001_R1_001.fastq.gz D21_S1_L002_R1_001.fastq.gz > /Home/Bermuda/output/concat/D21_R1_concat.fastq.gz
cat D22_S2_L001_R1_001.fastq.gz D22_S2_L002_R1_001.fastq.gz > /Home/Bermuda/output/concat/D22_R1_concat.fastq.gz
cat L7_S4_L001_R1_001.fastq.gz L7_S4_L002_R1_001.fastq.gz > /Home/Bermuda/output/concat/L7_R1_concat.fastq.gz
cat L8_S5_L001_R1_001.fastq.gz L8_S5_L002_R1_001.fastq.gz > /Home/Bermuda/output/concat/L8_R1_concat.fastq.gz
cat L9_S6_L001_R1_001.fastq.gz L9_S6_L002_R1_001.fastq.gz > /Home/Bermuda/output/concat/L9_R1_concat.fastq.gz
cat R13B_S11_L001_R1_001.fastq.gz R13B_S11_L002_R1_001.fastq.gz > /Home/Bermuda/output/concat/R13B_R1_concat.fastq.gz
cat R15B_S12_L001_R1_001.fastq.gz R15B_S12_L002_R1_001.fastq.gz > /Home/Bermuda/output/concat/R15B_R1_concat.fastq.gz
cat R33_S10_L001_R1_001.fastq.gz R33_S10_L002_R1_001.fastq.gz > /Home/Bermuda/output/concat/R33_R1_concat.fastq.gz
cat S4_S8_L001_R1_001.fastq.gz S4_S8_L002_R1_001.fastq.gz > /Home/Bermuda/output/concat/S4_R1_concat.fastq.gz
cat S5_S9_L001_R1_001.fastq.gz S5_S9_L002_R1_001.fastq.gz > /Home/Bermuda/output/concat/S5_R1_concat.fastq.gz
cat S7_S7_L001_R1_001.fastq.gz S7_S7_L002_R1_001.fastq.gz > /Home/Bermuda/output/concat/S7_R1_concat.fastq.gz
```

#### Concatenate reverse R3

```
cat D20_S3_L001_R3_001.fastq.gz D20_S3_L002_R3_001.fastq.gz > /Home/Bermuda/output/concat/D20_R3_concat.fastq.gz
cat D21_S1_L001_R3_001.fastq.gz D21_S1_L002_R3_001.fastq.gz > /Home/Bermuda/output/concat/D21_R3_concat.fastq.gz
cat D22_S2_L001_R3_001.fastq.gz D22_S2_L002_R3_001.fastq.gz > /Home/Bermuda/output/concat/D22_R3_concat.fastq.gz
cat L7_S4_L001_R3_001.fastq.gz L7_S4_L002_R3_001.fastq.gz > /Home/Bermuda/output/concat/L7_R3_concat.fastq.gz
cat L8_S5_L001_R3_001.fastq.gz L8_S5_L002_R3_001.fastq.gz > /Home/Bermuda/output/concat/L8_R3_concat.fastq.gz
cat L9_S6_L001_R3_001.fastq.gz L9_S6_L002_R3_001.fastq.gz > /Home/Bermuda/output/concat/L9_R3_concat.fastq.gz
cat R13B_S11_L001_R3_001.fastq.gz R13B_S11_L002_R3_001.fastq.gz > /Home/Bermuda/output/concat/R13B_R3_concat.fastq.gz
cat R15B_S12_L001_R3_001.fastq.gz R15B_S12_L002_R3_001.fastq.gz > /Home/Bermuda/output/concat/R15B_R3_concat.fastq.gz
cat R33_S10_L001_R3_001.fastq.gz R33_S10_L002_R3_001.fastq.gz > /Home/Bermuda/output/concat/R33_R3_concat.fastq.gz
cat S4_S8_L001_R3_001.fastq.gz S4_S8_L002_R3_001.fastq.gz > /Home/Bermuda/output/concat/S4_R3_concat.fastq.gz
cat S5_S9_L001_R3_001.fastq.gz S5_S9_L002_R3_001.fastq.gz > /Home/Bermuda/output/concat/S5_R3_concat.fastq.gz
cat S7_S7_L001_R3_001.fastq.gz S7_S7_L002_R3_001.fastq.gz > /Home/Bermuda/output/concat/S7_R3_concat.fastq.gz
```

#### Create and activate a new conda environment

```
conda create -n newrnapipeline
conda activate newrnapipeline
```

#### Install all necessary programs within your new conda environment

```
conda install -c bioconda hisat2
conda install -c bioconda samtools
conda install -c bioconda multiqc  
conda install -c bioconda/label/broken trimmomatic
conda install -c bioconda cutadapt
conda install -c bioconda fastqc 

wget -c http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.1.5.Linux_x86_64.tar.gz
tar xzf stringtie-2.1.5.Linux_x86_64.tar.gz
```

### Quality filtering 

Run [`fastq-filter-PE_Conda_1_RUN1.sh`]() that performs quality filtering. Argument 1 uses FastQC to perform the initial quality check of raw reads. Argument 2 calls for the script [`fastq-filter_Conda_job_1.sh`](), which contains all the cutadapt and trimmomatics commands, and removes the adapters using the file [`adapters4d.fa`](). Argument 3 uses FastQC to perform the quality check of the filtered reads. Argument 4 compiles the MultiQC report. 

### Alignment of the clean reads to _P. astreoides_ reference genome 



### Assemble aligned reads and quantify transcripts 



### Assess the performance of the assembly 


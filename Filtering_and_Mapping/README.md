
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

[Design table](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/Metadata/design_PE.csv) with samples list.

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

Run [`fastq-filter-PE_Conda_1_RUN1.sh`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/Filtering_and_Mapping/fastq-filter-PE_Conda_1_RUN1.sh) that performs quality filtering. Argument 1 uses FastQC to perform the initial quality check of raw reads. Argument 2 calls for the script [`fastq-filter_Conda_job_1.sh`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/Filtering_and_Mapping/fastq-filter_Conda_job_1.sh), which contains all the cutadapt and trimmomatics commands, and removes the adapters using the file [`adapters4d.fa`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/Filtering_and_Mapping/adapters4d.fa). Argument 3 uses FastQC to perform the quality check of the filtered reads.
Run [`MultiQC.sh`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/Filtering_and_Mapping/MultiQC.sh) for multiQC.

### Alignment of the clean reads to _P. astreoides_ reference genome 

Create a new directory for Hisat
```
mkdir HISAT
```
Run the script [`StyHISAT2_4.sh`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/Filtering_and_Mapping/StyHISAT2_4.sh) that does the following:
1) Symbolically links the filtered fastq files to the HISAT directory.
```
/data/home/mass/fscucchia/Bermuda/output/filtered/*gz.filtered ./
```
2) Index the _P. astreoides_ reference genome in the reference directory and performs the alignment using the script [`StyHISAT_withSummary.sh`](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/Filtering_and_Mapping/StyHISAT_withSummary.sh). 
```
hisat2-build -f /data/home/mass/fscucchia/databases/Pastreoides_genome_KW/past_filtered_assembly.fasta ./Past_ref
```
```
#StyHISAT_withSummary.sh:

   F="/data/home/mass/fscucchia/Bermuda/output/filtered"
   # Aligning paired end reads
   # The R1 in array1 is changed to R3 in the for loop. SAM files are of both forward and reverse reads
   array1=($(ls $F/*_R1_concat.fastq.gz.filtered))

   # Bam files are created and sorted, since Stringtie takes sorted file as input
   # The sam file is removed at the end since it is not needed anymore
   # The command --summary-file ${i}.txt reates a summary file per sample, which can be used by multiqc

   for i in ${array1[@]}; do
        hisat2 -p 8 --new-summary --rf --dta -q -x Past_ref -1 ${i} -2 $(echo ${i}|sed s/_R1/_R3/) -S ${i}.sam --summary-file ${i}.txt 
        samtools sort -@ 8 -o ${i}.bam ${i}.sam
    		echo "${i}_bam"
        rm ${i}.sam
        echo "HISAT2 PE ${i}" $(date)
   done
```
3) Runs multiqc at the end of mapping if the `--summary-file` option was enabled.

<details>
<summary>Troubleshooting tips!</summary>
<br>
After the alignment of the first sample, I got this message "samtools: error while loading shared libraries: libcrypto.so.1.0.0: cannot open shared object file: No such file or directory"
So I ran these 3 commands: $ conda config --add channels bioconda
                           $ conda config --add channels conda-forge
                           $ conda install samtools==1.11
</details>

- Alignment with Hisat took 7 hours (12 samples PE).

### Assemble aligned reads and quantify transcripts 



### Assess the performance of the assembly 


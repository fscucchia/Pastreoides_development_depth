
The following document contains the bioinformatic pipeline used for cleaning, aligning and assembling our raw RNA sequences.

---

**Tools used**  
Quality check: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [MultiQC](https://multiqc.info/)  
Quality trimming: [Fastp](https://github.com/OpenGene/fastp)  
Alignment to reference genome: [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml)  
Preparation of alignment for assembly: [SAMtools](http://www.htslib.org/doc/samtools.html)  
Transcript assembly and quantification: [StringTie](https://ccb.jhu.edu/software/stringtie/) 

---

### Concatenate reads from different lines 

# concatenate forward R1
cd /data/home/mass/tmass/rawdata20210721/210726_A00929_0394_BHG3TFDRXY/fastq/D20
cat D21_S1_L001_R1_001.fastq.gz D21_S1_L002_R1_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/D21_R1_concat.fastq.gz
cd /data/home/mass/tmass/rawdata20210721/210726_A00929_0394_BHG3TFDRXY/fastq/D22
cat D22_S2_L001_R1_001.fastq.gz D22_S2_L002_R1_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/D22_R1_concat.fastq.gz
cat L7_S4_L001_R1_001.fastq.gz L7_S4_L002_R1_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/L7_R1_concat.fastq.gz
cat L8_S5_L001_R1_001.fastq.gz L8_S5_L002_R1_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/L8_R1_concat.fastq.gz
cat L9_S6_L001_R1_001.fastq.gz L9_S6_L002_R1_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/L9_R1_concat.fastq.gz
cat R13B_S11_L001_R1_001.fastq.gz R13B_S11_L002_R1_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/R13B_R1_concat.fastq.gz
cat R15B_S12_L001_R1_001.fastq.gz R15B_S12_L002_R1_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/R15B_R1_concat.fastq.gz
cat R33_S10_L001_R1_001.fastq.gz R33_S10_L002_R1_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/R33_R1_concat.fastq.gz
cat S4_S8_L001_R1_001.fastq.gz S4_S8_L002_R1_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/S4_R1_concat.fastq.gz
cat S5_S9_L001_R1_001.fastq.gz S5_S9_L002_R1_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/S5_R1_concat.fastq.gz
cat S7_S7_L001_R1_001.fastq.gz S7_S7_L002_R1_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/S7_R1_concat.fastq.gz

# concatenate reverse R3
cd /data/home/mass/tmass/rawdata20210721/210726_A00929_0394_BHG3TFDRXY/fastq/D20
cat D20_S3_L001_R3_001.fastq.gz D20_S3_L002_R3_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/D20_R3_concat.fastq.gz
cat D21_S1_L001_R3_001.fastq.gz D21_S1_L002_R3_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/D21_R3_concat.fastq.gz
cat D22_S2_L001_R3_001.fastq.gz D22_S2_L002_R3_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/D22_R3_concat.fastq.gz
cat L7_S4_L001_R3_001.fastq.gz L7_S4_L002_R3_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/L7_R3_concat.fastq.gz
cat L8_S5_L001_R3_001.fastq.gz L8_S5_L002_R3_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/L8_R3_concat.fastq.gz
cat L9_S6_L001_R3_001.fastq.gz L9_S6_L002_R3_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/L9_R3_concat.fastq.gz
cat R13B_S11_L001_R3_001.fastq.gz R13B_S11_L002_R3_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/R13B_R3_concat.fastq.gz
cat R15B_S12_L001_R3_001.fastq.gz R15B_S12_L002_R3_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/R15B_R3_concat.fastq.gz
cat R33_S10_L001_R3_001.fastq.gz R33_S10_L002_R3_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/R33_R3_concat.fastq.gz
cat S4_S8_L001_R3_001.fastq.gz S4_S8_L002_R3_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/S4_R3_concat.fastq.gz
cat S5_S9_L001_R3_001.fastq.gz S5_S9_L002_R3_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/S5_R3_concat.fastq.gz
cat S7_S7_L001_R3_001.fastq.gz S7_S7_L002_R3_001.fastq.gz > /data/home/mass/fscucchia/Bermuda/output/concat/S7_R3_concat.fastq.gz

###Create and activate a conda environment
conda create -n newrnapipeline
conda activate newrnapipeline

###Install all necessary programs within your conda environment
conda install -c bioconda hisat2
conda install -c bioconda samtools
conda install -c bioconda multiqc # I installed multiqc in a new environment called "multiqc_env" to avoi conflicts with other softwares in "newrnapipeline"
conda install -c bioconda/label/broken trimmomatic
conda install -c bioconda cutadapt
conda install -c bioconda fastqc # I installed multiqc in the environment called "multiqc_env" to avoi conflicts with other softwares in "newrnapipeline"

wget -c http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.1.5.Linux_x86_64.tar.gz
tar xzf stringtie-2.1.5.Linux_x86_64.tar.gz
#cd ~/bin/
#ln -s ~/src/stringtie-2.1.5.Linux_x86_64/stringtie



###### Quality filtering ######

###Create and activate a conda environment
conda create -n qualityfilt python=3   #with python 3.9.6
conda activate qualityfilt #needed to use 'fastq-filter_Conda_job_1.sh'

###Install all necessary programs within your conda environment
conda install -c bioconda cutadapt  #needed to use 'fastq-filter_Conda_job_1.sh'
conda install -c bioconda/label/broken trimmomatic  #needed to use 'fastq-filter_Conda_job_1.sh'

Run 'fastq-filter-PE_Conda_1_RUN1.sh' that performs quality filtering (run starting from argument 0, then 1,2,3 and 5) #outside conda env
Argument 2 calls for 'fastq-filter_Conda_job_1.sh' (this script has all the cutadapt and trimmomatics commands) and removes the adapters using the file 'adapters4d.fa'.



###### Alignment of clean reads to reference genome ######

#run the script StyHISAT2_3.sh (takes 1H and 30 min for 3 samples PE) that does the following:
#1) Create a subdirectory within data for HISAT2 and symbolically link it your clean fastq files.
mkdir hisat2
cd hisat2
ln -s /data/home/mass/fscucchia/pH_exp_1/filtered/*filtered.gz* ./
#2) Index the reference genome in the reference directory and alignment (calling for the script 'StyHISAT2.sh')
cd /data/home/mass/fscucchia/pH_exp_1/ref_genome
hisat2-build -f /data/home/mass/fscucchia/pH_exp_1/ref_genome/smic_stylophoraP_concat.fasta ./Sty_ref

#StyHISAT2.sh:

   F="/data/home/mass/fscucchia/pH_exp_1/filtered"
   #Aligning paired end reads
   #Has the R1 in array1 because the sed in the for loop changes it to an R2. SAM files are of both forward and reverse reads
   array1=($(ls $F/*_R1_001.fastq.gz.filtered.gz))

   # This then makes it into a bam file
   # And then also sorts the bam file because Stringtie takes a sorted file for input
   # And then removes the sam file because I don't need it anymore

   for i in ${array1[@]}; do
        hisat2 -p 8 --rf --dta -q -x Sty_ref -1 ${i} -2 $(echo ${i}|sed s/_R1/_R2/) -S ${i}.sam
        samtools sort -@ 8 -o ${i}.bam ${i}.sam
    		echo "${i}_bam"
        rm ${i}.sam
        echo "HISAT2 PE ${i}" $(date)
   done

#after running the script I got this message "Warning: Unsupported file format", not sure what is the problem here

#after the alignment of the first sample, I got this message "samtools: error while loading shared libraries: libcrypto.so.1.0.0: cannot open shared object file: No such file or directory"
#so I ran these 3 commands: conda config --add channels bioconda
                           #conda config --add channels conda-forge
                           #conda install samtools==1.11

#created a new script StyHISAT2_4.sh to get the multiqc at the end of the mapping. This script calls for StyHISAT_withSummary.sh, which has the command --summary-file ${i}.txt 
#to make a summary file per sample, which will be used by multiqc. StyHISAT_withSummary.sh needs to be in the same directory as the ref_genome to work.



### Assemble aligned reads and quantify transcripts ###
#First, create and enter into StringTie directory. 
mkdir ../stringtie
cd stringtie

#run the script 'StySTRINGTIE.sh' (takes 18 min for 3 samples PE) that does the following:
#1) Create a symbolic link to our reference genome (in .gtf format) inside the stringtie directory
ln -s /data/home/mass/fscucchia/pH_exp_1/ref_genome/smic_stylophoraP_concat.gtf ./
#2) then it runs the script StyStringTie_assembly.sh:

  ##!/bin/bash

  #Specify working directory
  W="/data/home/mass/fscucchia/pH_exp_1/filtered"

  #StringTie reference-guided assembly
  #Has the R1 in array1 because of the naming convention in the former script. However, these BAM files contain both forward and reverse reads.
  array1=($(ls $W/*_R1_001.fastq.gz.filtered.gz.bam))

  for i in ${array1[@]}; do
        stringtie -A gene_abundance/{i}.gene_abund.tab -p 8 --rf -e -G smic_stylophoraP_concat.gtf -o ${i}.gtf ${i}
        mv /ref-guided-gtfs/${i}.gtf /data/home/mass/fscucchia/pH_exp_1/output/stringtie/BAM
        echo "StringTie-assembly-to-ref ${i}" $(date)
  done



### Assess the performance of the assembly ###
#First, install gffcompare within your conda environment
conda activate newrnapipeline
conda install -c bioconda gffcompare

#Then create the .txt file, list_to_merge.txt. This file lists all of the file names to be merged. This file is in the StringTie directory.
#Run the script 'Stringtie_merge_compare_count.sh', which does the following:
#1) merge the assembly-generated GTF files to assess how well the predicted transcripts track to the reference annotation file
stringtie --merge -p 8 -G /data/home/mass/fscucchia/pH_exp_1/ref_genome/smic_stylophoraP_concat.gtf -o ../stringtie_merged.gtf list_to_merge.txt
#2) use the program gffcompare to compare the merged GTF to our reference genome
gffcompare -r /data/home/mass/fscucchia/pH_exp_1/ref_genome/smic_stylophoraP_concat.gtf -o ../compared stringtie_merged.gtf
    #Manually move all of the gffcompare output files to a separate directory "/data/home/mass/fscucchia/pH_exp_1/output/Gffcompare"
#3) Compilation of GTF-files into gene and transcript count matrices
#The StringTie program includes a script, prepDE.py that compiles your assembly files into gene and transcript count matrices. This script requires 
#as input the list of sample names and their full file paths, sample_list.txt. This file will live in StringTie program directory.
./prepDE.py -g ../gene_count_matrix.csv -i ./sample_list.txt
    #Manually move the 'gene_count_matrix.csv' and 'transcript_count_matrix.csv' files to a separate directory "/data/home/mass/fscucchia/pH_exp_1/output/count_matrix"

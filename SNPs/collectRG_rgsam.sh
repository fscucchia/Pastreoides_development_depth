
#!/bin/sh

cd /data/home/mass/fscucchia/Bermuda/output/filtered
samtools view S7_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | rgsam collect -s S7_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam -o rg.txt
samtools view -h S7_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | 
  rgsam tag -r rg.txt |
  samtools view -b - > S7_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam	

cd /data/home/mass/fscucchia/Bermuda/output/filtered
samtools view S4_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | rgsam collect -s S4_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam -o rg.txt
samtools view -h S4_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | 
  rgsam tag -r rg.txt |
  samtools view -b - > S4_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam	

cd /data/home/mass/fscucchia/Bermuda/output/filtered
samtools view S5_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | rgsam collect -s S5_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam -o rg.txt
samtools view -h S5_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | 
  rgsam tag -r rg.txt |
  samtools view -b - > S5_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam	

cd /data/home/mass/fscucchia/Bermuda/output/filtered
samtools view R33_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | rgsam collect -s R33_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam -o rg.txt
samtools view -h R33_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | 
  rgsam tag -r rg.txt |
  samtools view -b - > R33_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam	

cd /data/home/mass/fscucchia/Bermuda/output/filtered
samtools view R15B_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | rgsam collect -s R15B_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam -o rg.txt
samtools view -h R15B_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | 
  rgsam tag -r rg.txt |
  samtools view -b - > R15B_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam	

cd /data/home/mass/fscucchia/Bermuda/output/filtered
samtools view R13B_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | rgsam collect -s R13B_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam -o rg.txt
samtools view -h R13B_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | 
  rgsam tag -r rg.txt |
  samtools view -b - > R13B_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam	

cd /data/home/mass/fscucchia/Bermuda/output/filtered
samtools view L9_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | rgsam collect -s L9_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam -o rg.txt
samtools view -h L9_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | 
  rgsam tag -r rg.txt |
  samtools view -b - > L9_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam	

cd /data/home/mass/fscucchia/Bermuda/output/filtered
samtools view L8_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | rgsam collect -s L8_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam -o rg.txt
samtools view -h L8_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | 
  rgsam tag -r rg.txt |
  samtools view -b - > L8_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam	

cd /data/home/mass/fscucchia/Bermuda/output/filtered
samtools view L7_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | rgsam collect -s L7_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam -o rg.txt
samtools view -h L7_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | 
  rgsam tag -r rg.txt |
  samtools view -b - > L7_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam	

cd /data/home/mass/fscucchia/Bermuda/output/filtered
samtools view D20_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | rgsam collect -s D20_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam -o rg.txt
samtools view -h D20_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | 
  rgsam tag -r rg.txt |
  samtools view -b - > D20_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam	

cd /data/home/mass/fscucchia/Bermuda/output/filtered
samtools view D21_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | rgsam collect -s D21_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam -o rg.txt
samtools view -h D21_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | 
  rgsam tag -r rg.txt |
  samtools view -b - > D21_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam	

cd /data/home/mass/fscucchia/Bermuda/output/filtered
samtools view D22_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | rgsam collect -s D22_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam -o rg.txt
samtools view -h D22_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.bam | 
  rgsam tag -r rg.txt |
  samtools view -b - > D22_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam	




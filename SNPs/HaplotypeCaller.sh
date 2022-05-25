
#!/bin/sh

F="/data/home/mass/fscucchia/Bermuda/output/SplitNCigarReads"
G="/data/home/mass/fscucchia/databases/Pastreoides_genome_Kevin/past_filtered_assembly.fasta"


if [ $1 -eq 1 ]; then
     sbatch --mem=300000 -N1 -n20 --ntasks-per-node=20 --workdir=$F --job-name "HaplotypeCaller" -o "$F/HaplotypeCaller.out" -e "$F/HaplotypeCaller.err" \
     -p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d \
  	 --wrap "cd /data/home/mass/fscucchia/Bermuda/output/SplitNCigarReads
             gatk HaplotypeCaller --reference $G --input S7_R1.SplitNCigarReads.split.bam --output S7.HaplotypeCaller.g.vcf.gz -dont-use-soft-clipped-bases -ERC GVCF
cd /data/home/mass/fscucchia/Bermuda/output/SplitNCigarReads
             gatk HaplotypeCaller --reference $G --input S4_R1.SplitNCigarReads.split.bam --output S4.HaplotypeCaller.g.vcf.gz -dont-use-soft-clipped-bases -ERC GVCF 
cd /data/home/mass/fscucchia/Bermuda/output/SplitNCigarReads
             gatk HaplotypeCaller --reference $G --input S5_R1.SplitNCigarReads.split.bam --output S5.HaplotypeCaller.g.vcf.gz -dont-use-soft-clipped-bases -ERC GVCF 
cd /data/home/mass/fscucchia/Bermuda/output/SplitNCigarReads
             gatk HaplotypeCaller --reference $G --input R33_R1.SplitNCigarReads.split.bam --output R33.HaplotypeCaller.g.vcf.gz -dont-use-soft-clipped-bases -ERC GVCF 
cd /data/home/mass/fscucchia/Bermuda/output/SplitNCigarReads
             gatk HaplotypeCaller --reference $G --input R15B_R1.SplitNCigarReads.split.bam --output R15B.HaplotypeCaller.g.vcf.gz -dont-use-soft-clipped-bases -ERC GVCF 
cd /data/home/mass/fscucchia/Bermuda/output/SplitNCigarReads
             gatk HaplotypeCaller --reference $G --input R13B_R1.SplitNCigarReads.split.bam --output R13B.HaplotypeCaller.g.vcf.gz -dont-use-soft-clipped-bases -ERC GVCF 
cd /data/home/mass/fscucchia/Bermuda/output/SplitNCigarReads
             gatk HaplotypeCaller --reference $G --input L7_R1.SplitNCigarReads.split.bam --output L7.HaplotypeCaller.g.vcf.gz -dont-use-soft-clipped-bases -ERC GVCF 
cd /data/home/mass/fscucchia/Bermuda/output/SplitNCigarReads
             gatk HaplotypeCaller --reference $G --input L8_R1.SplitNCigarReads.split.bam --output L8.HaplotypeCaller.g.vcf.gz -dont-use-soft-clipped-bases -ERC GVCF 
cd /data/home/mass/fscucchia/Bermuda/output/SplitNCigarReads
             gatk HaplotypeCaller --reference $G --input L9_R1.SplitNCigarReads.split.bam --output L9.HaplotypeCaller.g.vcf.gz -dont-use-soft-clipped-bases -ERC GVCF 
cd /data/home/mass/fscucchia/Bermuda/output/SplitNCigarReads
             gatk HaplotypeCaller --reference $G --input D20_R1.SplitNCigarReads.split.bam --output D20.HaplotypeCaller.g.vcf.gz -dont-use-soft-clipped-bases -ERC GVCF 
cd /data/home/mass/fscucchia/Bermuda/output/SplitNCigarReads
             gatk HaplotypeCaller --reference $G --input D21_R1.SplitNCigarReads.split.bam --output D21.HaplotypeCaller.g.vcf.gz -dont-use-soft-clipped-bases -ERC GVCF 
cd /data/home/mass/fscucchia/Bermuda/output/SplitNCigarReads
             gatk HaplotypeCaller --reference $G --input D22_R1.SplitNCigarReads.split.bam --output D22.HaplotypeCaller.g.vcf.gz -dont-use-soft-clipped-bases -ERC GVCF"

fi    

    
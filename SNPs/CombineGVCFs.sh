
#!/bin/sh

F="/data/home/mass/fscucchia/Bermuda/output/SplitNCigarReads"
REF="/data/home/mass/fscucchia/databases/Pastreoides_genome_Kevin/past_filtered_assembly.fasta"

if [ $1 -eq 1 ]; then
     sbatch -N1 -n24 --ntasks-per-node=24 --mem=300000 --workdir=$F --job-name "CombineGVCFs" -o "$F/CombineGVCFs.out" -e "$F/CombineGVCFs.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap "gatk CombineGVCFs -R "${REF}" -V D22.HaplotypeCaller.g.vcf.gz -V D20.HaplotypeCaller.g.vcf.gz -V D21.HaplotypeCaller.g.vcf.gz -V L7.HaplotypeCaller.g.vcf.gz -V L8.HaplotypeCaller.g.vcf.gz -V L9.HaplotypeCaller.g.vcf.gz -V R33.HaplotypeCaller.g.vcf.gz -V R13B.HaplotypeCaller.g.vcf.gz -V R15B.HaplotypeCaller.g.vcf.gz -V S4.HaplotypeCaller.g.vcf.gz -V S5.HaplotypeCaller.g.vcf.gz -V S7.HaplotypeCaller.g.vcf.gz -O cohort.g.vcf.gz"
fi




 

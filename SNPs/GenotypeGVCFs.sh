
#!/bin/sh

F="/data/home/mass/fscucchia/Bermuda/output/SplitNCigarReads"
REF="/data/home/mass/fscucchia/databases/Pastreoides_genome_Kevin/past_filtered_assembly.fasta"


if [ $1 -eq 1 ]; then
     sbatch -N1 -n24 --ntasks-per-node=24 --mem=128000 --workdir=$F --job-name "GenotypeGVCFs" -o "$F/GenotypeGVCFs.out" -e "$F/GenotypeGVCFs.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap "gatk --java-options "-Xmx4g" GenotypeGVCFs --reference "${REF}" --output "cohort_genotypes.vcf.gz" -V "cohort.g.vcf.gz" -stand-call-conf 30 --annotation AS_MappingQualityRankSumTest --annotation AS_ReadPosRankSumTest > GenotypeGVCFs.log 2>&1"
fi



 

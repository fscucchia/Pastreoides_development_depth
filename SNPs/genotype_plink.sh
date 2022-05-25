
#!/bin/sh

P="/data/home/mass/fscucchia/programs/plink_2"
O="/data/home/mass/fscucchia/Bermuda/output/genotype_analysis"

##########################################################################

if [ $1 -eq 1 ]; then
     sbatch -N1 -n24 --ntasks-per-node=24 --mem=128000 --workdir=$P --job-name "bcftools" -o "$P/bcftools.out" -e "$P/bcftools.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	--wrap ". $CONDA; conda activate functional_annotation
             bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' /data/home/mass/fscucchia/programs/plink_2/GVCFall_SNPs_VarScores_filterPASSED_DPfilterNoCall.vcf"

elif [ $1 -eq 2 ]; then
     sbatch -N1 -n24 --ntasks-per-node=24 --mem=128000 --workdir=$P --job-name "MakeBed" -o "$P/MakeBed.out" -e "$P/MakeBed.err" \
     -p hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap "./plink2 --vcf /data/home/mass/fscucchia/programs/plink_2/bcftools.out --make-bed --chr-set 28 no-xy --allow-extra-chr --double-id --max-alleles 2 --out VCF_annotated"

elif [ $1 -eq 3 ]; then
     sbatch -N1 -n24 --ntasks-per-node=24 --mem=128000 --workdir=$P --job-name "plink" -o "$P/plink.out" -e "$P/plink.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap "./plink2 --bfile /data/home/mass/fscucchia/programs/plink_2/VCF_annotated --allow-no-sex --allow-extra-chr --bad-ld --indep-pairwise 25 0.5"


elif [ $1 -eq 4 ]; then
     sbatch -N1 -n24 --ntasks-per-node=24 --mem=128000 --workdir=$P --job-name "plink" -o "$P/plink.out" -e "$P/plink.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap "./plink2 --vcf /data/home/mass/fscucchia/programs/plink_2/bcftools.out --extract plink2.prune.in --allow-extra-chr --make-bed --double-id --out P_astreoides_vcf_pruned"

fi



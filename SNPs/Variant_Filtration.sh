
#!/bin/sh

F="/data/home/mass/fscucchia/Bermuda/output/SplitNCigarReads"
REF="/data/home/mass/fscucchia/databases/Pastreoides_genome_Kevin/past_filtered_assembly.fasta"
OUT="GVCFall"
O="/data/home/mass/fscucchia/Bermuda/output/Variant_Filtration"
####
SNP_QD_MIN=2.00
SNP_MQ_MIN=50.00
SNP_FS_MAX=60.000
SNP_SOR_MAX=4.000

INDEL_QD_MIN=2.00
INDEL_MQ_MIN=45.00
INDEL_FS_MAX=200.000
INDEL_SOR_MAX=10.000

DP_MIN=2
####

if [ $1 -eq 1 ]; then
     sbatch -N1 -n24 --ntasks-per-node=24 --mem=128000 --workdir=$F --job-name "SelectVariants" -o "$O/SelectVariants.out" -e "$O/SelectVariants.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap "gatk SelectVariants --reference "${REF}" --variant "cohort_genotypes.vcf.gz" --output "/data/home/mass/fscucchia/Bermuda/output/Variant_Filtration/${OUT}_SNPs.vcf.gz"   -select-type SNP 1> "/data/home/mass/fscucchia/Bermuda/output/Variant_Filtration/${OUT}_SNPs.vcf.gz.log" 2>&1
             gatk SelectVariants --reference "${REF}" --variant "cohort_genotypes.vcf.gz" --output "/data/home/mass/fscucchia/Bermuda/output/Variant_Filtration/${OUT}_INDELs.vcf.gz" -select-type INDEL 1> "/data/home/mass/fscucchia/Bermuda/output/Variant_Filtration/${OUT}_INDELs.vcf.gz.log" 2>&1
            "

elif [ $1 -eq 2 ]; then
     sbatch -N1 -n20 --ntasks-per-node=20 --mem=128000 --workdir=$O --job-name "VariantsToTable" -o "$O/VariantsToTable.out" -e "$O/VariantsToTable.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap "gatk VariantsToTable --reference "${REF}" --variant "${OUT}_SNPs.vcf.gz"   --output "${OUT}_SNPs.table"   -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum 1> "${OUT}_SNPs.table.log" 2>&1
              gatk VariantsToTable --reference "${REF}" --variant "${OUT}_INDELs.vcf.gz" --output "${OUT}_INDELs.table" -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum 1> "${OUT}_INDELs.table.log" 2>&1
              "

elif [ $1 -eq 3 ]; then
     sbatch -N1 -n20 --ntasks-per-node=20 --mem=128000 --workdir=$O --job-name "variants_scores" -o "$O/variants_scores.out" -e "$O/variants_scores.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap ". $CONDA; conda activate R.settings2
                        Rscript plot_variants_scores.R "${OUT}.variants_scores" "${OUT}_SNPs.table" "${OUT}_INDELs.table" \
	                   $SNP_QD_MIN $SNP_MQ_MIN $SNP_FS_MAX $SNP_SOR_MAX \
	                   $INDEL_QD_MIN $INDEL_MQ_MIN $INDEL_FS_MAX $INDEL_SOR_MAX \
	                   1> "${OUT}.variants_scores.log" 2>&1
                        conda deactivate
                        "

elif [ $1 -eq 4 ]; then
     sbatch -N1 -n20 --ntasks-per-node=20 --mem=128000 --workdir=$O --job-name "VariantFiltration" -o "$O/VariantFiltration.out" -e "$O/VariantFiltration.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap "chmod u+x Variant_filt.sh; ./Variant_filt.sh"

elif [ $1 -eq 5 ]; then
     sbatch -N1 -n20 --ntasks-per-node=20 --mem=128000 --workdir=$O --job-name "VariantsToTable" -o "$O/VariantsToTable.out" -e "$O/VariantsToTable.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap "gatk VariantsToTable --reference "${REF}" --variant "${OUT}_SNPs_VarScores_filterPASSED.vcf"   --output "${OUT}_SNPs_VarScores_filterPASSED.table" \
	         -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum \
	         1> "${OUT}_SNPs_VarScores_filterPASSED.table.log" 2>&1
              gatk VariantsToTable --reference "${REF}" --variant "${OUT}_INDELs_VarScores_filterPASSED.vcf" --output "${OUT}_INDELs_VarScores_filterPASSED.table" \
	         -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum \
	         1> "${OUT}_INDELs_VarScores_filterPASSED.table.log" 2>&1
              . $CONDA; conda activate R.settings2
                        Rscript check_filtering.R "${OUT}_SNPs_VarScores_filterPASSED.table" "${OUT}_INDELs_VarScores_filterPASSED.table" \
                        $SNP_QD_MIN $SNP_MQ_MIN $SNP_FS_MAX $SNP_SOR_MAX \
                        $INDEL_QD_MIN $INDEL_MQ_MIN $INDEL_FS_MAX $INDEL_SOR_MAX \
                    	1> "${OUT}.check_filtering.log" 2>&1
                        Rscript plot_variants_scores2.R "${OUT}.variants_scores_afterFiltering" "${OUT}_SNPs_VarScores_filterPASSED.table" "${OUT}_INDELs_VarScores_filterPASSED.table" \
                        $SNP_QD_MIN $SNP_MQ_MIN $SNP_FS_MAX $SNP_SOR_MAX \
                        $INDEL_QD_MIN $INDEL_MQ_MIN $INDEL_FS_MAX $INDEL_SOR_MAX \
	                   1> "${OUT}.variants_scores_afterFiltering.log" 2>&1
                        conda deactivate"

elif [ $1 -eq 6 ]; then
     sbatch -N1 -n20 --ntasks-per-node=20 --mem=128000 --workdir=$O --job-name "VariantsToTable" -o "$O/VariantsToTable.out" -e "$O/VariantsToTable.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap "gatk VariantsToTable --reference "${REF}" --variant "/data/home/mass/fscucchia/Bermuda/output/SplitNCigarReads/cohort_genotypes.vcf.gz" --output "${OUT}.DP.table" -F CHROM -F POS -GF GT -GF DP 1> "${OUT}.DP.table.log" 2>&1
             "
      
elif [ $1 -eq 7 ]; then
     sbatch -N1 -n20 --ntasks-per-node=20 --mem=128000 --workdir=$O --job-name "Variant_filt_DP" -o "$O/Variant_filt_DP.out" -e "$O/Variant_filt_DP.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap "chmod u+x Variant_filt_DP.sh; ./Variant_filt_DP.sh"

elif [ $1 -eq 8 ]; then
     sbatch -N1 -n20 --ntasks-per-node=20 --mem=128000 --workdir=$O --job-name "SelectVariants" -o "$O/SelectVariants.out" -e "$O/SelectVariants.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap "gatk SelectVariants --reference "${REF}" --variant "${OUT}_SNPs_VarScores_filterPASSED_DPfilter.vcf" --output "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.vcf" --set-filtered-gt-to-nocall \
	1> "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.vcf.log" 2>&1
     gatk SelectVariants --reference "${REF}" --variant "${OUT}_INDELs_VarScores_filterPASSED_DPfilter.vcf" --output "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.vcf" --set-filtered-gt-to-nocall \
	1> "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.vcf.log" 2>&1"

elif [ $1 -eq 9 ]; then
     sbatch -N1 --workdir=$O --job-name "Extract_plot_DP" -o "$O/Extract_plot_DP.out" -e "$O/Extract_plot_DP.err" \
     -p hiveunlim \
  	 --wrap "chmod u+x Extract_plot_DP_info.sh; ./Extract_plot_DP_info.sh"

elif [ $1 -eq 10 ]; then
     sbatch -N1 -n24 --ntasks-per-node=24 --mem=128000 --workdir=$O --job-name "VariantsToTable" -o "$O/VariantsToTable.out" -e "$O/VariantsToTable.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	--wrap "gatk VariantsToTable --reference "${REF}" --variant "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.vcf" --output "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.table" -F CHROM -F POS -GF GT -GF AD -GF DP \
	1> "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.AD.table.log" 2>&1"

elif [ $1 -eq 11 ]; then
     sbatch -N1 -n24 --ntasks-per-node=24 --mem=128000 --workdir=$O --job-name "VariantsToTable" -o "$O/VariantsToTable.out" -e "$O/VariantsToTable.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
  	 --wrap "gatk VariantsToTable --reference "${REF}" --variant "${OUT}_SNPs_VarScores_filterPASSED_DPfilterNoCall.vcf"   --output "${OUT}_SNPs_filterPASSED_final.table"   -F CHROM -F POS -GF GT \
	1> "${OUT}_SNPs_filterPASSED_final.table.log" 2>&1
     gatk VariantsToTable --reference "${REF}" --variant "${OUT}_INDELs_VarScores_filterPASSED_DPfilterNoCall.vcf" --output "${OUT}_INDELs_filterPASSED_final.table" -F CHROM -F POS -GF GT \
	1> "${OUT}_INDELs_filterPASSED_final.table.log" 2>&1"

     
fi



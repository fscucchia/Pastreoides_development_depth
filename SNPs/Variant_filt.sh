#!/bin/sh

F="/data/home/mass/fscucchia/Bermuda/output/SplitNCigarReads"
REF="/data/home/mass/fscucchia/databases/Pastreoides_genome_Kevin/past_filtered_assembly.fasta"
OUT="GVCFall"
O="/data/home/mass/fscucchia/Bermuda/output/Variant_Filtration"

### 
SNP_QD_MIN=20.000
SNP_MQ_MIN=50.000
SNP_FS_MAX=5.000
SNP_SOR_MAX=2.500

INDEL_QD_MIN=20.000
INDEL_MQ_MIN=45.000
INDEL_FS_MAX=5.000
INDEL_SOR_MAX=2.500
####

gatk VariantFiltration --reference "${REF}" --variant "${OUT}_SNPs.vcf.gz"   --output "${OUT}_SNPs_VarScores_filter.vcf.gz" \
	--filter-name "VarScores_filter_QD"  --filter-expression "QD < $SNP_QD_MIN" \
	--filter-name "VarScores_filter_MQ"  --filter-expression "MQ < $SNP_MQ_MIN" \
	--filter-name "VarScores_filter_FS"  --filter-expression "FS > $SNP_FS_MAX" \
	--filter-name "VarScores_filter_SOR" --filter-expression "SOR > $SNP_SOR_MAX" \
	1> "${OUT}_SNPs_VarScores_filter.vcf.gz.log" 2>&1


# gatk VariantFiltration --reference "${REF}" --variant "${OUT}_INDELs.vcf.gz" --output "${OUT}_INDELs_VarScores_filter.vcf.gz" \
# 	--filter-name "VarScores_filter_QD"  --filter-expression "QD < $INDEL_QD_MIN" \
# 	--filter-name "VarScores_filter_MQ"  --filter-expression "MQ < $INDEL_MQ_MIN" \
# 	--filter-name "VarScores_filter_FS"  --filter-expression "FS > $INDEL_FS_MAX" \
# 	--filter-name "VarScores_filter_SOR" --filter-expression "SOR > $INDEL_SOR_MAX" \
# 	1> "${OUT}_INDELs_VarScores_filter.vcf.gz.log" 2>&1
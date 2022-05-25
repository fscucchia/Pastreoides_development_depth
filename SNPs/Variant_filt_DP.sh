#!/bin/sh

F="/data/home/mass/fscucchia/Bermuda/output/SplitNCigarReads"
REF="/data/home/mass/fscucchia/databases/Pastreoides_genome_Kevin/past_filtered_assembly.fasta"
OUT="GVCFall"
O="/data/home/mass/fscucchia/Bermuda/output/Variant_Filtration"
###
DP_MIN=20.000
###

gatk VariantFiltration --reference "${REF}" --variant "${OUT}_SNPs_VarScores_filterPASSED.vcf" --output "${OUT}_SNPs_VarScores_filterPASSED_DPfilter.vcf" \
	        --genotype-filter-name "DP_filter" --genotype-filter-expression "DP < $DP_MIN" \
	        1> "${OUT}_SNPs_VarScores_filterPASSED_DPfilter.vcf.log" 2>&1
gatk VariantFiltration --reference "${REF}" --variant "${OUT}_INDELs_VarScores_filterPASSED.vcf" --output "${OUT}_INDELs_VarScores_filterPASSED_DPfilter.vcf" \
	        --genotype-filter-name "DP_filter" --genotype-filter-expression "DP < $DP_MIN" \
	        1> "${OUT}_INDELs_VarScores_filterPASSED_DPfilter.vcf.log" 2>&1
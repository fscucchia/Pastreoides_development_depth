#!/bin/bash

db_dir="/data/home/mass/fscucchia/databases/symbiont_combined_prot"

		if [ $1 -eq 1 ]; then 
			sbatch -N1 -n24 --ntasks-per-node=24 --mem=128000 -o $db_dir/run2.out -e $db_dir/run2.err \
   			-p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
			--wrap ". $CONDA; conda activate functional_annotation; diamond makedb --in corals_symbionts_combined_prot.fa -d corals_symbionts_combined_prot.fa_makedb; conda deactivate"
			 

        elif [ $1 -eq 2 ]; then 
			sbatch -N1 -n24 --ntasks-per-node=24 --mem=128000 -o $db_dir/run2.out -e $db_dir/run2.err \
   			-p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
			--wrap ". $CONDA; conda activate functional_annotation; diamond dbinfo -d corals_symbionts_combined_prot.fa_makedb.dmnd; conda deactivate"

		
		fi
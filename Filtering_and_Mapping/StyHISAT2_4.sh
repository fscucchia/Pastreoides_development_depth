H="/data/home/mass/fscucchia/Bermuda/output/HISAT"
F="/data/home/mass/fscucchia/Bermuda/scripts"
G="/data/home/mass/fscucchia/databases/Pastreoides_genome_KW"
A="/data/home/mass/fscucchia/Bermuda/output/filtered"
########################

if [ $1 -eq 1 ]; then
   mkdir -p $H
   sbatch -N1 -n1 --ntasks-per-node=1 --workdir=$H \
   -p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d \
	 --wrap ". $CONDA; conda activate newrnapipeline; ln -s /data/home/mass/fscucchia/Bermuda/output/filtered/*gz.filtered ./; conda deactivate"

elif [ $1 -eq 2 ]; then
     mkdir -p $G
     sbatch --mem=128000 -N1 -n20 --ntasks-per-node=20 --workdir=$G \
     -p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d \
  	 --wrap ". $CONDA; conda activate newrnapipeline; hisat2-build -f /data/home/mass/fscucchia/databases/Pastreoides_genome_KW/past_filtered_assembly.fasta ./Past_ref
	                   mkdir -p $F
					   chmod u+x StyHISAT_withSummary.sh
					   ./StyHISAT_withSummary.sh 
					   conda deactivate"

elif [ $1 -eq 3 ]; then # multiqc
		mkdir -p $A
		sbatch -N1 -n1 --ntasks-per-node=1 --workdir=$A -o "$A/multiqc_hisat/multiqc".out -e "$A/multiqc_hisat/multiqc".err \
		-p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d \
		--wrap ". $CONDA; conda activate multiqc_env; mkdir -p $A/multiqc_hisat; multiqc $A; conda deactivate"			
						
						
fi                        
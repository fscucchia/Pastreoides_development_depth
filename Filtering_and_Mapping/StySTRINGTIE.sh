H="/data/home/mass/fscucchia/Bermuda/output/StringTie"

########################

if [ $1 -eq 1 ]; then
   mkdir -p $H
   sbatch --mem=128000 -N1 -n20 --ntasks-per-node=20 --workdir=$H \
   -p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d \
	--wrap ". $CONDA; conda activate newrnapipeline
                     ln -s /data/home/mass/fscucchia/databases/Pastreoides_genome_KW/Pastreoides_all_v1.gff ./
                     chmod u+x StyStringTie_assembly.sh; ./StyStringTie_assembly.sh; conda deactivate"

fi
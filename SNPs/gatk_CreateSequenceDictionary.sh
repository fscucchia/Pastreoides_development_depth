
O="/data/home/mass/fscucchia/databases/Pastreoides_genome_Kevin"
R="/data/home/mass/fscucchia/databases/Pastreoides_genome_Kevin/past_filtered_assembly.fasta"

########################

if [ $1 -eq 1 ]; then
     sbatch -N1 -n1 --ntasks-per-node=1 --workdir=$O --job-name "CreateSeqDict" -o "$O/CreateSequenceDictionary.out" -e "$O/CreateSequenceDictionary.err" \
     -p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d \
  	 --wrap "gatk CreateSequenceDictionary -REFERENCE $R -OUTPUT ${R%*.*}.dict" 

elif [ $1 -eq 2 ]; then
     sbatch -N1 -n1 --ntasks-per-node=1 --workdir=$O --job-name "CreateSeqDict" -o "$O/CreateSequenceDictionary.out" -e "$O/CreateSequenceDictionary.err" \
     -p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d \
  	 --wrap "samtools faidx $R" 
						
fi                        
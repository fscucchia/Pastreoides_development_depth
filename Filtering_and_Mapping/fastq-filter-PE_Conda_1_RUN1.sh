
path1='/data/home/mass/fscucchia/Bermuda/output/concat/*_R1_concat.fastq.gz'   
dir2='/data/home/mass/fscucchia/Bermuda/output/filtered'
sed1='s/_R1_concat.fastq.gz/_R3_concat.fastq.gz/'

##################################################################

# check jobs status: sacct --state r ; sacct --state cd
# run on 5 nodes, totally 50 threads, 10 threads per node: srun -N5 -n50 --ntasks-per-node=10

if [ ! -d $dir2"/fastqc" ]; then mkdir $dir2"/fastqc"; fi
if [ ! -d $dir2"/fastq_screen" ]; then mkdir $dir2"/fastq_screen"; fi

i=0
j=0
for f1 in $path1; do
	let i=$i+1
	let j=$j+1
	f2=$(echo $f1 | sed "$sed1" )
    #f2=$(echo $f1 | sed 's/_R1.fastq.gz/_R3.fastq.gz/')
	echo ">"$i" "$j
	echo $f1
	if [ -e $f2 ] && [ $f1 != $f2 ]; then
		echo $f2" found"
		f1_=$dir2"/"$(echo $f1 | sed 's/.*\///')".filtered"
		f2_=$dir2"/"$(echo $f2 | sed 's/.*\///')".filtered"
		echo $f1_
		echo $f2_
		f1_unp=$f1_".unpair"
		f2_unp=$f2_".unpair"
		title=$(echo $f1 | sed 's/.*\///' | sed 's/.1.cor.fq//')
		echo 'title = '$title
		
		if [ $1 -eq 1 ]; then
			echo '1'
			sbatch -N1 -n1 --ntasks-per-node=1 --workdir=$dir2/fastqc \
				-p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
                --wrap ". $CONDA; conda activate multiqc_env; fastqc -o $dir2/fastqc $f1; conda deactivate" &
				# --wrap "module load FastQC/0.11.5-Java-1.8.0_74; fastqc -o $dir2/fastqc $f1" &
			sbatch -N1 -n1 --ntasks-per-node=1 --workdir=$dir2/fastqc \
				-p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
                --wrap ". $CONDA; conda activate multiqc_env; fastqc -o $dir2/fastqc $f2; conda deactivate" &
				#--wrap "module load FastQC/0.11.5-Java-1.8.0_74; fastqc -o $dir2/fastqc $f2" &
		elif [ $1 -eq 2 ]; then
			echo '2'
			qual_threshold1=25
			qual_window1=10
			minlen=20
            threads=24
			sbatch -N1 -n24 --ntasks-per-node=24 --mem=128000 -o $dir2/$title.out -e $dir2/$title.err \
                -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
				fastq-filter_Conda_job_1.sh \
					'PE' 'adapters4d.fa' $qual_threshold1 $qual_window1 $minlen $threads \
					$f1 $f2 $f1_ $f2_
		elif [ $1 -eq 3 ]; then
			echo '3'
			sbatch -N1 -n1 --ntasks-per-node=1 --workdir=$dir2/fastqc \
				--wrap ". $CONDA; conda activate multiqc_env; fastqc -o $dir2/fastqc $f1_; conda deactivate" &  
			sbatch -N1 -n1 --ntasks-per-node=1 --workdir=$dir2/fastqc \
				--wrap ". $CONDA; conda activate multiqc_env; fastqc -o $dir2/fastqc $f2_; conda deactivate" &
		
		fi
	fi
	#break
	echo "-----------------------"
done


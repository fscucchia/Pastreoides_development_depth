script_params='design/design_diamond'
outdir="/data/home/mass/fscucchia/Bermuda/output/Diamond/Diamond_prot5"
db="/data/home/mass/fscucchia/databases/symbiont_combined_prot/corals_symbionts_combined_prot.fa_makedb.dmnd"

############################################

function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

. "$script_params"
assert_

for i in ${!title[@]}; do
	echo "${title[i]}"
	echo "${for1[i]}"
	if [ -f "${for1[i]}" ]; then
		outIdx="$outdir/${title[i]}"
		q="${for1[i]}"
		echo "$outIdx"
		if [ $1 -eq 1 ]; then
		    # --more-sensitive
			if [ ! -f $outIdx.m8 ]; then
				echo OK
				t=20
				sbatch -N 1 -n"$t" --ntasks-per-node="$t" --mem=100000 -o "$outIdx.out" -e "$outIdx.err" \
				-p hive7d,hiveunlim,queen,preempt7d,preempt31d,hive1d,preempt1d \
				--wrap "echo diamond
						function assert_ { rc=\$?; if [[ \$rc != 0 ]]; then echo 'exit !!!!'; exit \$rc; fi }
						. \"$CONDA\"; assert_
						conda activate functional_annotation
						diamond blastx -a \"$outIdx\" --tmpdir \"$outdir\" --threads \"$t\" -d \"$db\" -q \"$q\" --index-chunks 1 --top 10 --evalue 0.01
						assert_
						diamond view -a $outIdx.daa -o $outIdx.m8; assert_
						rm $outIdx.daa; assert_
						conda deactivate
						"
			fi
	fi

done




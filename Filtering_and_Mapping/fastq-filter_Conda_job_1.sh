#!/bin/sh

option="$1" # SE/PE
adapters_file="$2"
qual_threshold1="$3"
qual_window1="$4"
MINLEN="$5"
threads="$6"

# SE:
# fastq-filter_job_1.sh $option $adapters_file $qual_threshold1 $qual_window1 $trimommatic $MINLEN $fq_in1 $fq_out1
# fastq-filter_job_1.sh SE primers.fa 27 50 $trimommatic 200 $fq_in1
#
# PE:
# 
# fastq-filter_job_1.sh $option $adapters_file $qual_threshold1 $qual_window1 $trimommatic $MINLEN $fq_in1 $fq_in2 $fq_out1 $fq_out2
# fastq-filter_job_1.sh PE primers.fa 27 50 $trimommatic 200 $fq_in1 $fq_in2

function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }


. "$CONDA"; assert_
conda activate qualityfilt; assert_

# another option: --quality-base=64

echo "option=$option"

if [ $option == 'SE' ]; then
	# SE :
	fq_in1="${7}"
	fq_out1="${8}"
	if [ "$fq_in1" != "$fq_out1" ] && [ -e "$fq_in1" ]; then
		# --match-read-wildcards  --format=fastq is not supported for parallel
		cutadapt --cores="$threads" -b file:"$adapters_file" --times 2 -o "$fq_out1.cutadapt" "$fq_in1"
		assert_
		trimmomatic SE -threads 1 -phred33 \
			      "$fq_out1.cutadapt" \
	              "$fq_out1" \
				  SLIDINGWINDOW:"$qual_window1":"$qual_threshold1" MINLEN:"$MINLEN"
		assert_
		rm "$fq_out1.cutadapt"
		assert_
	fi
elif [ $option == 'PE' ]; then
	# PE :
	echo 'xxxxxxxxxxx'
	fq_in1="${7}"
	fq_in2="${8}"
	fq_out1="${9}"
	fq_out2="${10}"
	if [ "$fq_in1" != "$fq_out1" ] && [ -e "$fq_in1" ] && \
	   [ "$fq_in2" != "$fq_out2" ] && [ -e "$fq_in2" ]; then
		# --match-read-wildcards --format=fastq  is not supported for parallel
		cutadapt --cores="$threads" -b file:"$adapters_file" -B file:"$adapters_file" --times 2 -o "$fq_out1.cutadapt" -p "$fq_out2.cutadapt" "$fq_in1" "$fq_in2"
		assert_
		trimmomatic PE -threads 1 -phred33 \
			      "$fq_out1.cutadapt" \
		          "$fq_out2.cutadapt" \
			      "$fq_out1" "$fq_out1.unpair" \
			      "$fq_out2" "$fq_out2.unpair" \
				  SLIDINGWINDOW:"$qual_window1":"$qual_threshold1" MINLEN:"$MINLEN"
		assert_
		rm "$fq_out1.cutadapt"
		assert_
		rm "$fq_out2.cutadapt"
		assert_
	fi
fi

conda deactivate

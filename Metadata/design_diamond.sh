# files:
function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }
title[0]="D20"
for1[0]="/data/home/mass/fscucchia/Bermuda/output/filtered/D20_R1_concat.fastq.gz.filtered"
rev1[0]="/data/home/mass/fscucchia/Bermuda/output/filtered/D20_R3_concat.fastq.gz.filtered"
title[1]="D21"
for1[1]="/data/home/mass/fscucchia/Bermuda/output/filtered/D21_R1_concat.fastq.gz.filtered"
rev1[1]="/data/home/mass/fscucchia/Bermuda/output/filtered/D21_R3_concat.fastq.gz.filtered"
title[2]="D22"
for1[2]="/data/home/mass/fscucchia/Bermuda/output/filtered/D22_R1_concat.fastq.gz.filtered"
rev1[2]="/data/home/mass/fscucchia/Bermuda/output/filtered/D22_R3_concat.fastq.gz.filtered"
title[3]="L7"
for1[3]="/data/home/mass/fscucchia/Bermuda/output/filtered/L7_R1_concat.fastq.gz.filtered"
rev1[3]="/data/home/mass/fscucchia/Bermuda/output/filtered/L7_R3_concat.fastq.gz.filtered"
title[4]="L8"
for1[4]="/data/home/mass/fscucchia/Bermuda/output/filtered/L8_R1_concat.fastq.gz.filtered"
rev1[4]="/data/home/mass/fscucchia/Bermuda/output/filtered/L8_R3_concat.fastq.gz.filtered"
title[5]="L9"
for1[5]="/data/home/mass/fscucchia/Bermuda/output/filtered/L9_R1_concat.fastq.gz.filtered"
rev1[5]="/data/home/mass/fscucchia/Bermuda/output/filtered/L9_R3_concat.fastq.gz.filtered"
title[6]="R13B"
for1[6]="/data/home/mass/fscucchia/Bermuda/output/filtered/R13B_R1_concat.fastq.gz.filtered"
rev1[6]="/data/home/mass/fscucchia/Bermuda/output/filtered/R13B_R3_concat.fastq.gz.filtered"
title[7]="R15B"
for1[7]="/data/home/mass/fscucchia/Bermuda/output/filtered/R15B_R1_concat.fastq.gz.filtered"
rev1[7]="/data/home/mass/fscucchia/Bermuda/output/filtered/R15B_R3_concat.fastq.gz.filtered"
title[8]="R33"
for1[8]="/data/home/mass/fscucchia/Bermuda/output/filtered/R33_R1_concat.fastq.gz.filtered"
rev1[8]="/data/home/mass/fscucchia/Bermuda/output/filtered/R33_R3_concat.fastq.gz.filtered"
title[9]="S4"
for1[9]="/data/home/mass/fscucchia/Bermuda/output/filtered/S4_R1_concat.fastq.gz.filtered"
rev1[9]="/data/home/mass/fscucchia/Bermuda/output/filtered/S4_R3_concat.fastq.gz.filtered"
title[10]="S5"
for1[10]="/data/home/mass/fscucchia/Bermuda/output/filtered/S5_R1_concat.fastq.gz.filtered"
rev1[10]="/data/home/mass/fscucchia/Bermuda/output/filtered/S5_R3_concat.fastq.gz.filtered"
title[11]="S7"
for1[11]="/data/home/mass/fscucchia/Bermuda/output/filtered/S7_R1_concat.fastq.gz.filtered"
rev1[11]="/data/home/mass/fscucchia/Bermuda/output/filtered/S7_R3_concat.fastq.gz.filtered"



#run:
#for i in ${!title[@]}; do
#	echo $i
#	echo ${title[i]}
#	echo ${for1[i]}
#	echo ${rev1[i]}
#	echo '----------------'
#done

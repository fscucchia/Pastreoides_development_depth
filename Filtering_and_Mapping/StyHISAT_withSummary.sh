#!/bin/sh

F="/data/home/mass/fscucchia/Bermuda/output/filtered"
array1=($(ls $F/*_R1_concat.fastq.gz.filtered))
for i in ${array1[@]}; do
        hisat2 -p 8 --new-summary --rf --dta -q -x Past_ref -1 ${i} -2 $(echo ${i}|sed s/_R1/_R3/) -S ${i}.sam --summary-file ${i}.txt 
        samtools sort -@ 8 -o ${i}.bam ${i}.sam
    		echo "${i}_bam"
        rm ${i}.sam
        echo "HISAT2 PE ${i}" $(date)
done



#!/bin/sh

# F="/data/home/mass/fscucchia/Bermuda/output/filtered"
# array1=($(ls $F/*_R1_concat.fastq.gz.filtered))
# for i in ${array1[@]}; do
#         gatk FastqToSam \
# 				--FASTQ ${i} \
# 				--FASTQ2 $(echo ${i}|sed s/_R1/_R3/) \
# 				--OUTPUT ${i}.FastqToSam.unmapped.bam \
# 				--SAMPLE_NAME ${i}; touch ${i}.FastqToSam.done
# done

F="/data/home/mass/fscucchia/Bermuda/output/filtered"
array1=($(ls $F/*_R1_concat.fastq.gz.filtered))
for i in ${array1[@]}; do
        gatk FastqToSam \
				--FASTQ ${i} \
				--OUTPUT ${i}.FastqToSam.unmapped.bam \
				--SAMPLE_NAME ${i}; touch ${i}.FastqToSam.done
done
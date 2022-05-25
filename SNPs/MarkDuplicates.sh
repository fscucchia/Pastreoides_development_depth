
#!/bin/sh

G="/data/home/mass/fscucchia/Bermuda/output/filtered"

if [ $1 -eq 1 ]; then
     sbatch --mem=300000 -N1 -n20 --ntasks-per-node=20 --workdir=$G --job-name "MarkDuplicates" -o "$G/MarkDuplicates.out" -e "$G/MarkDuplicates.err" \
     -p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d \
  	 --wrap "cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk MarkDuplicates --INPUT S7_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --OUTPUT S7_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE S7_R1_MergeBamAlignment.MarkDuplicates.metrics
cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk MarkDuplicates --INPUT S4_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --OUTPUT S4_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE S4_R1_MergeBamAlignment.MarkDuplicates.metrics
cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk MarkDuplicates --INPUT S5_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --OUTPUT S5_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE S5_R1_MergeBamAlignment.MarkDuplicates.metrics
cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk MarkDuplicates --INPUT R33_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --OUTPUT R33_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE R33_R1_MergeBamAlignment.MarkDuplicates.metrics
cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk MarkDuplicates --INPUT R15B_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --OUTPUT R15B_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE R15B_R1_MergeBamAlignment.MarkDuplicates.metrics
cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk MarkDuplicates --INPUT R13B_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --OUTPUT R13B_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE R13B_R1_MergeBamAlignment.MarkDuplicates.metrics
cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk MarkDuplicates --INPUT L7_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --OUTPUT L7_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE L7_R1_MergeBamAlignment.MarkDuplicates.metrics
cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk MarkDuplicates --INPUT L8_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --OUTPUT L8_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE L8_R1_MergeBamAlignment.MarkDuplicates.metrics
cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk MarkDuplicates --INPUT L9_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --OUTPUT L9_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE L9_R1_MergeBamAlignment.MarkDuplicates.metrics
cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk MarkDuplicates --INPUT D20_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --OUTPUT D20_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE D20_R1_MergeBamAlignment.MarkDuplicates.metrics
cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk MarkDuplicates --INPUT D21_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --OUTPUT D21_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE D21_R1_MergeBamAlignment.MarkDuplicates.metrics
cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk MarkDuplicates --INPUT D22_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --OUTPUT D22_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE D22_R1_MergeBamAlignment.MarkDuplicates.metrics"

fi    
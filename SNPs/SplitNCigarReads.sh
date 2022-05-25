
#!/bin/sh

G="/data/home/mass/fscucchia/Bermuda/output/filtered"
O="/data/home/mass/fscucchia/Bermuda/output/SplitNCigarReads"
R="/data/home/mass/fscucchia/databases/Pastreoides_genome_Kevin/past_filtered_assembly.fasta"

if [ $1 -eq 1 ]; then
     sbatch -N1 -n24 --ntasks-per-node=24 --mem=300000 --workdir=$G --job-name "SplitNCigarReads" -o "$G/SplitNCigarReads.out" -e "$G/SplitNCigarReads.err" \
     -p hive1d,hive7d,preempt1d,preempt7d,preempt31d,hiveunlim,queen \
     --wrap "cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk SplitNCigarReads -R $R -I R13B_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam -O $O/R13B_R1.SplitNCigarReads.split.bam
cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk SplitNCigarReads -R $R -I L7_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam -O $O/L7_R1.SplitNCigarReads.split.bam
cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk SplitNCigarReads -R $R -I L8_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam -O $O/L8_R1.SplitNCigarReads.split.bam
cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk SplitNCigarReads -R $R -I L9_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam -O $O/L9_R1.SplitNCigarReads.split.bam
cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk SplitNCigarReads -R $R -I D20_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam -O $O/D20_R1.SplitNCigarReads.split.bam
cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk SplitNCigarReads -R $R -I D21_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam -O $O/D21_R1.SplitNCigarReads.split.bam
cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk SplitNCigarReads -R $R -I R15B_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam -O $O/R15B_R1.SplitNCigarReads.split.bam
cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk SplitNCigarReads -R $R -I R33_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam -O $O/R33_R1.SplitNCigarReads.split.bam
cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk SplitNCigarReads -R $R -I D22_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam -O $O/D22_R1.SplitNCigarReads.split.bam
cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk SplitNCigarReads -R $R -I S4_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam -O $O/S4_R1.SplitNCigarReads.split.bam
cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk SplitNCigarReads -R $R -I S5_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam -O $O/S5_R1.SplitNCigarReads.split.bam
cd /data/home/mass/fscucchia/Bermuda/output/filtered
             gatk SplitNCigarReads -R $R -I S7_R1_MergeBamAlignment.MarkDuplicates.dedupped.bam -O $O/S7_R1.SplitNCigarReads.split.bam

fi    
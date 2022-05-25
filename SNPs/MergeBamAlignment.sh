
#!/bin/sh

# F="/data/home/mass/fscucchia/Bermuda/output/filtered"
# G="/data/home/mass/fscucchia/databases/Pastreoides_genome_Kevin/past_filtered_assembly.fasta"
#   array1=($(ls $F/*_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam))
#   for i in ${array1[@]}; do
#                 gatk MergeBamAlignment --REFERENCE_SEQUENCE $G --UNMAPPED_BAM ${i}.FastqToSam.unmapped.rg.bam --ALIGNED_BAM ${i}_R1_concat.fastq.gz.filtered.bam --OUTPUT ${i}.MergeBamAlignment.merged.bam --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT;
# 				touch ${i}.MergeBamAlignment.done

#   done

G="/data/home/mass/fscucchia/databases/Pastreoides_genome_Kevin/past_filtered_assembly.fasta"

cd /data/home/mass/fscucchia/Bermuda/output/filtered
gatk MergeBamAlignment --REFERENCE_SEQUENCE $G --UNMAPPED_BAM S7_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam --ALIGNED_BAM S7_R1_concat.fastq.gz.filtered.bam --OUTPUT S7_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT;
touch S7_R1_concat.fastq.gz.filtered.MergeBamAlignment.done

cd /data/home/mass/fscucchia/Bermuda/output/filtered
gatk MergeBamAlignment --REFERENCE_SEQUENCE $G --UNMAPPED_BAM S4_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam --ALIGNED_BAM S4_R1_concat.fastq.gz.filtered.bam --OUTPUT S4_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT;
touch S4_R1_concat.fastq.gz.filtered.MergeBamAlignment.done

cd /data/home/mass/fscucchia/Bermuda/output/filtered
gatk MergeBamAlignment --REFERENCE_SEQUENCE $G --UNMAPPED_BAM S5_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam --ALIGNED_BAM S5_R1_concat.fastq.gz.filtered.bam --OUTPUT S5_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT;
touch S5_R1_concat.fastq.gz.filtered.MergeBamAlignment.done

cd /data/home/mass/fscucchia/Bermuda/output/filtered
gatk MergeBamAlignment --REFERENCE_SEQUENCE $G --UNMAPPED_BAM R33_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam --ALIGNED_BAM R33_R1_concat.fastq.gz.filtered.bam --OUTPUT R33_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT;
touch R33_R1_concat.fastq.gz.filtered.MergeBamAlignment.done

cd /data/home/mass/fscucchia/Bermuda/output/filtered
gatk MergeBamAlignment --REFERENCE_SEQUENCE $G --UNMAPPED_BAM R15B_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam --ALIGNED_BAM R15B_R1_concat.fastq.gz.filtered.bam --OUTPUT R15B_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT;
touch R15B_R1_concat.fastq.gz.filtered.MergeBamAlignment.done

cd /data/home/mass/fscucchia/Bermuda/output/filtered
gatk MergeBamAlignment --REFERENCE_SEQUENCE $G --UNMAPPED_BAM R13B_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam --ALIGNED_BAM R13B_R1_concat.fastq.gz.filtered.bam --OUTPUT R13B_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT;
touch R13B_R1_concat.fastq.gz.filtered.MergeBamAlignment.done

cd /data/home/mass/fscucchia/Bermuda/output/filtered
gatk MergeBamAlignment --REFERENCE_SEQUENCE $G --UNMAPPED_BAM L7_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam --ALIGNED_BAM L7_R1_concat.fastq.gz.filtered.bam --OUTPUT L7_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT;
touch L7_R1_concat.fastq.gz.filtered.MergeBamAlignment.done

cd /data/home/mass/fscucchia/Bermuda/output/filtered
gatk MergeBamAlignment --REFERENCE_SEQUENCE $G --UNMAPPED_BAM L8_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam --ALIGNED_BAM L8_R1_concat.fastq.gz.filtered.bam --OUTPUT L8_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT;
touch L8_R1_concat.fastq.gz.filtered.MergeBamAlignment.done

cd /data/home/mass/fscucchia/Bermuda/output/filtered
gatk MergeBamAlignment --REFERENCE_SEQUENCE $G --UNMAPPED_BAM L9_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam --ALIGNED_BAM L9_R1_concat.fastq.gz.filtered.bam --OUTPUT L9_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT;
touch L9_R1_concat.fastq.gz.filtered.MergeBamAlignment.done

cd /data/home/mass/fscucchia/Bermuda/output/filtered
gatk MergeBamAlignment --REFERENCE_SEQUENCE $G --UNMAPPED_BAM D20_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam --ALIGNED_BAM D20_R1_concat.fastq.gz.filtered.bam --OUTPUT D20_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT;
touch D20_R1_concat.fastq.gz.filtered.MergeBamAlignment.done

cd /data/home/mass/fscucchia/Bermuda/output/filtered
gatk MergeBamAlignment --REFERENCE_SEQUENCE $G --UNMAPPED_BAM D21_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam --ALIGNED_BAM D21_R1_concat.fastq.gz.filtered.bam --OUTPUT D21_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT;
touch D21_R1_concat.fastq.gz.filtered.MergeBamAlignment.done

cd /data/home/mass/fscucchia/Bermuda/output/filtered
gatk MergeBamAlignment --REFERENCE_SEQUENCE $G --UNMAPPED_BAM D22_R1_concat.fastq.gz.filtered.FastqToSam.unmapped.rg.bam --ALIGNED_BAM D22_R1_concat.fastq.gz.filtered.bam --OUTPUT D22_R1_concat.fastq.gz.filtered.MergeBamAlignment.merged.bam --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT;
touch D22_R1_concat.fastq.gz.filtered.MergeBamAlignment.done
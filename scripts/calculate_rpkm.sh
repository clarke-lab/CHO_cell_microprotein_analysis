#!/bin/bash
mkdir -p uORF_analysis/reference_counts

# count the CHX riboseq data
for i in nts_r1 nts_r2 nts_r3 nts_r4 ts_r1 ts_r2 ts_r3 ts_r4
do
cs count \
  --count_files data/riboseq_chx/mapped/individual/$i".bam" \
  --countfile_format BAM \
  --fiveprime \
  --offset 12 \
  diff_translation_analysis/reference_annotation_transcript.positions \
  uORF_analysis/reference_counts/riboseq_$i
done


# count the RNASeq data
for i in nts_r1 nts_r2 nts_r3 nts_r4 ts_r1 ts_r2 ts_r3 ts_r4
do
cs count \
  --count_files data/rnaseq_se/mapped/individual/$i".bam" \
  --countfile_format BAM \
  --fiveprime \
  ./diff_translation_analysis/reference_annotation_transcript.positions \
  ./uORF_analysis/reference_counts/rnaseq_$i
done
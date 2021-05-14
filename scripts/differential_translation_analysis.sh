#!/bin/bash

mkdir -p diff_translation_analysis/reference_annotation
gtf2bed < reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf > \
bed diff_translation_analysis/reference_annotation/cgr_reference.bed

python scripts/create_plastid_mask_reference.py

cs generate \
--annotation_files diff_translation_analysis/reference_annotation/cgr_reference.bed  \
--annotation_format BED \
--sorted \
./diff_translation_analysis/reference_annotation/reference_annotation \
--mask_annotation_files \
    diff_translation_analysis/reference_annotation/reference_start_codon_masks.bed \
    diff_translation_analysis/reference_annotation/reference_stop_codon_masks.bed \
--mask_annotation_format BED


mkdir diff_translation_analysis/reference_counts

# count the CHX riboseq data
for i in nts_r1 nts_r2 nts_r3 nts_r4 ts_r1 ts_r2 ts_r3 ts_r4
do
cs count \
  --count_files data/riboseq_chx/mapped/individual/$i".bam" \
  --countfile_format BAM \
  --fiveprime \
  --offset 12 \
  diff_translation_analysis/reference_annotation/reference_annotation_gene.positions \
  diff_translation_analysis/reference_counts/riboseq_$i
done

# count the RNASeq data
for i in nts_r1 nts_r2 nts_r3 nts_r4 ts_r1 ts_r2 ts_r3 ts_r4
do
cs count \
  --count_files data/rnaseq_se/mapped/individual/$i".bam" \
  --countfile_format BAM \
  --fiveprime \
  ./diff_translation_analysis/reference_annotation/reference_annotation_gene.positions \
  ./diff_translation_analysis/reference_counts/rnaseq_$i
done



mkdir -p diff_translation_analysis/orf_annotation

python scripts/create_plastid_mask_orfrater.py

cs generate \
--annotation_files ./orfrater_analysis/orfrater_predictions.reference.bed  \
--annotation_format BED \
--sorted \
./diff_translation_analysis/orf_annotation/orfrater_annotation \
--mask_annotation_files \
    diff_translation_analysis/orf_annotation/orfrater_start_codon_masks.bed \
    diff_translation_analysis/orf_annotation/orfrater_stop_codon_masks.bed \
--mask_annotation_format BED

mkdir diff_translation_analysis/orfrater_counts

# count the CHX riboseq data
for i in nts_r1 nts_r2 nts_r3 nts_r4 ts_r1 ts_r2 ts_r3 ts_r4
do
cs count \
  --count_files data/riboseq_chx/mapped/individual/$i".bam" \
  --countfile_format BAM \
  --fiveprime \
  --offset 12 \
  diff_translation_analysis/orf_annotation/orfrater_annotation_gene.positions \
  diff_translation_analysis/orfrater_counts/riboseq_$i
done

# count the RNASeq data
for i in nts_r1 nts_r2 nts_r3 nts_r4 ts_r1 ts_r2 ts_r3 ts_r4
do
cs count \
  --count_files data/rnaseq_se/mapped/individual/$i".bam" \
  --countfile_format BAM \
  --fiveprime \
  ./diff_translation_analysis/orf_annotation/orfrater_annotation_gene.positions \
  ./diff_translation_analysis/orfrater_counts/rnaseq_$i
done

python scripts/combine_counts.py 

Rscript scripts/deseq_analysis.R
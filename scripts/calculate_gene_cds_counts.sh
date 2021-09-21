#!/bin/bash
#### Description: Gene level differential translation counting   
####              1. Count the CHX Ribo-seq RPFs in individual 
####                 samples for the TS and NTS groups
####              2. Count the reads for individual 
####                 samples for the TS and NTS groups
####              3. Combine for import into R    
#### 
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

# count the CHX riboseq data
mkdir diff_translation_analysis/reference_counts

for i in nts_r1 nts_r2 nts_r3 nts_r4 ts_r1 ts_r2 ts_r3 ts_r4
do
cs count \
  --count_files data/riboseq_chx/mapped/individual/$i".bam" \
  --countfile_format BAM \
  --fiveprime \
  --offset 12 \
  --min_length 28 \
  --max_length 31 \
  diff_translation_analysis/reference_annotation_gene.positions \
  diff_translation_analysis/reference_counts/riboseq_$i
done

# count the RNASeq data
for i in nts_r1 nts_r2 nts_r3 nts_r4 ts_r1 ts_r2 ts_r3 ts_r4
do
cs count \
  --count_files data/rnaseq_se/mapped/individual/$i".bam" \
  --countfile_format BAM \
  --fiveprime \
  ./diff_translation_analysis/reference_annotation_gene.positions \
  ./diff_translation_analysis/reference_counts/rnaseq_$i
done

# conda deactivate

# merge the RNASeq and Ribo-seq counts into a single table
python scripts/combine_counts.py 
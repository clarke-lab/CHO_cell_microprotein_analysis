#!/bin/bash

mkdir -p diff_translation_analysis

# parse the annotated protein coding genes from the NCBI reference
grep 'protein_coding\|XM_\|NM_\|NP_' \
reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf > \
diff_translation_analysis/protein_coding_reference.tmp.gtf

# remove peakfrac transcripts from the count
awk -F'\t' '$3 == "peakfrac" {print $1}' orfrater_analysis/tid_removal_summary.txt > \
diff_translation_analysis/peakfrac.txt 

grep -Fv diff_translation_analysis/peakfrac.txt diff_translation_analysis/protein_coding_reference.tmp.gtf > \
diff_translation_analysis/protein_coding_reference.gtf 

rm diff_translation_analysis/protein_coding_reference.tmp.gtf

# convert to bed format
gtf2bed --gtf diff_translation_analysis/protein_coding_reference.gtf \
--bed diff_translation_analysis/cgr_reference.bed

# selected novel ofs comes from R - filtered for potential false positive "new" annotation
 grep -f diff_translation_analysis/selected_novel_orfs.txt \
 orfrater_analysis/orfrater_predictions.reference.gtf > \
 diff_translation_analysis/cgr_new_orfs.gtf

# convert to bed format
grep -f diff_translation_analysis/selected_novel_orfs.txt \
orfrater_analysis/orfrater_predictions.reference.bed > \
 diff_translation_analysis/cgr_new_orfs.bed

# join the reference and novel annotations
 cat diff_translation_analysis/cgr_new_orfs.gtf \
 diff_translation_analysis/protein_coding_reference.gtf > \
 diff_translation_analysis/diff_translation.gtf

cat diff_translation_analysis/cgr_new_orfs.bed diff_translation_analysis/cgr_reference.bed > \
diff_translation_analysis/diff_translation.bed

# use plastid to mask and generate the count reference
source ~/miniconda3/etc/profile.d/conda.sh
conda activate plastid
python scripts/create_plastid_mask.py

sort -k1,1 -k4,4n diff_translation_analysis/diff_translation.gtf | \
bgzip > diff_translation_analysis/diff_translation.gtf.gz 
tabix -p gff diff_translation_analysis/diff_translation.gtf.gz

sort -k1,1 -k2,2n -k3,3n  diff_translation_analysis/start_codon_masks.bed | \
bgzip  > diff_translation_analysis/start_codon_masks.bed.gz
tabix -p bed diff_translation_analysis/start_codon_masks.bed.gz

sort -k1,1 -k2,2n -k3,3n  diff_translation_analysis/stop_codon_masks.bed | \
bgzip  > diff_translation_analysis/stop_codon_masks.bed.gz
tabix -p bed diff_translation_analysis/stop_codon_masks.bed.gz


 cs generate \
 --annotation_files diff_translation_analysis/diff_translation.gtf  \
 --annotation_format GTF2 \
 ./diff_translation_analysis/reference_annotation \
 --mask_annotation_files \
     quantitation/plastid_reference/gene/start_codon_masks.bed \
     quantitation/plastid_reference/gene/stop_codon_masks.bed \
--mask_annotation_format BED 

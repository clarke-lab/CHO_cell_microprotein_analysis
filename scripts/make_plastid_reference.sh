#!/bin/bash
#### Description: Make a plastid reference for transcript quantitation
####              and gene-level counting   
####              1. Retain only annotated protein coding genes 
####                 from the NCBI PICR-H reference
####              2. Remove transcripts where reads accumulate at 
####                 small number of positions
####              3. Join the annotated PGCs with new ORF identifications 
####                 from ORF-RATER
####              4. Make a mask to elimate codons that may be biased due to 
####                 CHX
####              5. Apply mask and make count reference    
#### 
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

# parse the annotated protein coding genes from the NCBI reference
grep 'protein_coding\|XM_\|NM_\|NP_' \
reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf > \
plastid_reference/protein_coding_reference.tmp.gtf

# remove peakfrac transcripts from the count
awk -F'\t' '$3 == "peakfrac" {print $1}' orfrater_analysis/tid_removal_summary.txt > \
plastid_reference/peakfrac.txt 

grep -Fv plastid_reference/peakfrac.txt plastid_reference/protein_coding_reference.tmp.gtf > \
plastid_reference/protein_coding_reference.gtf 

rm plastid_reference/protein_coding_reference.tmp.gtf

# convert to bed format
#gtf2bed --gtf diff_translation_analysis/protein_coding_reference.gtf \
#--bed diff_translation_analysis/cgr_reference.bed

# selected novel ofs comes from R - filtered for potential false positive "new" annotation
 grep -f orf_lists/non_coding_orfs.txt \
 orfrater_analysis/orfrater_predictions.reference.gtf > \
 plastid_reference/cgr_new_orfs.gtf

# convert to bed format
# grep -f diff_translation_analysis/selected_novel_orfs.txt \
# orfrater_analysis/orfrater_predictions.reference.bed > \
#  diff_translation_analysis/cgr_new_orfs.bed

# join the reference and novel annotations
 cat plastid_reference/cgr_new_orfs.gtf \
 plastid_reference/protein_coding_reference.gtf > \
 plastid_reference/plastid_annotations.gtf

# cat plastid_reference/cgr_new_orfs.bed diff_translation_analysis/cgr_reference.bed > \
# diff_translation_analysis/diff_translation.bed

# use plastid to mask and generate the count reference
source ~/miniconda3/etc/profile.d/conda.sh
conda activate plastid

# mask codons
python scripts/create_plastid_mask.py

# sort -k1,1 -k4,4n plastid_reference/plastid_annotations.gtf | \
# bgzip > plastid_reference/plastid_annotations.gtf.gz 
# tabix -p gff plastid_reference/plastid_annotations.gtf.gz 

# sort -k1,1 -k2,2n -k3,3n  plastid_reference/start_codon_masks.bed | \
# bgzip  > plastid_reference/start_codon_masks.bed.gz
# tabix -p bed plastid_reference/start_codon_masks.bed.gz

# sort -k1,1 -k2,2n -k3,3n  plastid_reference/stop_codon_masks.bed | \
# bgzip  > plastid_reference/stop_codon_masks.bed.gz
# tabix -p plastid_reference/stop_codon_masks.bed.gz

 cs generate \
 --annotation_files plastid_reference/plastid_annotations.gtf  \
 --annotation_format GTF2 \
 ./plastid_reference/annotation \
 --mask_annotation_files \
     plastid_reference/start_codon_masks.bed \
     plastid_reference/stop_codon_masks.bed \
--mask_annotation_format BED 

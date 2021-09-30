#!/bin/bash
#### Description: calculate offsets for Ribo-seq   
####              data using Plastid              
####              the RNA-seq data is also analyse for comparison 
####              this analysis is carried out for 
####              1. individual 
####              2. merged samples   
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

source ~/miniconda3/etc/profile.d/conda.sh
conda activate plastid

mkdir rpf_psite_analysis

# make a Plastid annotation for the known CDS start positions
# across the genome
metagene generate rpf_psite_analysis/cgr_orfs \
--landmark cds_start \
--annotation_files reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf

# 1. individual sample psite offsets
for seqtype in riboseq_chx riboseq_harr riboseq_nd
  do
  
  mkdir -p rpf_psite_analysis/individual_files/$seqtype/offset
  mkdir rpf_psite_analysis/individual_files/$seqtype/periodicity

  while read -ra a ;
    do
      psite rpf_psite_analysis/cgr_orfs_rois.txt  \
      rpf_psite_analysis/individual_files/$seqtype/offset/${a[1]} \
      --min_length 28 \
      --max_length 31 \
      --require_upstream \
      --count_files data/$seqtype/mapped/individual/${a[1]}".bam"

      phase_by_size rpf_psite_analysis/cgr_orfs_rois.txt \
      rpf_psite_analysis/individual_files/$seqtype/periodicity/${a[1]} \
      --count_files data/$seqtype/mapped/individual/${a[1]}".bam" \
      --fiveprime  \
      --offset 12 \
      --codon_buffer 5 \
      --min_length 28 \
      --max_length 31

    done < data/"$seqtype".txt
done

# 2. merged samples psite offset
for seqtype in riboseq_chx riboseq_harr riboseq_nd rnaseq_se
do
    psite rpf_psite_analysis/cgr_orfs_rois.txt \
    rpf_psite_analysis/$seqtype \
    --min_length 28 \
    --max_length 31 \
    --require_upstream \
    --min_counts 10 \
    --count_files data/$seqtype/mapped/merged/$seqtype".bam"

done

# Determine the Proportion of reads found to be in frame at the 12 nt offset
for seqtype in riboseq_chx riboseq_harr riboseq_nd rnaseq_se
do
  phase_by_size rpf_psite_analysis/cgr_orfs_rois.txt \
  rpf_psite_analysis/$seqtype \
  --count_files data/$seqtype/mapped/merged/$seqtype".bam" \
  --fiveprime \
  --offset 12 \
  --codon_buffer 5 \
  --min_length 28 \
  --max_length 31
done

conda deactivate

#!/bin/bash
#### Description: Trancript level RPKM  
####              1. Get CHX Ribo-seq CDS RPKM values for 
####                 individual samples for the TS and NTS groups
####              2. Get RNA-seq CDS values individual 
####                 samples for the TS and NTS groups
####              3. Combine for import into R    
#### 
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

source ~/miniconda3/etc/profile.d/conda.sh
conda activate plastid

# count the CHX riboseq data
for i in nts_r1 nts_r2 nts_r3 nts_r4 ts_r1 ts_r2 ts_r3 ts_r4
do
cs count \
  --count_files data/riboseq_chx/mapped/individual/$i".bam" \
  --countfile_format BAM \
  --fiveprime \
  --offset 12 \
  plastid_reference/annotation_transcript.positions \
  quantitation/transcript_cds_rpkm/riboseq_$i
done


# count the RNASeq data
for i in nts_r1 nts_r2 nts_r3 nts_r4 ts_r1 ts_r2 ts_r3 ts_r4
do
cs count \
  --count_files data/rnaseq_se/mapped/individual/$i".bam" \
  --countfile_format BAM \
  --fiveprime \
  plastid_reference/annotation_transcript.positions \
  quantitation/transcript_cds_rpkm/rnaseq_$i
done

# merge the Ribo-seq and RNA-seq CDS RPKM
python scripts/combine_transcript_cds_rpkm.py

conda deactivate
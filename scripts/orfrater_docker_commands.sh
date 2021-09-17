
#!/bin/bash
#### Description: ORF-RATER analysis run via docker   
####                           
#### 
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

# set the working directory
work_dir=ribosome_footprint_profiling/

threads=50

prune_transcripts.py \
--inbed $work_dir/orfrater_analysis/cgr.orfrater.annotation.reference.bed \
--summarytable $work_dir/orfrater_analysis/tid_removal_summary.txt \
-p $threads \
--minlen 28 \
--maxlen 31 \
$work_dir/reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna \
$work_dir/data/riboseq_chx/mapped/merged/riboseq_chx.bam \
$work_dir/data/riboseq_nd/mapped/merged/riboseq_nd.bam \
--pseudogenes $work_dir/orfrater_analysis/pseudogene_transcript_ids.txt \
--outbed $work_dir/orfrater_analysis/transcripts.bed \
-v \
--force > $work_dir/orfrater_analysis/1prune.$source.log

make_tfams.py --force \
--inbed $work_dir/orfrater_analysis/transcripts.bed \
--tfamstem $work_dir/orfrater_analysis/tfams >   \
$work_dir/orfrater_analysis/2make_tfams.log

# # echo "identifying orf candidates"
FILE=$work_dir/orfrater_analysis/orf.h5
if test -f "$FILE"; then
    rm $FILE # need to remove this as problem with indexing if exists
fi    

find_orfs_and_types.py \
$work_dir/reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna \
--codons NTG \
-p $threads \
--force \
-v \
--tfamstem $work_dir/orfrater_analysis/tfams \
--inbed $work_dir/orfrater_analysis/transcripts.bed \
--orfstore $work_dir/orfrater_analysis/orf.h5 > \
$work_dir/orfrater_analysis/3find_ORF.log

mkdir $work_dir/orfrater_analysis/chx
psite_trimmed.py \
$work_dir/data/riboseq_chx/mapped/merged/riboseq_chx.bam \
--minrdlen 28 \
--maxrdlen 31 \
--subdir $work_dir/orfrater_analysis/chx \
--tallyfile tallies.txt \
-v \
--cdsbed $work_dir/orfrater_analysis/cgr.orfrater.annotation.reference.bed \
-p $threads \
--force > $work_dir/orfrater_analysis/psite_chx.log

mkdir $work_dir/orfrater_analysis/nd
psite_trimmed.py \
$work_dir/data/riboseq_nd/mapped/merged/riboseq_nd.bam \
--minrdlen 28 \
--maxrdlen 31 \
--subdir $work_dir/orfrater_analysis/nd \
--tallyfile tallies.txt \
--cdsbed $work_dir/orfrater_analysis/cgr.orfrater.annotation.reference.bed \
-p $threads \
-v \
--force > $work_dir/orfrater_analysis/psite_nd.log

mkdir $work_dir/orfrater_analysis/harr
psite_trimmed.py \
$work_dir/data/riboseq_harr/mapped/merged/riboseq_harr.bam \
--minrdlen 28 \
--maxrdlen 31 \
--subdir $work_dir/orfrater_analysis/harr \
--tallyfile tallies.txt \
--cdsbed $work_dir/orfrater_analysis/cgr.orfrater.annotation.reference.bed \
-p $threads \
-v \
--force > $work_dir/orfrater_analysis/psite_har.log

regress_orfs.py \
$work_dir/data/riboseq_harr/mapped/merged/riboseq_harr.bam \
--startonly \
--subdir $work_dir/orfrater_analysis/harr \
--orfstore $work_dir/orfrater_analysis/orf.h5 \
--inbed $work_dir/orfrater_analysis/transcripts.bed \
-p 50 \
-v \
--startcount 1 \
--force > $work_dir/orfrater_analysis/regress_start.log

regress_orfs.py \
$work_dir/data/riboseq_chx/mapped/merged/riboseq_chx.bam \
--subdir $work_dir/orfrater_analysis/chx \
--orfstore $work_dir/orfrater_analysis/orf.h5 \
--inbed $work_dir/orfrater_analysis/transcripts.bed \
--restrictbystarts $work_dir/orfrater_analysis/harr \
--startcount 1 \
-p $threads \
-v \
--force > chx_regress_stop.log

 regress_orfs.py \
 $work_dir/data/riboseq_nd/mapped/merged/riboseq_nd.bam \
 --subdir $work_dir/orfrater_analysis/nd \
 --orfstore $work_dir/orfrater_analysis/orf.h5 \
 --inbed $work_dir/orfrater_analysis/transcripts.bed \
 --restrictbystarts $work_dir/orfrater_analysis/harr \
 --startcount 1 \
 -p $threads \
 --force > $work_dir/orfrater_analysis/nodrug_regress_stop.log

rate_regression_output.py \
$work_dir/orfrater_analysis/harr \
$work_dir/orfrater_analysis/chx \
$work_dir/orfrater_analysis/nd \
--orfstore $work_dir/orfrater_analysis/orf.h5 \
--ratingsfile $work_dir/orfrater_analysis/orfratings.h5 \
-p $threads \
-v \
--CSV $work_dir/orfrater_analysis/rate_regression.ncbi.csv \
--force > $work_dir/orfrater_analysis/rate.regression.log

 make_orf_bed.py \
 --inbed $work_dir/orfrater_analysis/transcripts.bed \
 --ratingsfile $work_dir/orfrater_analysis/orfratings.h5 \
 --minlen 5 \
 --force \
 --outbed $work_dir/orfrater_analysis/orfrater_predictions.reference.bed

#  quantify_orfs.py \
#  $work_dir/data/riboseq_chx/mapped/individual/nts_r1.bam \
#  $work_dir/data/riboseq_chx/mapped/individual/nts_r2.bam \
#  $work_dir/data/riboseq_chx/mapped/individual/nts_r3.bam \
#  $work_dir/data/riboseq_chx/mapped/individual/nts_r4.bam \
#  $work_dir/data/riboseq_chx/mapped/individual/ts_r1.bam \
#  $work_dir/data/riboseq_chx/mapped/individual/ts_r2.bam \
#  $work_dir/data/riboseq_chx/mapped/individual/ts_r3.bam \
#  $work_dir/data/riboseq_chx/mapped/individual/ts_r4.bam \
#  --subdir $work_dir/orfrater_analysis/chx \
#  --inbed $work_dir/orfrater_analysis/transcripts.bed \
#  --startmask 1 2 \
#  --startmask 0 3 \
#  --ratingsfile $work_dir/orfrater_analysis/orfratings.h5 \
#  --minrating 0.3 \
#  --minlen 0 \
#  --quantfile test_rpkm \
#  --CSV  $work_dir/orfrater_analysis/orfrater_rpf_rpkm.csv \
#  --force \
#  -p 32

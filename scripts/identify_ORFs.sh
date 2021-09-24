#!/bin/bash
#### Description: Identify ORFs in the Chinese hamster genome
####              Using ORF-RATER
####              1. Remove transcripts from Pseudogenes
####              data using Plastid              
####              2. convert GTF to genePred format 
####              3. Remove CGR chromomosome that causes error
####              4. Run ORF-RATER in docker
####              5. convert ORF-RATER bed to GTF  
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

export PATH=$PATH:/mnt/HDD2/colin/bin/ORF-RATER
kent_path=/mnt/HDD2/colin/bin/kentUtils/bin/linux.x86_64/

mkdir orfrater_analysis

gtf=reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf

# 1. remove pseudogenes 
grep '; pseudo'  $gtf | \
awk -F'\'t '{print $9}' | awk -F';' '{print $1}' | sed s/gene_id//g | \
sed s/\"//g | sed 's/[[:space:]]//g'| sort | uniq > orfrater_analysis/pseudogene_gene_ids

grep -f orfrater_analysis/pseudogene_gene_ids $gtf | \
awk -F'\'t '{print $9}' | awk -F';' '{print $2}' | sed s/transcript_id//g | \
sed s/\"//g | sed 's/[[:space:]]//g'| uniq > orfrater_analysis/pseudogene_transcript_ids.txt

# 2. convert NCBI gtf to genePred format 
$kent_path/gtfToGenePred \
-ignoreGroupsWithoutExons -allErrors $gtf stdout | \
$kent_path/genePredToBed stdin \
orfrater_analysis/cgr.orfrater.annotation.tmp.bed

# 3. remove chromosome that cause error 
grep -v NW_023277000.1  orfrater_analysis/cgr.orfrater.annotation.tmp.bed > \
orfrater_analysis/cgr.orfrater.annotation.reference.bed

rm orfrater_analysis/cgr.orfrater.annotation.tmp.bed

# 4. run the ORF-RATER docker with the commands file
docker run --rm \
-v /mnt/HDD2/colin/ribosome_footprint_profiling/:/ribosome_footprint_profiling \
-t orfrater:final \
bash "ribosome_footprint_profiling/scripts/orfrater_docker_commands.sh"

# 5. convert orfrater BED to GTF
/mnt/HDD2/colin/bin/kentUtils/bin/linux.x86_64/bedToGenePred \
orfrater_analysis/orfrater_predictions.reference.bed stdout | \
/mnt/HDD2/colin/bin/kentUtils/bin/linux.x86_64/genePredToGtf file stdin \
orfrater_analysis/orfrater_predictions.reference.gtf
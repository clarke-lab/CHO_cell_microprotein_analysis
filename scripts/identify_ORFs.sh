#!/bin/bash

export PATH=$PATH:/mnt/HDD2/colin/bin/ORF-RATER
kent_path=/mnt/HDD2/colin/bin/kentUtils/bin/linux.x86_64/

mkdir orfrater_analysis

gtf=reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf

# remove pseudogenes 
grep '; pseudo'  $gtf | \
awk -F'\'t '{print $9}' | awk -F';' '{print $1}' | sed s/gene_id//g | \
sed s/\"//g | sed 's/[[:space:]]//g'| sort | uniq > orfrater_analysis/pseudogene_gene_ids

grep -f pseudogene_gene_ids $gtf | \
awk -F'\'t '{print $9}' | awk -F';' '{print $2}' | sed s/transcript_id//g | \
sed s/\"//g | sed 's/[[:space:]]//g'| uniq > orfrater_analysis/pseudogene_transcript_ids.txt

# convert NCBI gtf to genePred format 
$kent_path/gtfToGenePred \
-ignoreGroupsWithoutExons -allErrors $gtf stdout | \
$kent_path/genePredToBed stdin \
orfrater_analysis/cgr.orfrater.annotation.tmp.bed

# remove chromosome that cause error 
grep -v NW_023277000.1  orfrater_analysis/cgr.orfrater.annotation.tmp.bed > \
orfrater_analysis/cgr.orfrater.annotation.reference.bed

rm orfrater_analysis/cgr.orfrater.annotation.tmp.bed

# run the ORF-RATER docker with the commands file
docker run --rm \
-v /mnt/HDD2/colin/ribosome_footprint_profiling/:/ribosome_footprint_profiling \
-t orfrater:final \
bash "ribosome_footprint_profiling/scripts/orfrater_docker_commands.sh"

# convert orfrater
/mnt/HDD2/colin/bin/kentUtils/bin/linux.x86_64/bedToGenePred \
orfrater_analysis/orfrater_predictions.reference.bed stdout | \
/mnt/HDD2/colin/bin/kentUtils/bin/linux.x86_64/genePredToGtf file stdin \
orfrater_analysis/orfrater_predictions.reference.gtf
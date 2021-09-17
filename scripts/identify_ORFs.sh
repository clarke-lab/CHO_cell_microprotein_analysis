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

# filter the ORFs - list gerated in R [orfscore > 0.6, AAlen >= 10, no annotated, no truncation] 

 cat ../riboseq_work/selected_novel_orfs.txt | awk '{print $1}' > orfrater_analysis/selected_orfnames.txt

grep -f orfrater_analysis/selected_orfnames.txt  \
orfrater_analysis/orfrater_predictions.reference.bed > \
orfrater_analysis/novel_proteoforms.bed 

# # Make a fasta file for proteomics analysis
 

# # make a nucleotide fasta for the CDS gregion - from the orf-rater
cat orfrater_analysis/novel_proteoforms.bed | awk -v type=CDS -f scripts/bed12toAnnotation.awk > orfrater_analysis/novel_orfrater_predictions.bed12

bedtools getfasta -fi reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna \
 -bed orfrater_analysis/novel_orfrater_predictions.bed12 \
 -fo orfrater_analysis/orfrater_predictions.nucleotide.fasta -split -s -name

# # convert to protein fasta
 transeq -seq orfrater_analysis/orfrater_predictions.nucleotide.fasta \
 -frame 1 \
 -outseq orfrater_analysis/orfrater_predictions.protein.fasta

# # Download the UNIPROT proteins and join to novel ORF-RATER for MS
wget  \
ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000001075/UP000001075_10029.fasta.gz -P reference_genome \
&& gunzip reference_genome/UP000001075_10029.fasta.gz
 cat orfrater_analysis/orfrater_predictions.protein.fasta reference_genome/UP000001075_10029.fasta > orfrater_analysis/extended_cgr_database.fasta


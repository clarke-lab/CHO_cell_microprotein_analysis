#!/usr/bin/env bash
star_path=../bin/STAR-2.7.2d/bin/Linux_x86_64



for seqtype in riboseq_chx total_rnaseq_se
do
mkdir -p data/$seqtype/mapped/individual
  while read -ra a ;
  do
    echo ${a[1]}
    $star_path/STAR \
      --outSAMtype BAM SortedByCoordinate \
      --runThreadN 32 \
      --outFilterMismatchNmax 2 \
      --genomeDir reference_genome/star_index_stringtie \
      --readFilesIn  data/$seqtype/preproccessed_data/${a[1]}".fastq" \
      --outFileNamePrefix data/$seqtype/mapped/individual/${a[1]} \
      --outFilterMultimapNmax 16 \
      --outFilterMatchNmin 16 \
      --alignEndsType EndToEnd \
      --outMultimapperOrder Random \
      --outSAMattributes All

      mv data/$seqtype/mapped/individual/${a[1]}Aligned.sortedByCoord.out.bam >\
      data/$seqtype/mapped/individual/${a[1]}.bam

      samtools index data/$seqtype/mapped/individual/${a[1]}.bam
  done <  data/"$seqtype".txt
done

while read -ra a ;
do
  echo ${a[0]}
  $star_path/STAR \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN 8 \
    --outFilterMismatchNmax 2 \
    --genomeDir reference_genome/star_index_ncbi \
    --readFilesIn  data/rnaseq_se/preprocessed_data/${a[0]} \
    --outFileNamePrefix data/rnaseq_se/mapped/individual/${a[1]} \
    --outFilterMultimapNmax 16 \
    --outFilterMatchNmin 16 \
    --alignEndsType EndToEnd \
    --outMultimapperOrder Random \
    --readFilesCommand zcat \
    --quantMode TranscriptomeSAM \
    --outSAMattributes All

    mv data/$seqtype/rnaseq_se/individual/${a[1]}Aligned.sortedByCoord.out.bam data/$seqtype/mapped/individual/${a[1]}.bam
  
   samtools index data/rnaseq_se/mapped/individual/${a[1]}.bam
done <  data/total_rnaseq_se.txt

# Count the SE data
mkdir xtail_analysis/gene/counts/rna
while read -ra a ; do
  htseq-count data/rnaseq_se/mapped/individual/${a[1]}"Aligned.sortedByCoord.out.bam" \
  orfrater_analysis/cgr_combined_annotation.clean.gtf \
  -m union \
  -i gene_id \
  --nonunique all \
  -s no \
  -f bam > xtail_analysis/gene/counts/rna/${a[1]}".count" &
done  < data/total_rnaseq_se.txt

# Count the CHX RPFs
mkdir -p xtail_analysis/gene/counts/rpf
while read -ra a ; do
   python scripts/xtail_count_rpf.py data/riboseq_chx/mapped/individual/${a[1]}"Aligned.sortedByCoord.out.bam" \
   orfrater_analysis/cgr_combined_annotation.clean.gtf > \
   xtail_analysis/gene/counts/rpf/${a[1]}".count" &
done  < data/riboseq_chx.txt

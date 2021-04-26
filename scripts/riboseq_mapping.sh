#!/bin/bash

star_path=../bin/STAR-2.7.2d/bin/Linux_x86_64

for treatment in riboseq_chx riboseq_harr riboseq_nodrug
  mkdir data/$treatment/mapped
  mkdir data/$treatment/bigwig

while read -ra a ;
  $star_path/STAR \
  --outFilterType BySJout \
  --runThreadN 32 \
  --outFilterMismatchNmax 2 \
  --genomeDir reference_genome/star_index_stringtie \
  --readFilesIn data/$treatment/preprocessed/merged/$treatment.fastq \
  --outFileNamePrefix data/$treatment/preprocessed/merged/$treatment \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode TranscriptomeSAM GeneCounts \
  --outFilterMultimapNmax 1 \
  --outFilterMatchNmin 16 \
  --alignEndsType EndToEnd \
  --outSAMattributes All


bamCoverage \
-b merged_bam_files/chx.bam \
-o data/$treatment/bigwig/chx.p.site.genome.bw \
--outFileFormat bigwig \
--normalizeUsing BPM \
--Offset 12 \
-bs 1 \
-p 32 &


done < data/"$seqtype".txt

#!/bin/bash 

star_path=../bin/STAR-2.7.8a/bin/Linux_x86_64

mkdir -p data/rnaseq_se/preprocessed_data/complete
mkdir -p data/rnaseq_se/mapped/individual

while read -ra a ;
 do
     # run cutadapt
     cutadapt  -q 30 -m 20 --report=full \
     -a AGATCGGAAGAGCACACGTCT -j 32 \
     --minimum-length 1 \
     -o data/rnaseq_se/preprocessed_data/complete/${a[0]} data/rnaseq_se/raw_data/${a[0]}

    $star_path/STAR \
       --outSAMtype BAM SortedByCoordinate \
       --runThreadN 16 \
       --outFilterMismatchNmax 2 \
       --genomeDir reference_genome/star_index_ncbi \
       --readFilesIn  data/rnaseq_se/preprocessed_data/complete/${a[0]} \
       --outFileNamePrefix data/rnaseq_se/mapped/individual/${a[1]} \
       --outFilterMultimapNmax 16 \
       --outFilterMatchNmin 16 \
       --alignEndsType EndToEnd \
       --readFilesCommand zcat \
       --outMultimapperOrder Random \
       --outSAMattributes All

     mv data/rnaseq_se/mapped/individual/${a[1]}Aligned.sortedByCoord.out.bam data/rnaseq_se/mapped/individual/${a[1]}.bam
     samtools index data/rnaseq_se/mapped/individual/${a[1]}.bam
done < data/rnaseq_se.txt
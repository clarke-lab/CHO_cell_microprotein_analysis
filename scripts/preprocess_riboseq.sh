#!/bin/bash

# num_cores=$(getconf _NPROCESSORS_ONLN)

star_path=../bin/STAR-2.7.2d/bin/Linux_x86_64

# build indexes for contamination filtering

# 1. rRNA
contam_dir=reference_genome/contaminant_sequences/
mkdir $contam_dir/star_rRNA
$star_path/STAR --runMode genomeGenerate \
--genomeSAindexNbases 6  \
--genomeDir $contam_dir/star_rRNA \
--runThreadN 32 \
--genomeFastaFiles $contam_dir/cgr_rRNA_rnacentral_v17.fasta

# 2. tRNA
mkdir $contam_dir/star_tRNA
$star_path/STAR --runMode genomeGenerate \
--genomeSAindexNbases 6  \
--genomeDir $contam_dir/star_tRNA \
--runThreadN 32 \
--genomeFastaFiles $contam_dir/cgr_tRNA_rnacentral_v17.fasta

# 2. snoRNA
mkdir $contam_dir/star_snoRNA
$star_path/STAR --runMode genomeGenerate \
--genomeSAindexNbases 6  \
--genomeDir $contam_dir/star_snoRNA \
--runThreadN 32 \
--genomeFastaFiles $contam_dir/cgr_snoRNA_rnacentral_v17.fasta

# trim and filter
for seqtype in riboseq_chx riboseq_harr riboseq_nodrug
do
  mkdir -p data/$seqtype/preproccessed_data/trimmed
  mkdir  data/$seqtype/preproccessed_data/rRNA_filter
  mkdir  data/$seqtype/preproccessed_data/tRNA_filter
  mkdir  data/$seqtype/preproccessed_data/snoRNA_filter

while read -ra a ;
  do
    # run cutadapt
    if [ $seqtype == 'riboseq_chx' ]
    then
      cutadapt  --report=full -a AGATCGGAAGAGCACACGTCT -j 32 \
      -o data/$seqtype/preproccessed_data/trimmed/${a[0]} data/$seqtype/raw_data/${a[0]}
    else
      cutadapt  --discard-untrimmed -m 20 --report=full -a AGATCGGAAGAGCACACGTCT -j 32 \
      -o data/$seqtype/preproccessed_data/trimmed/${a[0]} data/$seqtype/raw_data/${a[0]}
    fi
    # filter rRNA
    $star_path/STAR \
    --genomeLoad NoSharedMemory \
    --seedSearchStartLmaxOverLread .5 \
    --outFilterMultimapNmax 1000 \
    --outFilterMismatchNmax 2 \
    --genomeDir $contam_dir/star_rRNA \
    --runThreadN 32 \
    --outFileNamePrefix data/$seqtype/preproccessed_data/rRNA_filter/${a[1]}.rRNA. \
    --readFilesIn data/$seqtype/preproccessed_data/trimmed/${a[0]} \
    --readFilesCommand zcat \
    --outReadsUnmapped Fastx

    # filter snoRNA
    $star_path/STAR \
    --genomeLoad NoSharedMemory \
    --seedSearchStartLmaxOverLread .5 \
    --outFilterMultimapNmax 1000 \
    --outFilterMismatchNmax 2 \
    --genomeDir $contam_dir/star_snoRNA \
    --runThreadN 32 \
    --outFileNamePrefix  data/$seqtype/preproccessed_data/snoRNA_filter/${a[1]}.snoRNA. \
    --readFilesIn data/$seqtype/preproccessed_data/rRNA_filter/${a[1]}.rRNA.Unmapped.out.mate1 \
    --outReadsUnmapped Fastx

    # filter tRNA
    $star_path/STAR \
    --genomeLoad NoSharedMemory \
    --seedSearchStartLmaxOverLread .5 \
    --outFilterMultimapNmax 1000 \
    --outFilterMismatchNmax 2 \
    --genomeDir $contam_dir/star_tRNA \
    --runThreadN 32 \
    --outFileNamePrefix  data/$seqtype/preproccessed_data/tRNA_filter/${a[1]}.tRNA. \
    --readFilesIn data/$seqtype/preproccessed_data/snoRNA_filter/${a[1]}.snoRNA.Unmapped.out.mate1 \
    --outReadsUnmapped Fastx

    # # filter for 29 to 31nt
    # awk 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 29 && length(seq) <= 31) {print header, seq, qheader, qseq}}' \
    # < data/$seqtype/preproccessed_data/tRNA_filter/${a[1]}.tRNA.Unmapped.out.mate1 > \
    # data/$seqtype/preproccessed_data/${a[1]}.fastq
    #
    #
    #
    # mv data/$seqtype/preproccessed_data/tRNA_filter/${a[1]}.tRNA.Unmapped.out.mate1 data/$seqtype/preproccessed_data/${a[1]}.fastq

  done < data/"$seqtype".txt
done

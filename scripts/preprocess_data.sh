#!/bin/bash

# 1. Total RNASeq (Paired-end)

data/rnaseq_pe/preprocessed
IN_DIR=rnaseq_pe/raw_data
OUT_DIR=rnaseq_pe/preprocessed_data
mkdir $OUT_DIR
while read -ra a ; do
  cutadapt  \
  -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
  -o $OUT_DIR/"$SAMPLE_ID"_1.fastq.gz  -p $OUT_DIR/"$SAMPLE_ID"_2.fastq.gz \
  $IN_DIR/"$SAMPLE_ID"_1.fastq.gz $IN_DIR/"$SAMPLE_ID"_2.fastq.gz \
  -j 32

  java -jar  ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
  -threads $NUM_THREADS \
  $IN_DIR/"$SAMPLE_ID"_1.fastq.gz $IN_DIR/"$SAMPLE_ID"_2.fastq.gz \
  $OUT_DIR/paired/"$SAMPLE_ID"_1.fastq.gz $OUT_DIR/unpaired/"$SAMPLE_ID"_1.fastq.gz \
  $OUT_DIR/paired/"$SAMPLE_ID"_2.fastq.gz $OUT_DIR/unpaired/"$SAMPLE_ID"_2.fastq.gz \
  SLIDINGWINDOW:4:20 MINLEN:36 -trimlog $OUT_DIR/"$SAMPLE_ID".trimmomatic.log
done  < data/total_rnaseq_pe.txt

# 2. Total RNASeq (Single-end)
mkdir -p data/rnaseq_pe/preprocessed

# 3. CHX riboseq
mkdir -p data/rnaseq_pe/preprocessed

# 4. Harr riboseq
mkdir -p data/rnaseq_pe/preprocessed

# 5. No drug riboseq
mkdir -p data/rnaseq_pe/preprocessed

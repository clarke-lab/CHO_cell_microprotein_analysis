#!/bin/bash

# 1. Total RNASeq (Paired-end)
IN_DIR=data/rnaseq_pe/raw_data
OUT_DIR=data/rnaseq_pe/preprocessed_data
mkdir $OUT_DIR
mkdir $OUT_DIR/paired $OUT_DIR/unpaired
while read -ra a ; do
  file_name=${a[0]}
  sample_name=${a[1]}
   cutadapt  \
   -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
   -o $OUT_DIR/$sample_name"_1_trimmed.fastq.gz"  -p $OUT_DIR/$sample_name"_2_trimmed.fastq.gz"  \
      $IN_DIR/$file_name"_1.fastq.gz" $IN_DIR/$file_name"_2.fastq.gz" \
   -j 32

   java -jar  ../Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 32 \
    $OUT_DIR/$sample_name"_1_trimmed.fastq.gz" $OUT_DIR/$sample_name"_2_trimmed.fastq.gz" \
    $OUT_DIR/paired/$sample_name"_1.fastq.gz" $OUT_DIR/unpaired/$sample_name"_1.fastq.gz" \
    $OUT_DIR/paired/$sample_name"_2.fastq.gz" $OUT_DIR/unpaired/$sample_name"_2.fastq.gz" \
    SLIDINGWINDOW:4:20 MINLEN:36 -trimlog $OUT_DIR/$sample_name".trimmomatic.log"

    rm $OUT_DIR/$sample_name"_1_trimmed.fastq.gz" $OUT_DIR/$sample_name"_2_trimmed.fastq.gz"
    mv $OUT_DIR/paired/$sample_name"_1.fastq.gz" $OUT_DIR
    mv $OUT_DIR/paired/$sample_name"_2.fastq.gz" $OUT_DIR
 done  < data/total_rnaseq_pe_delete.txt
rm -r $OUT_DIR/paired $OUT_DIR/unpaired

# 1. Total RNASeq (Single-end)
mkdir data/rnaseq_se/preprocessed_data/
while read -ra a ;
 do
     # run cutadapt
     cutadapt  -q 30 -m 20 --report=full \
     -a AGATCGGAAGAGCACACGTCT -j 32 \
     -o data/rnaseq_se/preprocessed_data/${a[0]} data/rnaseq_se/raw_data/${a[0]}

 done < data/total_rnaseq_se.txt

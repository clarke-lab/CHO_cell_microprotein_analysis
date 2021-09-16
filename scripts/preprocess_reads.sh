#!/bin/bash
#### Description: Prepocesses raw sequencing data for   
####              1. Ribo-seq
####                  a) adapter trimming
####                  b) ncRNA contamination
####                  c) Read length based on phasing
####              2. RNA-seq  
####                  a) adapter trimming
#### 
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

# num_cores=$(getconf_NPROCESSORS_ONLN)

star_path=../bin/STAR-2.7.8a/bin/Linux_x86_64

# Ribo-seq
# build indexes for contamination filtering

# 1. rRNA
contam_dir=reference_genome/contaminant_sequences/
mkdir $contam_dir/star_rRNA
$star_path/STAR --runMode genomeGenerate \
--genomeSAindexNbases 6  \
--genomeDir $contam_dir/star_rRNA \
--runThreadN 32 \
--genomeFastaFiles $contam_dir/rna_central_v18_rrna.fasta

# 2. tRNA
$star_path/STAR --runMode genomeGenerate \
--genomeSAindexNbases 6  \
--genomeDir $contam_dir/star_tRNA \
--runThreadN 32 \
--genomeFastaFiles $contam_dir/rna_central_v18_trna.fasta

# # # 2. snoRNA
mkdir $contam_dir/star_snoRNA
$star_path/STAR --runMode genomeGenerate \
--genomeSAindexNbases 6  \
--genomeDir $contam_dir/star_snoRNA \
--runThreadN 32 \
--genomeFastaFiles $contam_dir/rna_central_v18_snrna.fasta

# trim and filter contamination from each sample 
# for each type of Ribo-seq data
for seqtype in riboseq_chx riboseq_harr riboseq_nd
do
   mkdir -p data/$seqtype/preproccessed_data/trimmed
   mkdir  data/$seqtype/preproccessed_data/rRNA_filter
   mkdir  data/$seqtype/preproccessed_data/tRNA_filter
   mkdir  data/$seqtype/preproccessed_data/snoRNA_filter
   mkdir data/$seqtype/preproccessed_data/complete
   mkdir -p data/$seqtype/mapped/individual

while read -ra a ;
  do
    # run cutadapt; note chx data has adapter removed by sequencing provider
    if [ $seqtype == 'riboseq_chx' ]
    then
      cutadapt  --report=full -a AGATCGGAAGAGCACACGTCT -j 50 --minimum-length 1 \
      -o data/$seqtype/preproccessed_data/trimmed/${a[0]} data/$seqtype/raw_data/${a[0]}
    else
      cutadapt  --discard-untrimmed -m 20 --report=full -a AGATCGGAAGAGCACACGTCT -j 50 --minimum-length 1 \
      -o data/$seqtype/preproccessed_data/trimmed/${a[0]} data/$seqtype/raw_data/${a[0]}
    fi
  
    # filter rRNA
    $star_path/STAR \
    --genomeLoad NoSharedMemory \
    --seedSearchStartLmaxOverLread .5 \
    --outFilterMultimapNmax 1000 \
    --outFilterMismatchNmax 2 \
    --genomeDir $contam_dir/star_rRNA \
    --runThreadN 50 \
    --outFileNamePrefix data/$seqtype/preproccessed_data/rRNA_filter/${a[1]}.rRNA. \
    --readFilesIn data/$seqtype/preproccessed_data/trimmed/${a[0]} \
    --readFilesCommand zcat \
    --outReadsUnmapped Fastx
    
    # # filter snoRNA
    $star_path/STAR \
    --genomeLoad NoSharedMemory \
    --seedSearchStartLmaxOverLread .5 \
    --outFilterMultimapNmax 1000 \
    --outFilterMismatchNmax 2 \
    --genomeDir $contam_dir/star_snoRNA \
    --runThreadN 50 \
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
    --runThreadN 50 \
    --outFileNamePrefix  data/$seqtype/preproccessed_data/tRNA_filter/${a[1]}.tRNA. \
    --readFilesIn data/$seqtype/preproccessed_data/snoRNA_filter/${a[1]}.snoRNA.Unmapped.out.mate1 \
    --outReadsUnmapped Fastx
  
    # remove unecessary files
    rm data/$seqtype/preproccessed_data/*_filter/${a[1]}.*.Aligned.out.sam
    rm data/$seqtype/preproccessed_data/*_filter/${a[1]}.*.Log.progress.out
    rm data/$seqtype/preproccessed_data/*_filter/${a[1]}.*.SJ.out.tab

    # filter the read lengths based on phasing qc (28nt to 31nt)
    awk 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 28 && length(seq) <= 31) {print header, seq, qheader, qseq}}' \
    < data/$seqtype/preproccessed_data/tRNA_filter/${a[1]}.tRNA.Unmapped.out.mate1 > data/$seqtype/preproccessed_data/complete/${a[1]}.fastq
  
    # map the retained reads to the PICR reference genome
    $star_path/STAR \
    --outFilterType BySJout \
    --runThreadN 16 \
    --outFilterMismatchNmax 2 \
    --genomeDir reference_genome/star_index_ncbi \
    --readFilesIn data/$seqtype/preproccessed_data/complete/${a[1]}.fastq \
    --outFileNamePrefix data/$seqtype/mapped/individual/${a[1]} \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --outFilterMultimapNmax 1 \
    --outFilterMatchNmin 16 \
    --alignEndsType EndToEnd \
    --outSAMattributes All

    # rename
    mv data/$seqtype/mapped/individual/${a[1]}"Aligned.sortedByCoord.out.bam" \
    data/$seqtype/mapped/individual/${a[1]}".bam"

    # index bam
    samtools index data/$seqtype/mapped/individual/${a[1]}".bam"

  done < data/"$seqtype".txt
done

# merge and map the files for each treatment for ORF-RATER
for seqtype in riboseq_harr riboseq_nd riboseq_chx
  do
    mkdir data/$seqtype/mapped/merged

    # merge the preprocessed reads 
    cat data/$seqtype/preproccessed_data/complete/*.fastq > \
    data/$seqtype/mapped/merged/$seqtype.fastq

    $star_path/STAR \
    --outFilterType BySJout \
    --runThreadN 16 \
    --outFilterMismatchNmax 2 \
    --genomeDir reference_genome/star_index_ncbi \
    --readFilesIn data/$seqtype/mapped/merged/$seqtype.fastq \
    --outFileNamePrefix data/$seqtype/mapped/merged/$seqtype \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --outFilterMultimapNmax 1 \
    --outFilterMatchNmin 16 \
    --alignEndsType EndToEnd \
    --outSAMattributes All

    # rename
    mv data/$seqtype/mapped/merged/$seqtype"Aligned.sortedByCoord.out.bam" \
    data/$seqtype/mapped/merged/$seqtype".bam"

    # index bam
    samtools index data/$seqtype/mapped/merged/$seqtype.bam
done

# 2. RNASeq data 
mkdir -p data/rnaseq_se/preprocessed_data/complete
mkdir -p data/rnaseq_se/mapped/individual

# individual samples
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

# merged RNA-seq 
mkdir -p data/rnaseq_se/mapped/merged 

zcat data/rnaseq_se/preprocessed_data/complete/* > \
data/rnaseq_se/mapped/merged/rnaseq_se.fastq 

$star_path/STAR \
    --outFilterType BySJout \
    --runThreadN 16 \
    --outFilterMismatchNmax 2 \
    --genomeDir reference_genome/star_index_ncbi \
    --readFilesIn data/rnaseq_se/mapped/merged/rnaseq_se.fastq \
    --outFileNamePrefix data/rnaseq_se/mapped/merged/rnaseq_se \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --outFilterMultimapNmax 1 \
    --outFilterMatchNmin 16 \
    --alignEndsType EndToEnd \
    --outSAMattributes All
# end
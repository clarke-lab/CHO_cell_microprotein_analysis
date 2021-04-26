#!/bin/bash

# create the star index
star_path=../bin/STAR-2.7.2d/bin/Linux_x86_64
mkdir -p transcriptome_assembly/star_index

$star_path/STAR --runThreadN 32 \
     --runMode genomeGenerate \
     --sjdbOverhang 124 \
     --genomeChrBinNbits 16 \
     --genomeDir transcriptome_assembly/star_index \
     --genomeFastaFiles reference_genome/cgr_ensembl_v103.fa \
     --sjdbGTFfile reference_genome/cgr_ensembl_v103.gtf

# map PE RNASeq data to the reference genome
mkdir data/rnaseq_pe/mapped
bam_dir=data/rnaseq_pe/mapped
fastq_dir=data/rnaseq_pe/preprocessed_data/

while read -ra a ; do
sample_name=${a[1]}
$star_path/STAR \
--runThreadN 30 \
--readFilesIn $fastq_dir/$sample_name"_1.fastq.gz" $fastq_dir/$sample_name"_2.fastq.gz" \
--genomeDir reference_genome/star_index_ncbi \
--readFilesCommand gunzip -c \
--outFileNamePrefix $bam_dir/$sample_name \
--outSAMtype BAM SortedByCoordinate \
--twopassMode Basic
mv $bam_dir/$sample_name"Aligned.sortedByCoord.out.bam" $bam_dir/$sample_name".bam"
samtools index $bam_dir/$sample_name".bam"
done  < data/total_rnaseq_pe.txt

mkdir -p transcriptome_assembly/individual_gtfs
stringtie_path=/mnt/HDD2/colin/bin/stringtie-2.1.5.Linux_x86_64
while read -ra a ; do
  sample_name=${a[1]}
  $stringtie_path/stringtie \
  -p 4 \
  -G reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gff \
  --rf \
  -j 5 \
  -o transcriptome_assembly/individual_gtfs/$sample_name".gtf" \
  $bam_dir/$sample_name".bam" &
done  < data/total_rnaseq_pe.txt

readlink -f transcriptome_assembly/individual_gtfs/*.gtf >> \
transcriptome_assembly/mergelist.txt

$stringtie_path/stringtie \
--merge transcriptome_assembly/mergelist.txt \
-o transcriptome_assembly/cgr_stringtie_assembly_raw.ncbi.gtf \
-G reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gff \
-f 0.1  \
-c 10

wget https://gist.githubusercontent.com/gpertea/b83f1b32435e166afa92a2d388527f4b/raw/2b6ef946f33ebfa0b7ab354fe9f118bec5f0f112/mstrg_prep.pl \
-P scripts

perl scripts/mstrg_prep.pl transcriptome_assembly/cgr_stringtie_assembly_raw.gtf > \
transcriptome_assembly/cgr_stringtie_assembly.gtf
rm scripts/mstrg_prep.pl transcriptome_assembly/cgr_stringtie_assembly_raw.gtf

# rebuild the star index with the stringtie gtf
mkdir reference_genome/star_index_stringtie
$star_path/STAR --runThreadN 32 \
     --runMode genomeGenerate \
     --sjdbOverhang 124 \
     --genomeChrBinNbits 16 \
     --genomeDir reference_genome/star_index_stringtie \
     --genomeFastaFiles reference_genome/cgr_ensembl_v103.fa \
     --sjdbGTFfile transcriptome_assembly/cgr_stringtie_assembly.gtf

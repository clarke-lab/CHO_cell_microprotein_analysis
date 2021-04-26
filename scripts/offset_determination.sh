#!/bin/bash

star_path=../bin/STAR-2.7.2d/bin/Linux_x86_64

awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' \
plastid_analysis/merged_data/riboseq_chx.fastq


source ~/.bashrc
conda activate plastid
conda activate ~/anaconda2/envs/plastid
metagene generate plastid_analysis/merlin_orfs \
                     --landmark cds_start \
                     --annotation_files reference_genome/cgr_ensembl_v103_protein_coding.gtf

for seqtype in riboseq_chx riboseq_harr riboseq_nodrug
do

  psite plastid_analysis/merlin_orfs_rois.txt plastid_analysis/$seqtype \
                          --min_length 28 \
                          --max_length 31 \
                          --require_upstream \
                          --min_counts 10 \
                          --count_files plastid_analysis/mapped/$seqtype"Aligned.sortedByCoord.out.bam"
done



metagene generate plastid_analysis/merlin_orfs_cds_start \
                    --landmark cds_start \
                    --annotation_files reference_genome/cgr_ensembl_v103_protein_coding.gtf

for seqtype in riboseq_chx riboseq_harr riboseq_nodrug
do
  phase_by_size plastid_analysis/merlin_orfs_cds_start_rois.txt plastid_analysis/$seqtype \
                  --count_files plastid_analysis/mapped/$seqtype"Aligned.sortedByCoord.out.bam" \
                  --fiveprime \
                  --offset 12 \
                  --codon_buffer 5 \
                  --min_length 28 \
                  --max_length 31

done


# plastid individual files
for seqtype in riboseq_chx riboseq_harr riboseq_nodrug
do
  mkdir -p plastid_analysis/individual_files/$seqtype/mapped
  mkdir plastid_analysis/individual_files/$seqtype/offset
  mkdir plastid_analysis/individual_files/$seqtype/periodicity

while read -ra a ;
do
$star_path/STAR \
  --outSAMtype BAM SortedByCoordinate \
  --runThreadN 32 \
  --outFilterMismatchNmax 2 \
  --genomeDir reference_genome/star_index_stringtie \
  --readFilesIn  data/$seqtype/preproccessed_data/tRNA_filter/${a[1]}".tRNA.Unmapped.out.mate1" \
  --outFileNamePrefix plastid_analysis/individual_files/$seqtype/mapped/${a[1]} \
  --outFilterMultimapNmax 16 \
  --outFilterMatchNmin 16 \
  --alignEndsType EndToEnd \
  --outMultimapperOrder Random \
  --outSAMattributes All
#
  samtools index plastid_analysis/individual_files/$seqtype/mapped/${a[1]}Aligned.sortedByCoord.out.bam

  psite plastid_analysis/merlin_orfs_rois.txt  plastid_analysis/individual_files/$seqtype/offset/${a[1]} \
      --min_length 28 \
      --max_length 31 \
      --require_upstream \
      --count_files plastid_analysis/individual_files/$seqtype/mapped/${a[1]}Aligned.sortedByCoord.out.bam

  phase_by_size plastid_analysis/merlin_orfs_cds_start_rois.txt plastid_analysis/individual_files/$seqtype/periodicity/${a[1]} \
                  --count_files plastid_analysis/individual_files/$seqtype/mapped/${a[1]}Aligned.sortedByCoord.out.bam \
                  --fiveprime  \
                  --offset 12 \
                  --codon_buffer 5 \
                  --min_length 28 \
                  --max_length 31

done < data/"$seqtype".txt
done

# filter reads for 29nt and 30nt reads
for seqtype in riboseq_chx riboseq_harr riboseq_nodrug
do
while read -ra a ;
do
 awk 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 28 && length(seq) <= 31) {print header, seq, qheader, qseq}}' \
 < data/$seqtype/preproccessed_data/tRNA_filter/${a[1]}.tRNA.Unmapped.out.mate1 > \
 data/$seqtype/preproccessed_data/${a[1]}.fastq
done < data/"$seqtype".txt
done

# zcat my.fastq.gz | echo $((`wc -l`/4))

# cat my.fastq.gz | echo $((`wc -l`/4))
# echo $(cat riboseq_chx/preproccessed_data/*.fastq |wc -l)/4| bc


  for seqtype in riboseq_chx riboseq_harr riboseq_nodrug
  do

    mkdir -p plastid_analysis/merged_files/mapped
    rm data/$seqtype/preproccessed_data/$seqtype"_merged.fastq"
    # merge each seqtype
     cat data/$seqtype/preproccessed_data/*.fastq > \
     data/$seqtype/preproccessed_data/$seqtype"_merged.fastq"

    # # read distributions
    # awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' \
    # plastid_analysis/merged_data/$seqtype.fastq > plastid_analysis/read_dist/$seqtype.txt &

    $star_path/STAR \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN 32 \
    --outFilterMismatchNmax 2 \
    --genomeDir reference_genome/star_index_stringtie \
    --readFilesIn  data/$seqtype/preproccessed_data/$seqtype"_merged.fastq" \
    --outFileNamePrefix plastid_analysis/merged_files/mapped/$seqtype \
    --outFilterMultimapNmax 16 \
    --outFilterMatchNmin 16 \
    --alignEndsType EndToEnd \
    --outMultimapperOrder Random \
    --outSAMattributes All

    samtools index plastid_analysis/merged_files/mapped/$seqtype"Aligned.sortedByCoord.out.bam"

     psite plastid_analysis/merlin_orfs_rois.txt  plastid_analysis/merged_files/$seqtype \
           --min_length 28 \
           --max_length 31 \
           --require_upstream \
           --count_files plastid_analysis/merged_files/mapped/$seqtype"Aligned.sortedByCoord.out.bam"

     phase_by_size plastid_analysis/merlin_orfs_cds_start_rois.txt plastid_analysis/merged_files/$seqtype \
                     --count_files plastid_analysis/merged_files/mapped/$seqtype"Aligned.sortedByCoord.out.bam" \
                     --fiveprime \
                     --offset 12 \
                     --codon_buffer 5 \
                     --min_length 28 \
                     --max_length 31

  done


conda deactivate plastid

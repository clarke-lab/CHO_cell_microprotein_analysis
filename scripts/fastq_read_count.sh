#!/bin/bash
#### Description: Count the reads from each stage of preprocessing for the Ribo-seq data   
####              1. Raw reads
####              2. Trimmed
####              3. rRNA, tRNA, snoRNA
####              3. Final RPFs 28-31nt   
#### 
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

# 1. count the raw reads
mkdir results/read_counts
for seqtype in riboseq_chx riboseq_harr riboseq_nd
do
    for f in data/$seqtype/raw_data/*fastq.gz; 
        do echo -n "$f | "  && echo $(($(zcat $f | echo $((`wc -l`/4))))); 
    done > results/read_counts/$seqtype.raw.counts && sed  -i '1i file | raw_read_number' results/read_counts/$seqtype.raw.counts
done

# 2. count the reads removed by trimming 
for seqtype in riboseq_chx riboseq_harr riboseq_nd
do
    for f in data/$seqtype/preproccessed_data/trimmed/ribo*fastq.gz; 
        do echo -n "$f | "  && echo $(($(zcat $f | echo $((`wc -l`/4))))); 
    done > results/read_counts/$seqtype.trimmed.counts && sed  -i '1i file | raw_read_number' results/read_counts/$seqtype.trimmed.counts
done

# 3. count the reads remaining after filtering non-coding RNA
for seqtype in riboseq_chx riboseq_harr riboseq_nd
do
    for step in rRNA_filter snoRNA_filter tRNA_filter
    do
        for f in data/"$seqtype"/preproccessed_data/"$step"/*Unmapped.out.mate1; 
            do echo -n "$f | "  && echo $(($(cat $f | echo $((`wc -l`/4))))); 
        done > results/read_counts/$seqtype.$step.counts && sed  -i '1i file | raw_read_number' results/read_counts/$seqtype.$step.counts
    done
done

# 4. count the reads removed after preprocessing
for seqtype in riboseq_chx riboseq_harr riboseq_nd
do
    for f in data/$seqtype/preproccessed_data/complete/*fastq; 
        do echo -n "$f | "  && echo $(($(cat $f | echo $((`wc -l`/4))))); 
    done > results/read_counts/$seqtype.final.counts && sed  -i '1i file | raw_read_number' results/read_counts/$seqtype.final.counts
done







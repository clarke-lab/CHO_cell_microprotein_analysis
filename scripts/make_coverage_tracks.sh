#!/bin/bash
#### Description: Make coverage tracks for the manuscript figures
####             1. Merged Data for Figure 2B
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie 

### 1 merged transcriptome tracks for new annotations 


mkdir -p alignment_tracks/merged
bg_dir=alignment_tracks/merged

# 1. Merged Data for Figure 2B
## Demonstrate the coverage at transcript-level
## a) CHX full coverage
## b) CHX p-site offset
## c) CHX a-site offset
## d) HARR-ND p-site

## sort and index transcriptome aligned BAMs
for seqtype in riboseq_harr riboseq_chx riboseq_nd
do
    samtools sort \
    data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.bam" \
    -o data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam"

    samtools index \
    data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam"
done

## a) CHX full coverage
chx_trancript_bam_sorted=data/riboseq_chx/mapped/merged/riboseq_chxAligned.toTranscriptome.out.sort.bam

bamCoverage -b $chx_trancript_bam_sorted -o $bg_dir/chx.fullcov.site.transcriptome_bs0.bedgraph \
 --outFileFormat bedgraph -bs 1 -p 70 --normalizeUsing BPM  --smoothLength 25

## b) CHX p-site offset
mkdir -p alignment_tracks/merged
bg_dir=alignment_tracks/merged

## bug in deeptools: offset needs to be +1 actual (12)
bamCoverage \
-b data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam" \
-o alignment_tracks/merged/$seqtype".p.site.transcriptome.bedgraph" \
--outFileFormat bedgraph --Offset 13 -bs 1 -p 70 --normalizeUsing BPM 

## c) CHX a-site offset
bamCoverage \
-b data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam" \
-o alignment_tracks/merged/$seqtype".a.site.transcriptome.bedgraph" \
--outFileFormat bedgraph --Offset 16 -bs 1 -p 70 --normalizeUsing BPM 

## d) HARR-ND p-site offset
# make bigwigs to do the substraction
for seqtype in riboseq_harr riboseq_nd
do
    bamCoverage \
    -b data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam" \
    -o alignment_tracks/merged/$seqtype".p.site.transcriptome.bigwig" \
    --outFileFormat bigwig --Offset 13 -bs 1 -p 70 --normalizeUsing BPM 
done

bigwigCompare \
--bigwig1 alignment_tracks/merged/riboseq_harr.p.site.transcriptome.bigwig \
--bigwig2 alignment_tracks/merged/riboseq_nd.p.site.transcriptome.bigwig \
--outFileName alignment_tracks/merged/harr-nd.bigwig \
--outFileFormat bigwig \
--numberOfProcessors 70 \
-bs 1 \
--skipNonCoveredRegions \
--skipZeroOverZero \
--operation subtract

# convert to bedgraph for R
../bin/kentUtils/bin/linux.x86_64/bigWigToBedGraph \
alignment_tracks/merged/harr-nd.bigwig \
alignment_tracks/merged/harr-nd.bedgraph


for seqtype in riboseq_harr riboseq_chx riboseq_nd
do
    samtools sort \
    data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.bam" \
    -o data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam"

    samtools index \
    data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam"

bamCoverage \
-b data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam" \
-o alignment_tracks/merged/$seqtype".p.site.transcriptome.bedgraph" \
--outFileFormat bedgraph --Offset 13 -bs 1 -p 70 --normalizeUsing BPM 

bamCoverage \
-b data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam" \
-o alignment_tracks/merged/$seqtype".a.site.transcriptome.bedgraph" \
--outFileFormat bedgraph --Offset 16 -bs 1 -p 70 --normalizeUsing BPM 
done



chx_trancript_bam=data/riboseq_chx/mapped/merged/riboseq_chxAligned.toTranscriptome.out.bam
chx_trancript_bam_sorted=data/riboseq_chx/mapped/merged/riboseq_chxAligned.toTranscriptome.out.sort.bam

samtools sort $chx_trancript_bam -o $chx_trancript_bam_sorted
samtools index $chx_trancript_bam_sorted
bamCoverage -b $chx_trancript_bam_sorted -o $bg_dir/chx.fullcov.site.transcriptome.bedgraph \
 --outFileFormat bedgraph -bs 1 -p 70 --normalizeUsing BPM

bamCoverage -b $chx_trancript_bam_sorted -o $bg_dir/chx.fullcov.site.transcriptome_bs2.bedgraph \
 --outFileFormat bedgraph -bs 2 -p 70 --normalizeUsing BPM
bamCoverage -b $chx_trancript_bam_sorted -o $bg_dir/chx.fullcov.site.transcriptome_bs0.bedgraph \
 --outFileFormat bedgraph -bs 1 -p 70 --normalizeUsing BPM  --smoothLength 25

# mkdir -p alignment_tracks/individual/riboseq_chx/
# mkdir -p alignment_tracks/individual/rnaseq/

# for seqtype in riboseq_chx 
# do
#     for replicate in  nts_r1 nts_r2 nts_r3 nts_r4 ts_r1 ts_r2 ts_r3 ts_r4
#     do
#         # samtools sort data/$seqtype/mapped/individual/$replicate"Aligned.toTranscriptome.out.bam" \
#         # -o data/$seqtype/mapped/individual/$replicate"Aligned.toTranscriptome.out.sorted.bam"
        
#         samtools index data/$seqtype/mapped/individual/$replicate"Aligned.toTranscriptome.out.sorted.bam"

        
        
#         bamCoverage \
#         -b data/$seqtype/mapped/individual/$replicate"Aligned.toTranscriptome.out.sorted.bam" \
#         -o alignment_tracks/individual/$seqtype/$replicate".coverage.bigwig" \
#         --outFileFormat bigwig -bs 1 -p 32 --normalizeUsing BPM
#     done
# done 


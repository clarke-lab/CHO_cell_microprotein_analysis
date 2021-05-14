mkdir -p alignment_tracks/merged

bg_dir=alignment_tracks/merged


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
--numberOfProcessors 60 \
-bs 1 \
--skipNonCoveredRegions \
--skipZeroOverZero \
--operation subtract


../bin/kentUtils/bin/linux.x86_64/bigWigToBedGraph \
alignment_tracks/merged/harr-nd.bigwig \
alignment_tracks/merged/harr-nd.bedgraph


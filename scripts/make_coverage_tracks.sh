#!/bin/bash
#### Description: Make coverage tracks for the manuscript figures
####             1. Merged Data for Figure 2B
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie 

### 1 merged transcriptome tracks for new annotations 




# # 1. Merged Data for Figure 2B
# ## Demonstrate the coverage at transcript-level
# ## a) CHX full coverage
# ## b) CHX p-site offset
# ## c) CHX a-site offset
# ## d) HARR-ND p-site

mkdir -p alignment_tracks/merged
merged_dir=alignment_tracks/merged

## sort and index transcriptome aligned BAMs
for seqtype in riboseq_harr riboseq_chx riboseq_nd
do
    samtools sort \
    data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.bam" \
    -o data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam"

    samtools index \
    data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam"
done

 ## a) CHX
chx_trancript_bam_sorted=data/riboseq_chx/mapped/merged/riboseq_chxAligned.toTranscriptome.out.sort.bam

# full coverage
bamCoverage -b $chx_trancript_bam_sorted -o $merged_dir/chx.fullcov.transcriptome.bedgraph \
  --outFileFormat bedgraph -bs 1 -p 70 --normalizeUsing BPM --smoothLength 25

## p-site
bamCoverage -b $chx_trancript_bam_sorted -o $merged_dir/chx.p.site.transcriptome.bedgraph \
  --outFileFormat bedgraph -bs 1 -p 70 --normalizeUsing BPM  --Offset 13

 ## d) HARR-ND p-site offset
 # make bigwigs to do the substraction
 for seqtype in riboseq_harr riboseq_nd
 do
     bamCoverage \
     -b data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam" \
     -o $merged_dir/$seqtype".p.site.transcriptome.bigwig" \
     --outFileFormat bigwig --Offset 13 -bs 1 -p 70 --normalizeUsing BPM 
done

 bigwigCompare \
 --bigwig1 $merged_dir/riboseq_harr.p.site.transcriptome.bigwig \
 --bigwig2 $merged_dir/riboseq_nd.p.site.transcriptome.bigwig \
 --outFileName $merged_dir/harr-nd.bigwig \
 --outFileFormat bigwig \
 --numberOfProcessors 70 \
 -bs 1 \
 --skipNonCoveredRegions \
 --skipZeroOverZero \
 --operation subtract

# convert to bedgraph for R
 ../bin/kentUtils/bin/linux.x86_64/bigWigToBedGraph \
 $merged_dir/harr-nd.bigwig \
 $merged_dir/harr-nd.bedgraph

# ## 2. Individual genome coverage tracks for 
# ##  a) RNA-seq
# ##  b) CHX Ribo-seq Full coverage
# ##  c) CHX Ribo-seq A-site coverage
# ## Note: Track scaled by DESeq2 size factor
# ## Note: Deep tools bug for offset = psite + 3=Asite=15+1(to correct bug)

mkdir -p alignment_tracks/individual/riboseq/
mkdir -p alignment_tracks/individual/rnaseq/

rnaseq_indiv_dir=alignment_tracks/individual/rnaseq/
riboseq_indiv_dir=alignment_tracks/individual/riboseq/

# set sample names
sample_names=("nts_r1" "nts_r2" "nts_r3" "nts_r4" "ts_r1" "ts_r2" "ts_r3" "ts_r4")

# DESeq2 size factors
readarray -t rna_scale <results/section2.4/rnaseq_size_factors.txt
readarray -t ribo_scale <results/section2.4/rnaseq_size_factors.txt

for i in "${!sample_names[@]}";
do

# No norm, no offset, scaled by DESeq2
    bamCoverage -b data/rnaseq_se/mapped/individual/${sample_names[i]}.bam \
    -o $rnaseq_indiv_dir/"${sample_names[i]}".rnaseq_deseq.bw \
    --outFileFormat bigwig -bs 1 -p 70 --normalizeUsing None  --skipNonCoveredRegions --scaleFactor ${rna_scale[i]}

    bamCoverage -b data/riboseq_chx/mapped/individual/${sample_names[i]}.bam \
    -o $riboseq_indiv_dir/"${sample_names[i]}".chx_deseq.bw \
    --outFileFormat bigwig -bs 1 -p 70 --normalizeUsing None --skipNonCoveredRegions --scaleFactor ${ribo_scale[i]}

# DESeq A-site offset for riboseq
    
    bamCoverage -b data/riboseq_chx/mapped/individual/${sample_names[i]}.bam \
    -o $riboseq_indiv_dir/"${sample_names[i]}".chx_deseq_asite.bw \
    --outFileFormat bigwig -bs 1 -p 70 --normalizeUsing None --Offset 16 --skipNonCoveredRegions --scaleFactor ${ribo_scale[i]}

done
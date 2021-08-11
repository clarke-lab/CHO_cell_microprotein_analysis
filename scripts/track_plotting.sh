#!/bin/bash

mkdir -p alignment_tracks/merged


bg_dir=alignment_tracks/merged


# chx_trancript_bam=data/riboseq_chx/mapped/merged/riboseq_chxAligned.toTranscriptome.out.bam
# chx_trancript_bam_sorted=data/riboseq_chx/mapped/merged/riboseq_chxAligned.toTranscriptome.out.sort.bam

# samtools sort $chx_trancript_bam -o $chx_trancript_bam_sorted
# samtools index $chx_trancript_bam_sorted
# bamCoverage -b $chx_trancript_bam_sorted -o $bg_dir/chx.fullcov.site.transcriptome.bedgraph \
# --outFileFormat bedgraph -bs 1 -p 70 --normalizeUsing BPM



# harr_trancript_bam=data/riboseq_harr/mapped/merged/riboseq_harrAligned.toTranscriptome.out.bam
# harr_trancript_bam_sorted=data/riboseq_harr/mapped/merged/riboseq_harrAligned.toTranscriptome.out.sort.bam

# samtools sort $harr_trancript_bam -o $harr_trancript_bam_sorted
# samtools index $harr_trancript_bam_sorted
# bamCoverage -b $harr_trancript_bam_sorted -o $bg_dir/harr.fullcov.site.transcriptome.bedgraph \
# --outFileFormat bedgraph -bs 1 -p 70 --normalizeUsing BPM


# nd_trancript_bam=data/riboseq_nd/mapped/merged/riboseq_ndAligned.toTranscriptome.out.bam
# nd_trancript_bam_sorted=data/riboseq_nd/mapped/merged/riboseq_ndAligned.toTranscriptome.out.sort.bam

# samtools sort $nd_trancript_bam -o $nd_trancript_bam_sorted
# samtools index $nd_trancript_bam_sorted
# bamCoverage -b $nd_trancript_bam_sorted -o $bg_dir/nd.fullcov.site.transcriptome.bedgraph \
# --outFileFormat bedgraph -bs 1 -p 70 --normalizeUsing BPM


# --scaleFactor "${rpf_scale[i]}"  --scaleFactor "${rna_scale[i]}"

 sample_names=("nts_r1" "nts_r2" "nts_r3" "nts_r4" "ts_r1" "ts_r2" "ts_r3" "ts_r4")

 rpf_scale=("0.8754236" "1.0244444" "1.5854728" "1.2567031" "1.2301475" "1.0193872" "2.5516824" "1.5839588")

 rna_scale=("0.7477485" "0.6223411" "0.6120025" "0.6309549" "1.0389609" "0.6670906" "0.9084043" "0.8499545")


track_dir=alignment_tracks/individual/diff_trans
mkdir $track_dir

for i in "${!sample_names[@]}";
do
# Bins per million - no offset

    bamCoverage -b data/rnaseq_se/mapped/individual/${sample_names[i]}.bam \
    -o $track_dir/"${sample_names[i]}".rnaseq_bpm.bw \
    --outFileFormat bigwig -bs 1 -p 70 --normalizeUsing BPM --skipNonCoveredRegions 

    bamCoverage -b data/riboseq_chx/mapped/individual/${sample_names[i]}.bam \
    -o $track_dir/"${sample_names[i]}".chx_bpm.bw \
    --outFileFormat bigwig -bs 1 -p 70 --normalizeUsing BPM --skipNonCoveredRegions

# BPM A-site offset for riboseq
    
    bamCoverage -b data/riboseq_chx/mapped/individual/${sample_names[i]}.bam \
    -o $track_dir/"${sample_names[i]}".chx_bpm_offset.bw \
    --outFileFormat bigwig -bs 1 -p 70 --normalizeUsing BPM --Offset 16 --skipNonCoveredRegions

# No norm, no offset, scaled by DESeq2
    bamCoverage -b data/rnaseq_se/mapped/individual/${sample_names[i]}.bam \
    -o $track_dir/"${sample_names[i]}".rnaseq_deseq.bw \
    --outFileFormat bigwig -bs 1 -p 70 --normalizeUsing None  --skipNonCoveredRegions --scaleFactor ${rna_scale[i]}

    bamCoverage -b data/riboseq_chx/mapped/individual/${sample_names[i]}.bam \
    -o $track_dir/"${sample_names[i]}".chx_deseq.bw \
    --outFileFormat bigwig -bs 1 -p 70 --normalizeUsing None --skipNonCoveredRegions --scaleFactor ${rpf_scale[i]}

# DESeq A-site offset for riboseq
    
    bamCoverage -b data/riboseq_chx/mapped/individual/${sample_names[i]}.bam \
    -o $track_dir/"${sample_names[i]}".chx_deseq_offset.bw \
    --outFileFormat bigwig -bs 1 -p 70 --normalizeUsing None --Offset 16 --skipNonCoveredRegions --scaleFactor ${rpf_scale[i]}


    # ribo_bigwig=$track_dir/"${sample_names[i]}".chx.bw
    # rna_bigwig=$track_dir/"${sample_names[i]}".rnaseq.bw

    # bigwigCompare \
    # --bigwig1 $ribo_bigwig \
    # --bigwig2 $rna_bigwig \
    # --operation ratio \
    # -bs 1 \
    # -p32 \
    # --outFileFormat bigwig \
    # -o $track_dir/"${sample_names[i]}".rpf_density.bw   

done
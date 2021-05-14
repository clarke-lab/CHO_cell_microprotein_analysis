mkdir -p alignment_tracks/merged

bg_dir=alignment_tracks/merged

chx_bam=data/riboseq_chx/mapped/merged/riboseq_chx.bam
nd_bam=data/riboseq_nd/mapped/merged/riboseq_nd.bam
harr_bam=data/riboseq_harr/mapped/merged/riboseq_harr.bam

harr_transcript_bam=data/riboseq_harr/mapped/merged/riboseq_harrAligned.toTranscriptome.out.sorted.bam
bamCoverage -b $harr_transcript_bam -o $bg_dir/harr.p.site.transcriptome.bedgraph \
--outFileFormat bedgraph --Offset 13 -bs 1 -p 70 --normalizeUsing BPM

bamCoverage -b $nd_bam -o $bg_dir/nd.p.site.bigwig \
--outFileFormat bigwig --Offset 13 -bs 1 -p 60 --normalizeUsing BPM

bigwigCompare --bigwig1 $bg_dir/harr.p.site.bigwig --bigwig2 $bg_dir/nd.p.site.bigwig \
--outFileName $bg_dir/nd.harr.p.site.compared.bigwig --outFileFormat bigwig \
--numberOfProcessors 60 -bs 1 --skipNonCoveredRegions --skipZeroOverZero --operation subtract

bamCoverage -b $chx_bam -o $bg_dir/chx.p.site.bigwig \
--outFileFormat bigwig --Offset 13 -bs 1 -p 60 --normalizeUsing BPM

svist4get \
-fa reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna \
-gtf reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf  \
-bg $bg_dir/harr.p.site.bedgraph \
-bl 'Harr Ribo-Seq (P-site)' \
-g Cert1_1 \
-hi \
-bul max \
-o test.png

kent_path=../bin/kentUtils/bin/linux.x86_64

bg_dir=publication_results/track_plotting/bedgraph/merged

chx_bam=plastid_analysis/merged_files/mapped/riboseq_chx.ncbi.Aligned.sortedByCoord.out.bam
nodrug_bam=plastid_analysis/merged_files/mapped/riboseq_nodrug.ncbi.Aligned.sortedByCoord.out.bam
harr_bam=plastid_analysis/merged_files/mapped/riboseq_harr.ncbi.Aligned.sortedByCoord.out.bam

bamCoverage -b $harr_bam -o $bg_dir/harr.p.site.bedgraph \
--outFileFormat bedgraph --Offset 13 -bs 1 -p 32 --normalizeUsing BPM
bamCoverage -b $chx_bam -o $bg_dir/chx.p.site.bedgraph \
--outFileFormat bedgraph --Offset 13 -bs 1 -p 32 --normalizeUsing BPM
bamCoverage -b $nodrug_bam -o $bg_dir/nodrug.p.site.bedgraph \
--outFileFormat bedgraph --Offset 13 -bs 1 -p 32 --normalizeUsing BPM

bamCoverage -b $chx_bam -o $bg_dir/chx.a.site.bg \
--outFileFormat bedgraph --Offset 16 -bs 1 -p 32 --normalizeUsing BPM
bamCoverage -b $nodrug_bam -o $bg_dir/nodrug.a.site.bg \
--outFileFormat bedgraph --Offset 16 -bs 1 -p 32 --normalizeUsing BPM

# merged a site
bamCoverage -b data/ -b $chx_bam -o $bg_dir/chx.fullcov.bg \
--outFileFormat bedgraph -bs 1 -p 32 --normalizeUsing BPM
bamCoverage -b $nodrug_bam -o $bg_dir/nodrug.fullcov.bg \
--outFileFormat bedgraph -bs 1 -p 32 --normalizeUsing BPM

for seqtype in riboseq_harr riboseq_nodrug riboseq_chx
do
mkdir publication_results/track_plotting/bedgraph/individual_psite/$seqtype
while read -ra a ;
do
  samtools index plastid_analysis/individual_files/$seqtype/mapped/${a[1]}"Aligned.sortedByCoord.out.bam"
  bamCoverage -b plastid_analysis/individual_files/$seqtype/mapped/${a[1]}"Aligned.sortedByCoord.out.bam" \
  -o publication_results/track_plotting/bedgraph/individual_psite/$seqtype/${a[1]}.bedgraph \
  --outFileFormat bedgraph -bs 1 -p 32 --normalizeUsing BPM --Offset 13

  $kent_path/bedGraphToBigWig publication_results/track_plotting/bedgraph/individual_psite/$seqtype/${a[1]}.bedgraph \
  reference_genome/chrom.sizes \
  publication_results/track_plotting/bedgraph/individual_psite/$seqtype/${a[1]}.bigwig

done < data/"$seqtype".txt
done


for sample in nts_r1 nts_r2 nts_r3 nts_r4 ts_r1 ts_r2 ts_r3 ts_r4
do
bigwigCompare \
--bigwig1 publication_results/track_plotting/bedgraph/individual_psite/riboseq_harr/backup_copy/$sample.bigwig \
--bigwig2 publication_results/track_plotting/bedgraph/individual_psite/riboseq_nodrug/backup_copy/$sample.bigwig \
--outFileName orfrater_analysis/truncation_filtering/$sample.bigwig \
--skipNonCoveredRegions \
--outFileFormat bigwig \
--operation subtract \
--binSize 1 \
--numberOfProcessors 32
done

bigwigCompare -b1 publication_results/track_plotting/differential_translation/individual/rpf/${a[3]}.bigwig \
-b2 publication_results/track_plotting/differential_translation/individual/rna/${a[3]}.bigwig \



# compare
while read -ra a ;
do
echo ${a[3]}
bigwigCompare -b1 publication_results/track_plotting/differential_translation/individual/rpf/${a[3]}.bigwig \
-b2 publication_results/track_plotting/differential_translation/individual/rna/${a[3]}.bigwig \
--skipZeroOverZero --operation first --skipNAs -bs 1 -p 32 --outFileFormat bedgraph \
-o $rpf_dir/${a[3]}.te.bedgraph
done < publication_results/te_analysis_tracks.txt



svist4get \
-gtf reference_genome/cgr_ensembl_v103.gtf \
-fa reference_genome/cgr_ensembl_v103.fa \
-bg $bg_dir/harr_psite_riborf.bg $bg_dir/chx.fullcov.site.bg $bg_dir/nodrug.a.site.bg \
-bl 'Harr Ribo-Seq (P-site)' 'CHX Ribo-Seq (A-Site)' 'Nodrug Ribo-Seq (A-Site)' \
-g ENSCGRG00015023911 \
-hi \
-bul max

    svist4get \
  -gtf reference_genome/cgr_ensembl_v103.gtf \
  -fa reference_genome/cgr_ensembl_v103.fa \
  -bg  $bg_dir/harr.p.site.bg $bg_dir/chx.p.site.bg $bg_dir/nodrug.p.site.bg  \
  -bl 'P-Site Harr Ribo-Seq' 'P-Site CHX Ribo-Seq' 'P-Site no drug Ribo-Seq' \
  -t ENSCGRT00015035352 \
  -w tis 100 30
-it 'ribosomal RNA-processing protein 7 homolog A'



  svist4get \
-gtf transcriptome_assembly/cgr_stringtie_assembly.gtf \
-fa reference_genome/cgr_ensembl_v103.fa \
-bg  $bg_dir/harr.p.site.bg  $bg_dir/harr.p.site.off.bg $bg_dir/harr.p.site.off11.bg \
-bl 'Harr Ribo-Seq OffSet=12' 'Harr Ribo-Seq OffSet=13' 'Harr Ribo-Seq OffSet=11' \
-t MSTRG.20297.1

samtools view -h -o out.sam plastid_analysis/merged_files/mapped/riboseq_harrAligned.sortedByCoord.out.bam > $bg_dir/harr_full.sam
samtools view -bS $bg_dir/harr_psite.sam > $bg_dir/harr_psite.bam
samtools index $bg_dir/harr_psite.bam
bamCoverage -b $bg_dir/harr_psite_corrected_sorted.bam -o $bg_dir/harr_psite_riborf.bg \
--outFileFormat bedgraph -bs 1 -p 32 --normalizeUsing BPM

svist4get \
-gtf reference_genome/cgr_ensembl_v103.gtf \
-fa reference_genome/cgr_ensembl_v103.fa \
-bg  $bg_dir/harr.p.site.bg  $bg_dir/harr.p.site.off.bg $bg_dir/harr_psite_riborf.bg \
-bl 'Harr Ribo-Seq OffSet=12' 'Harr Ribo-Seq OffSet=13' 'RibOrf offset 12' \
-g MSTRG.11212.1 \
-w tis 100 30

gene_name='MSTRG.5352|ENSCGRG00015015354'
svist4get \
-gtf transcriptome_assembly/cgr_stringtie_assembly.gtf \
-fa reference_genome/cgr_ensembl_v103.fa \
-bg $bg_dir/harr_psite_riborf.bg $bg_dir/chx.p.site.bg $bg_dir/nodrug.p.site.bg \
-bl 'Harr Ribo-Seq OffSet=12' 'CHX Ribo-Seq' 'Nodrug Ribo-Seq' \
-w RAZU01000061.1 4213345 4213457 \
-rc \
-bul

svist4get \
-gtf transcriptome_assembly/cgr_stringtie_assembly.gtf \
-fa reference_genome/cgr_ensembl_v103.fa \
-bg $bg_dir/harr_psite_riborf.bg $bg_dir/chx.p.site.bg $bg_dir/nodrug.p.site.bg \
-bl 'Harr Ribo-Seq OffSet=12' 'CHX Ribo-Seq' 'Nodrug Ribo-Seq' \
-t MSTRG.4534.2 \
-bul max

# make an individaul bigWig for each individal rnaseq and rpf bam file
# use the DESeq size factor to normalise
# predefinted file for each samle "indvidual TE bedgraphs"

while read -ra a ;
do
  echo ${a[0]}
  echo ${a[1]}
  echo ${a[2]}
  echo ${a[3]}
rpf_bam=data/riboseq_chx/mapped/individual/${a[0]}
rna_bam=data/rnaseq_se/mapped/individual/${a[0]}

rpf_outdir=publication_results/track_plotting/differential_translation/individual/rpf/${a[3]}.bigWig
rna_outdir=publication_results/track_plotting/differential_translation/individual/rna/${a[3]}.bigWig

# samtools index $rpf_bam &
# samtools index $rna_bam

 bamCoverage -b $rpf_bam -o $rpf_outdir \
 --outFileFormat bigwig -bs 1 -p 10 --scaleFactor ${a[1]}

 bamCoverage -b $rna_bam -o $rna_outdir \
 --outFileFormat bigwig -bs 1 -p 10 --scaleFactor ${a[2]}

done < publication_results/te_analysis_tracks.txt

rpf_dir=publication_results/track_plotting/differential_translation/individual/rpf
rna_dir=publication_results/track_plotting/differential_translation/individual/rna






kent_path=../bin/kentUtils/bin/linux.x86_64/

$kent_path/bigWigToBedGraph $rpf_outdir/nts_r1Aligned.sortedByCoord.out.bigwig $rpf_outdir/nts_r1.bg


svist4get \
-gtf reference_genome/cgr_ensembl_v103.gtf \
-fa reference_genome/cgr_ensembl_v103.fa \
-bg $rpf_dir/nts_r1.te.bedgraph \
-bl 'Harr Ribo-Seq OffSet=12' \
-w RAZU01000234.1 23146955 23272578 \
-bul max

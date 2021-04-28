#!/bin/bash
export PATH=$PATH:/mnt/HDD2/colin/bin/ORF-RATER

# merge and align riboseq datasets for ORF-RATER
star_path=../bin/STAR-2.7.2d/bin/Linux_x86_64

mkdir orfrater_analysis
# take the stringtie assembly and include annotated CDS/start/stop annotation
# from ENSEMBL

grep -E "CDS|start_codon|stop_codon" reference_genome/cgr_ensembl_v103.gtf > \
reference_genome/cgr_ensembl_v103.cds.annotations.gtf

cat transcriptome_assembly/cgr_stringtie_assembly.gtf reference_genome/cgr_ensembl_v103.cds.annotations.gtf > \
orfrater_analysis/cgr_stringtie_assembly_cds.tmp.gtf

# remove transcripts from 1 gene causing errors
grep -v ENSCGRT00015032876 completed_annotation_final.gtf | \
grep -v ENSCGRT00015032832 | \
grep -v ENSCGRT00015032897 | \
grep -v ENSCGRT00015032946 > cgr_combined_annotation.clean.gtf

gtf=cgr_combined_annotation.clean.gtf
fasta=../reference_genome/cgr_ensembl_v103.fa

/mnt/HDD2/colin/bin/kentUtils/bin/linux.x86_64/gtfToGenePred \
-ignoreGroupsWithoutExons $gtf stdout | \
/mnt/HDD2/colin/bin/kentUtils/bin/linux.x86_64/genePredToBed stdin \
cgr.orfrater.annotation.bed

grep pseudogene reference_genome/cgr_ensembl_v103.gtf | \
awk -F'\'t '{print $9}' | awk -F';' '{print $1}' | sed s/gene_id//g | \
sed s/\"//g | sed 's/[[:space:]]//g'| uniq > orfrater_analysis/pseudogene_gene_ids

grep -f orfrater_analysis/pseudogene_gene_ids orfrater_analysis/cgr_stringtie_assembly_cds.gtf | \
awk -F'\'t '{print $9}' | awk -F';' '{print $2}' | sed s/transcript_id//g | \
sed s/\"//g | sed 's/[[:space:]]//g'| uniq > orfrater_analysis/pseudogene_transcript_ids

# # modify stringtie
# grep -P "0\t0\t0" merged_bam_files/no_drugReadsPerGene.out.tab | awk '{print $1}' > reference_genome/no_count_genes
# grep -Fvf reference_genome/no_count_genes reference_genome/Cricetulus_griseus_picr.CriGri-PICR.103.gtf > reference_genome/modified1.gtf
#
# # bed file for orfrater
# gtf=reference_genome/Cricetulus_griseus_picr.CriGri-PICR.103.gtf
# fasta=reference_genome/Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.fa
#
# grep -v RAZU01001824.1 reference_genome/stringtie.gtf  > reference_genome/modified.gtf
# gtf=reference_genome/modified.gtf

chx_bam=../plastid_analysis/merged_files/mapped/riboseq_chxAligned.sortedByCoord.out.bam
nodrug_bam=../plastid_analysis/merged_files/mapped/riboseq_nodrugAligned.sortedByCoord.out.bam
harr_bam=../plastid_analysis/merged_files/mapped/riboseq_harrAligned.sortedByCoord.out.bam

prune_transcripts.py \
--inbed  cgr.orfrater.annotation.bed \
--summarytable tid_removal_summary.txt \
-p 32 \
--minlen 28 \
--maxlen 31 \
-v $fasta \
$chx_bam \
$nodrug_bam \
--pseudogenes pseudogene_transcript_ids \
--outbed transcripts.bed \
--force > 1prune.log

make_tfams.py -v --force --inbed transcripts.bed \
--tfamstem tfams >   2make_tfams.log

# potential start codon list from Kearse et al. Genes Dev 2017
rm orf.h5 # need to remove this as problem with indexing if exists
find_orfs_and_types.py $fasta --codons NTG -p 32 -v --force \
--tfamstem tfams --inbed transcripts.bed \
--orfstore orf.h5 > 3find_ORF.log

psite_trimmed.py $chx_bam \
--minrdlen 28 \
--maxrdlen 31 \
--subdir chx \
--tallyfile tallies.txt \
--cdsbed cgr.orfrater.annotation.bed \
-p 25 \
-v \
--force > psite_chx.log

psite_trimmed.py $nodrug_bam \
--minrdlen 28 \
--maxrdlen 31 \
--subdir no_drug \
--tallyfile tallies.txt \
--cdsbed cgr.orfrater.annotation.bed \
-p 32 \
-v \
--force > psite_nogrug.log

psite_trimmed.py $harr_bam \
--minrdlen 28 \
--maxrdlen 31 \
--subdir harr \
--tallyfile tallies.txt \
--cdsbed cgr.orfrater.annotation.bed \
-p 32 \
-v \
--force > psite_harr.log


regress_orfs.py $harr_bam \
--startonly \
--subdir harr \
--orfstore orf.h5 \
--inbed transcripts.bed \
--regressfile regression.h5 \
-p 32 \
-v \
--startcount 1 \
--force > regress_start.log

regress_orfs.py $chx_bam  \
--subdir chx \
--orfstore orf.h5 \
--inbed transcripts.bed \
--restrictbystarts harr \
--startcount 1 \
-p 32 \
-v \
--force > chx_regress_stop.log

regress_orfs.py $nodrug_bam \
--subdir no_drug \
--orfstore orf.h5 \
--inbed transcripts.bed \
--restrictbystarts harr \
--startcount 1 \
-p 32 \
-v \
--force > nodrug_regress_stop.log

rate_regression_output.py \
harr \
chx \
no_drug \
--orfstore orf.h5 \
--ratingsfile orfratings.h5 \
-p 32 \
--CSV rate_regression_increased_start_codons.csv \
-v \
--force > rate.log

make_orf_bed.py --minlen 10 --outbed orfrater_predictions.bed

/mnt/HDD2/colin/bin/kentUtils/bin/linux.x86_64/bedToGenePred orfrater_predictions.bed stdout | \
/mnt/HDD2/colin/bin/kentUtils/bin/linux.x86_64/genePredToGtf file stdin orfrater_predictions.gtf


# make a nucleotide fasta for the CDS gregion - from the orf-rater
cat orfrater_predictions.bed | awk -v type=CDS -f ../scripts/bed12toAnnotation.awk > orfrater_predictions.bed12

bedtools getfasta -fi ../reference_genome/cgr_ensembl_v103.fa -bed orfrater_predictions.bed12 -fo orfrater_predictions.nucleotide.fasta

# convert to protein fasta
transeq -seq orfrater_predictions.nucleotide.fasta -frame 1 -outseq orfrater_predictions.protein.fasta

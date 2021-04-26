#!/bin/bash

 grep '; pseudo'  ../reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf | \
 awk -F'\'t '{print $9}' | awk -F';' '{print $1}' | sed s/gene_id//g | \
 sed s/\"//g | sed 's/[[:space:]]//g'| sort | uniq > pseudogene_gene_ids

   grep -f pseudogene_gene_ids cgr_stringtie_assembly_cds.gtf  | \
 awk -F'\'t '{print $9}' | awk -F';' '{print $1}' | sed s/gene_id//g | \
 sed s/\"//g | sed 's/[[:space:]]//g'| sort | uniq | wc -l


 grep -f pseudogene_gene_ids ../reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf | \
 awk -F'\'t '{print $9}' | awk -F';' '{print $2}' | sed s/transcript_id//g | \
 sed s/\"//g | sed 's/[[:space:]]//g'| uniq > pseudogene_transcript_ids


export PATH=$PATH:/mnt/HDD2/colin/bin/ORF-RATER

# merge and align riboseq datasets for ORF-RATER
star_path=../bin/STAR-2.7.2d/bin/Linux_x86_64

mkdir orfrater_analysis
# take the stringtie assembly and include annotated CDS/start/stop annotation
# from ENSEMBL

grep -E "CDS|start_codon|stop_codon" ../reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf > \
ncbi.cds.annotations.gtf

cat cgr_stringtie_assembly_raw.ncbi.gtf ncbi.cds.annotations.gtf > \
cgr_stringtie_assembly_cds.gtf

sed -i s/rna-//g cgr_stringtie_assembly_cds.gtf
sed -i s/gene-//g cgr_stringtie_assembly_cds.gtf
sed -i 's/[[:blank:]]*$//' cgr_stringtie_assembly_cds.gtf

chx_bam=../plastid_analysis/merged_files/mapped/riboseq_chx.ncbi.Aligned.sortedByCoord.out.bam
nodrug_bam=../plastid_analysis/merged_files/mapped/riboseq_nodrug.ncbi.Aligned.sortedByCoord.out.bam
harr_bam=../plastid_analysis/merged_files/mapped/riboseq_harr.ncbi.Aligned.sortedByCoord.out.bam
fasta=../reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna

 
for source in reference stringtie
do
  if [[ $source -ge "reference" ]]
    then
    gtf=../reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf
    echo $gtf
   else
     gtf=cgr_stringtie_assembly_cds.gtf
     echo $gtf
    fi

  /mnt/HDD2/colin/bin/kentUtils/bin/linux.x86_64/gtfToGenePred \
  -ignoreGroupsWithoutExons -allErrors $gtf stdout | \
  /mnt/HDD2/colin/bin/kentUtils/bin/linux.x86_64/genePredToBed stdin \
  cgr.orfrater.annotation1.bed

  grep -v NW_023277000.1  cgr.orfrater.annotation1.bed > cgr.orfrater.annotation.$source.bed

prune_transcripts.py \
--inbed  cgr.orfrater.annotation.$source.bed \
--summarytable tid_removal_summary.txt \
-p 32 \
--minlen 28 \
--maxlen 31 \
-v $fasta \
$chx_bam \
$nodrug_bam \
--pseudogenes pseudogene_gene_ids \
--outbed transcripts.bed \
--force > 1prune.$source.log

make_tfams.py -v --force --inbed transcripts.bed \
--tfamstem tfams >   2make_tfams.$source.log

# potential start codon list from Kearse et al. Genes Dev 2017
rm orf.h5 # need to remove this as problem with indexing if exists
find_orfs_and_types.py $fasta --codons NTG -p 32 -v --force \
--tfamstem tfams --inbed transcripts.bed \
--orfstore orf.h5 > 3find_ORF.$source.log

psite_trimmed.py $chx_bam \
--minrdlen 28 \
--maxrdlen 31 \
--subdir chx \
--tallyfile tallies.txt \
--cdsbed cgr.orfrater.annotation.$source.bed \
-p 25 \
-v \
--force > psite_chx.$source.log

psite_trimmed.py $nodrug_bam \
--minrdlen 28 \
--maxrdlen 31 \
--subdir no_drug \
--tallyfile tallies.txt \
--cdsbed cgr.orfrater.annotation.$source.bed \
-p 32 \
-v \
--force > psite_nogrug.log

psite_trimmed.py $harr_bam \
--minrdlen 28 \
--maxrdlen 31 \
--subdir harr \
--tallyfile tallies.txt \
--cdsbed cgr.orfrater.annotation.$source.bed \
-p 32 \
-v \
--force > psite_harr.$source.log

regress_orfs.py $harr_bam \
--startonly \
--subdir harr \
--orfstore orf.h5 \
--inbed transcripts.bed \
--regressfile regression.h5 \
-p 32 \
-v \
--startcount 1 \
--force > regress_start.$source.log

regress_orfs.py $chx_bam  \
--subdir chx \
--orfstore orf.h5 \
--inbed transcripts.bed \
--restrictbystarts harr \
--startcount 1 \
-p 32 \
-v \
--force > chx_regress_stop.$source.log

regress_orfs.py $nodrug_bam \
--subdir no_drug \
--orfstore orf.h5 \
--inbed transcripts.bed \
--restrictbystarts harr \
--startcount 1 \
-p 32 \
-v \
--force > nodrug_regress_stop.$source.log

rate_regression_output.py \
harr \
chx \
no_drug \
--orfstore orf.h5 \
--ratingsfile orfratings.h5 \
-p 32 \
--CSV rate_regression.$source.csv \
-v \
--force > rate.$source.log

make_orf_bed.py --minlen 7 --force --outbed orfrater_predictions.$source.bed

done


/mnt/HDD2/colin/bin/kentUtils/bin/linux.x86_64/bedToGenePred orfrater_predictions.reference.bed stdout | \
/mnt/HDD2/colin/bin/kentUtils/bin/linux.x86_64/genePredToGtf file stdin orfrater_predictions.gtf

grep -f novel_proteoforms.txt orfrater_predictions.reference.bed > novel_proteoforms.bed 
# make a nucleotide fasta for the CDS gregion - from the orf-rater
cat novel_proteoforms.bed | awk -v type=CDS -f ../scripts/bed12toAnnotation.awk > novel_orfrater_predictions.bed12

bedtools getfasta -fi ../reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna -bed novel_orfrater_predictions.bed12 -fo orfrater_predictions.nucleotide.fasta -split -s -name

# convert to protein fasta
transeq -seq orfrater_predictions.nucleotide.fasta -frame 1 -outseq orfrater_predictions.protein.fasta

# join to uniprot for MS
cat orfrater_predictions.protein.fasta uniprot-proteome_UP000001075.fasta > extended_cgr_database_2.fasta


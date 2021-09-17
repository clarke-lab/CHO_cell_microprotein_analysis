#!/bin/bash
#### Description: Create a fasta file with new ORF predictions 
####              for Mass spec HCP analysis 
#### 
#### Written by: Clarke Lab. NIBRT - colin.clarke@nibrt.ie 

mkdir proteomics

### 1 Dowload the uniprot database for CGR
uniprot_url=ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000001075
wget $uniprot_url/UP000001075_10029.fasta.gz -P proteomics && gunzip proteomics/UP000001075_10029.fasta.gz


kent_path=/mnt/HDD2/colin/bin/kentUtils/bin/linux.x86_64

$kent_path/bedToGenePred \
orfrater_analysis/orfrater_predictions.reference.bed stdout | \

$kent_path/genePredToGtf file stdin \
proteomics/unfiltered_orfrater_predictions.gtf

grep -f proteomics/orfs_for_proteomics.txt \
 proteomics/unfiltered_orfrater_predictions.gtf > \
 proteomics/orfrater_predictions.gtf

# convert the ORFRATER GTF output to fasta file
agat_sp_extract_sequences.pl \
--gff proteomics/orfrater_predictions.gtf \
--fasta reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna \
-t cds  \
-o proteomics/orfrater_cds_protein.fa \
-p

# join the uniprot and orfrater protein sequences 
cat proteomics/UP000001075_10029.fasta proteomics/orfrater_cds_protein.fa > \
proteomics/combined_reference_proteome.fa
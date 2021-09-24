#!/bin/bash
#### Description: Create a fasta file with new ORF predictions 
####              for Mass spec HCP analysis 
#### 
#### Written by: Clarke Lab. NIBRT - colin.clarke@nibrt.ie 

### 1 Dowload the uniprot database for CGR
uniprot_url=ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000001075
wget $uniprot_url/UP000001075_10029.fasta.gz -P proteomics_db && gunzip proteomics_db/UP000001075_10029.fasta.gz

grep -f orf_lists/orfs_for_proteomics.txt \
orfrater_analysis/orfrater_predictions.reference.gtf > \
proteomics_db/novel_orfs.gtf

# convert the ORFRATER GTF output to fasta file
agat_sp_extract_sequences.pl \
--gff proteomics_db/novel_orfs.gtf \
--fasta reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna \
-t cds  \
-o proteomics_db/novel_orfs.fa \
-p

# join the uniprot and orfrater protein sequences 
cat proteomics_db/UP000001075_10029.fasta proteomics_db/novel_orfs.fa > \
proteomics_db/proteomics_db.fa
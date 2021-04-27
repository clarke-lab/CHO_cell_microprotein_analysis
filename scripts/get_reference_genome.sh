#!/bin/bash

# retrieve the PICR reference genome and annotation from ensembl
mkdir -p reference_genome

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/668/045/GCF_003668045.3_CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna.gz \
-P reference_genome

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/668/045/GCF_003668045.3_CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf.gz \
-P reference_genome

gunzip reference_genome/*.gz

# # make a gtf file for annotated protein coding genes for initial xtail count
# grep protein_coding reference_genome/cgr_ensembl_v103.gtf > \
# reference_genome/cgr_ensembl_v103_protein_coding.gtf



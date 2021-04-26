#!/bin/bash

# retrieve the PICR reference genome and annotation from ensembl
mkdir -p reference_genome

wget http://ftp.ensembl.org/pub/release-103/fasta/cricetulus_griseus_picr/dna/Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.fa.gz \
-P reference_genome

mv  reference_genome/Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.fa.gz \
reference_genome/cgr_ensembl_v103.fa.gz

wget http://ftp.ensembl.org/pub/release-103/gtf/cricetulus_griseus_picr/Cricetulus_griseus_picr.CriGri-PICR.103.gtf.gz \
-P reference_genome

mv reference_genome/Cricetulus_griseus_picr.CriGri-PICR.103.gtf.gz \
reference_genome/cgr_ensembl_v103.gtf.gz

gunzip reference_genome/*.gz


# make a gtf file for annotated protein coding genes for initial xtail count
grep protein_coding reference_genome/cgr_ensembl_v103.gtf > \
reference_genome/cgr_ensembl_v103_protein_coding.gtf

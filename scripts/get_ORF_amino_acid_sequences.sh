
#!/bin/bash
#### Description: output protein sequences for selected ORFs   
####                1. Annotated ORF-RATER ORFs
####                2. Upstream ORFs
####                3. Non-coding ORFs
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

# 1. Annotated ORF-RATER ORFs
grep -f orf_lists/long_orfs_for_AA_freq.txt \
orfrater_analysis/orfrater_predictions.reference.gtf > amino_acid_analysis/long_orfs_for_AA_freq.gtf

agat_sp_extract_sequences.pl \
--gff amino_acid_analysis/long_orfs_for_AA_freq.gtf \
--fasta reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna \
-t cds  \
-o amino_acid_analysis/long_orfs_for_AA_freq.fasta \
-p 

# 2. Upstream ORFs
grep -f orf_lists/upstream_short_orfs_for_AA_freq.txt \
orfrater_analysis/orfrater_predictions.reference.gtf > amino_acid_analysis/upstream_short_orfs_for_AA_freq.gtf

agat_sp_extract_sequences.pl \
--gff amino_acid_analysis/upstream_short_orfs_for_AA_freq.gtf \
--fasta reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna \
-t cds  \
-o amino_acid_analysis/upstream_short_orfs_for_AA_freq.fasta \
-p 

# 3. Non-coding ORFs
grep -f orf_lists/non_coding_orfs.txt \
orfrater_analysis/orfrater_predictions.reference.gtf > amino_acid_analysis/non_coding_orfs.gtf

agat_sp_extract_sequences.pl \
--gff amino_acid_analysis/non_coding_orfs.gtf \
--fasta reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna \
-t cds  \
-o amino_acid_analysis/non_coding_orfs.fasta \
-p 
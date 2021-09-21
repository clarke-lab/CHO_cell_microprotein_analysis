

grep -f proteomics/long_orfs_for_AA_freq.txt \
orfrater_analysis/orfrater_predictions.reference.gtf > proteomics/long_orfs_for_AA_freq.gtf

agat_sp_extract_sequences.pl \
--gff proteomics/long_orfs_for_AA_freq.gtf \
--fasta reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna \
-t cds  \
-o proteomics/long_orfs_for_AA_freq.fasta \
-p 

grep -f proteomics/short_orfs_for_AA_freq.txt \
orfrater_analysis/orfrater_predictions.reference.gtf > proteomics/short_orfs_for_AA_freq.gtf

agat_sp_extract_sequences.pl \
--gff proteomics/short_orfs_for_AA_freq.gtf \
--fasta reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna \
-t cds  \
-o proteomics/short_orfs_for_AA_freq.fasta \
-p 

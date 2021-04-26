#!/bin/bash

export PATH=$PATH:/mnt/HDD2/colin/bin/TransDecoder/
transdecoder_path=/mnt/HDD2/colin/bin/TransDecoder/

$transdecoder_path/util/gtf_genome_to_cdna_fasta.pl \
transcriptome_assembly/cgr_stringtie_assembly_raw.gtf \
reference_genome/cgr_ensembl_v103.fa > \
transdecoder/transcripts.fasta

$transdecoder_path/util/gtf_to_alignment_gff3.pl \
transcriptome_assembly/cgr_stringtie_assembly_raw.gtf > \
transcripts.gff3

TransDecoder.LongOrfs -t transcripts.fasta -S

mkdir transdecoder/swiss_prot

wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz \
-P transdecoder/swiss_prot
gunzip transdecoder/swiss_prot/uniprot_sprot.fasta.gz

makeblastdb -in  transdecoder/swiss_prot/uniprot_sprot.fasta -dbtype prot 

mkdir transdecoder/pfam
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz \
-P transdecoder/pfam
gunzip transdecoder/pfam/Pfam-A.hmm.gz

hmmpress transdecoder/pfam/Pfam-A.hmm

blastp \
-query transdecoder/longest_orfs.pep \
-db transdecoder/swiss_prot/uniprot_sprot.fasta  \
-max_target_seqs 1 \
-outfmt 6 \
-evalue 1e-5 \
-num_threads 16 \
> transdecoder/swiss_prot/blastp.outfmt6

hmmscan \
--cpu 16 \
--domtblout transdecoder/pfam/pfam.domtblout \
transdecoder/pfam/Pfam-A.hmm \
transdecoder/longest_orfs.pep


TransDecoder.Predict -t transcripts.fasta \
--retain_pfam_hits pfam/pfam.domtblout \
--retain_blastp_hits swiss_prot/blastp.outfmt6 \
--single_best_only

$transdecoder_path/util/cdna_alignment_orf_to_genome_orf.pl \
     transcripts.fasta.transdecoder.gff3 \
     transcripts.gff3 \
     transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3

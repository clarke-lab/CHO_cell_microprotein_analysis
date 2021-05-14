gmap_build -D ./ -d NCBI ../reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna

gmap -D ./ \
     -d NCBI \
     -n 1 \
     --no-chimeras \
     -f gff3_gene \
     -t 32 \
      lncrna.transcripts.fasta > gmap.gff3
#!/bin/bash

export PATH=$PATH:/mnt/HDD2/colin/bin/Price_1.0.3b
export PATH=$PATH:/mnt/HDD2/colin/bin/STAR-2.7.2d/bin/Linux_x86_64
export PATH=$PATH:/mnt/HDD2/colin/bin/kallisto/
export PATH=$PATH:/mnt/HDD2/colin/bin/kraken2/kraken2
export PATH=$PATH:/mnt/HDD2/colin/bin/reaper/


# get a list of ensembl transcripts from the stringtie cgr_stringtie_assembly
grep ENSCGRT cgr_stringtie_assembly_raw.gtf | awk -F'\t' '{if($3 == "transcript")print $0}' |
awk -F'\t' '{print $9}' | \
awk -F';' '{print $1,$2}' | \
sed s/transcript_id//g | \
sed s/gene_id[[:space:]]//g | \
sed s/\"//g > ens_transcript_stringtie_gene_map

# replace the ENSCGRG ID for each CDS annotatoion with the MSTRG gene ID
# takes some time
while read -ra a ;
do
grep ${a[1]} cgr_ensembl_v103.cds.annotations.gtf > temp.file
mstrg_id=${a[0]}
sed  "s/ENSCGRG.[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]/${mstrg_id}/;" temp.file >> annotated_cds_with_mstrg.gtf
done < ens_transcript_stringtie_gene_map

cat cgr_stringtie_assembly_raw.gtf \
../ribosome_footprint_profiling/transcriptome_assembly/annotateDd_cds_with_mstrg.gtf > cgr_price.gtf


# agat_convert_sp_gff2gtf.pl --in cgr_price.gtf --gtf_version 2.5 -o agat_test2.gtf

gedi Nashorn -e 'new GtfFileReader("price_reference/.gtf","CDS").readIntoMemoryThrowOnNonUnique().ei().print()'


  gedi \
  -e IndexGenome \
  -s price_reference/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna \
  -a price_reference/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf \
  -n cgr_price \
  -p

    gedi \
  -e IndexGenome \
  -s price_rRNA/cgr_rRNA_rnacentral_v17.fasta \
  -n cgr.rrna \
  -p


  gedi -e Price -reads merged.bamlist \
  -prefix orf_predictions \
  -genomic cgr_price \
  -filter 28:31 \
  -novelTranscripts \
  -nthreads 32 \
  -fdr 0.05 \
  -D

    gedi -e Price -reads individual.bamlist \
    -prefix chx_onlu/orf_predictions \
    -genomic cgr_price \
    -nthreads 20 \
    -fdr 0.05 \
    -D




gedi Nashorn -e 'load("orf_predictions.orfs.cit").ei().map(function(o) new BedEntry(o.data.getStartStop(o,true).toMutable().setData(new NameAnnotation(o.data.getTranscript()+"_"+o.data.getType()+"_"+o.data.getOrfid())))).print()' > unfiltered.bed.file


gedi -e Pipeline -r serial --tmp=/mnt/HDD2/colin/ribosome_footprint_profiling/price_analysis/price_tmp -j chx_test.json shortread_mapping.sh resolve_ambi.sh report.sh price.sh

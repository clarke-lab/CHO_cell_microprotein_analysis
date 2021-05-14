import numpy
from plastid import *
transcript_reader = GTF2_TranscriptAssembler("reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf")
fout_start = open("diff_translation_analysis/reference_annotation/reference_start_codon_masks.bed","w")
fout_stop = open("diff_translation_analysis/reference_annotation/reference_stop_codon_masks.bed","w")

for feature in transcript_reader:
  cds = feature.get_cds()
  if cds.get_length() <= 100:
    
    start_codon_mask = list(cds.get_subchain(0,1))
    stop_codon_mask  = list(cds.get_subchain(cds.length-1,cds.length))
  else:
      start_codon_mask = list(cds.get_subchain(0,15))
      stop_codon_mask  = list(cds.get_subchain(cds.length-5,cds.length))

  for mask in start_codon_mask:
     fout_start.write(SegmentChain(mask).as_bed())
  for mask in stop_codon_mask:
     fout_stop.write(SegmentChain(mask).as_bed())

fout_start.close()
fout_stop.close()

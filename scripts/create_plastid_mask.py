import numpy
from plastid import * 

transcript_reader = GTF2_TranscriptAssembler(open("plastid_reference/plastid_annotations.gtf"))
fout_start = open("plastid_reference/start_codon_masks.bed","w")
fout_stop = open("plastid_reference/stop_codon_masks.bed","w")

for feature in transcript_reader:
  cds = feature.get_cds()
  if cds.get_length() <= 300: # for small proteins exclude first and last codon
    # print(feature.get_name(),feature.attr["type"],str(feature)) # do something
    # print(cds.get_length()) 
    start_codon_mask = list(cds.get_subchain(0,3))
    stop_codon_mask  = list(cds.get_subchain(cds.length-3,cds.length))
  else: # for other proteins the first 15 and last 5 codons are not counted 
      start_codon_mask = list(cds.get_subchain(0,45))
      stop_codon_mask  = list(cds.get_subchain(cds.length-15,cds.length))
  for mask in start_codon_mask:
     fout_start.write(SegmentChain(mask).as_bed())
  for mask in stop_codon_mask:
     fout_stop.write(SegmentChain(mask).as_bed())

fout_start.close()
fout_stop.close()

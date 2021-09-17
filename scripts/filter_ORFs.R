#!/usr/bin/env Rscript --vanilla
#### Description: Filters the raw ORF-RATER output
####              
#### 
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie


package_list <- c("tidyverse","writexl")

lapply(package_list, require, character.only = TRUE)

root_dir <-c("/mnt/HDD2/colin/ribosome_footprint_profiling/")

results_dir <- paste0(root_dir,"results/section3.2/")
if (!dir.exists(results_dir)) {
  dir.create(results_dir,recursive = TRUE)
}

orftable <- read.csv(file = paste0(root_dir, "/orfrater_analysis/rate_regression.ncbi.csv"), header = T)

# 1. filter by score and length
# 2. Remove truncation, internal, LOOF types
# 3. Modify ORF names
orftable_filtered <- orftable %>%
  filter(orfrating >= 0.5 & AAlen >= 5) %>%
  filter(orftype != "truncation" & 
           orftype != "internal" & 
           orftype != "LOOF") %>% 
  mutate(orftype = case_when(
    orftype ==  "Giso" ~ "Isoform",
    orftype ==  "NCiso" ~ "Isoform",
    orftype ==  "new_iso" ~ "Isoform",
    orftype ==  "Niso" ~ "Isoform",
    orftype ==  "Siso" ~ "Isoform",
    orftype ==  "Niso" ~ "Isoform",
    orftype ==  "Xiso" ~ "Isoform",
    orftype ==  "Ciso" ~ "Isoform",
    orftype == "annotated" ~ "Annotated",
    orftype == "new" ~ "New",
    orftype == "upstream" ~ "Upstream",
    orftype == "downstream" ~ "Downstream",
    orftype == "start_overlap" ~ "Start overlap",
    orftype == "stop_overlap" ~ "Stop overlap",
    orftype == "extension" ~ "Extension"
  ))

table(orftable_filtered$orftype)

# Filter potential false positives
non_truncates <- orftable_filtered %>% 
  filter(orftype != "Upstream") %>% 
  filter(orftype != "Start overlap") %>%
  filter(orftype != "New") 

# filter uORFs
single_proteoform_transcripts <-orftable_filtered %>% 
  filter(orftype == "Upstream" | orftype == "Start overlap" | orftype == "New") %>% 
  group_by(tid) %>% 
  mutate(count=n()) %>% 
  filter(count == 1) 

overlapping <- function(x){
  if (length(unique(x$tstop)) < length(x$tstop)) {
    n_occur <- data.frame(table(x$tstop))
    x_d <- x[duplicated(x$tstop),] 
    x_d <- x_d[which.max(x_d$tstop-x_d$tcoord),] 
    if (length(n_occur$Freq == 1) > 0) { 
      nd_tstop <- n_occur[n_occur$Freq <2,1] 
      x_nd <- x[x$tstop %in% nd_tstop,] 
      x_final <- rbind(x_nd, x_d) } 
    else { 
      x_final <- x_d
      }
    } else  {
      x_final <- x
      }
  x_final
  }

uORF_filtered <- orftable_filtered %>%
  filter(orftype == "Upstream") %>%
  group_by(tid) %>%
  mutate(count=n()) %>%
  filter(count > 1) %>%
  group_modify(~overlapping(.x)) %>%
  dplyr::select(-count) %>%
  dplyr::select(colnames(orftable_filtered))

ouORF_filtered <- orftable_filtered %>%
  filter(orftype == "Start overlap") %>%
  group_by(tid) %>%
  mutate(count=n()) %>%
  filter(count > 1) %>%
  group_modify(~overlapping(.x)) %>%
  dplyr::select(-count) %>%
  dplyr::select(colnames(orftable_filtered))

new_filtered <- orftable_filtered %>%
  filter(orftype == "New") %>%
  group_by(tid) %>%
  mutate(count=n()) %>%
  filter(count > 1) %>%
  group_modify(~overlapping(.x)) %>%
  dplyr::select(-count) %>%
  dplyr::select(colnames(orftable_filtered))

# join the filter ORFs 
orftable_filtered <- bind_rows(non_truncates, 
                               single_proteoform_transcripts, 
                               uORF_filtered, 
                               ouORF_filtered, 
                               new_filtered)

# write the ids to file for transcript level quantation using plastid
write(orftable_filtered$orfname, file=paste0(root_dir,"quantitation/transcript/filtered_orfs.txt"))

# where there are multiple new-transcripts take the longest
# this is for plastid counting requirment for 1 ORF per transcript
# select new ORFs on non-coding RNAs onle
selected_new_for_de <- orftable_filtered %>% 
  filter(orftype == "New") %>% 
  filter(str_detect(orfname, "NR|XR")) %>%
  #filter(AAlen >= 5 & AAlen <= 100) %>%
  group_by(tid) %>%
  top_n(AAlen, n=1) 

# write the ids to allow filtering of the GTF file
write(selected_new_for_de$orfname, file=paste0(root_dir,"quantitation/gene/non_coding_encoded_orfs.txt"))

pcg_orfs <- orftable_filtered %>% 
  filter(orftype == "Annotated" | orftype == "Extension" | orftype == "Isoform" )

write(pcg_orfs$orfname, file=paste0(root_dir,"quantitation/transcript/pcg_orfs.txt"))

# selecting sORFs for the proteomics analysis
selected_for_proteomics <- orftable_filtered %>%
  filter(orftype != "Annotated" & orftype != "Isoform" & orftype != "Extension")

write(selected_for_proteomics$orfname, file=paste0(root_dir,"proteomics/orfs_for_proteomics.txt"))

  selected_annotated <- orftable_filtered %>% 
    filter(orftype == "Annotated") 

write(selected_annotated$orfname, file=paste0(root_dir,"quantitation/transcript/selected_annotated.txt"))

selected_extension <- orftable_filtered %>% 
  filter(orftype == "Extension") %>%
  group_by(tid) %>%
  top_n(AAlen, n=1) 
  
write(selected_extension$orfname, file=paste0(root_dir,"quantitation/transcript/selected_extension.txt"))
  
selected_isoform <- orftable_filtered %>% filter(orftype == "Isoform") %>%
  group_by(tid) %>%
  top_n(AAlen, n=1) 
  
write(selected_isoform$orfname, file=paste0("quantitation/transcript/selected_isoform.txt"))
  
selected_upstream <- orftable_filtered %>% filter(orftype == "Upstream") %>%
  group_by(tid) %>%
  top_n(AAlen, n=1) 
   
write(selected_upstream$orfname, file=paste0("quantitation/transcript/selected_upstream.txt"))
    
selected_overlap <- orftable_filtered %>% filter(orftype == "Start overlap") %>%
  group_by(tid) %>%
  top_n(AAlen, n=1) 

write(selected_overlap$orfname, file=paste0("quantitation/transcript/selected_overlap.txt"))

write(selected_new_for_de$orfname, file=paste0("quantitation/transcript/selected_new.txt"))
   
CriGri_PICRH_1_0_annotation <- read_delim(paste0("reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_feature_table.txt"), 
          "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
          mutate(tfam= `product_accession`) 

table_s2 <- left_join(orftable_filtered, CriGri_PICRH_1_0_annotation, by="tfam") 

table_s2 <-table_s2 %>%
  dplyr::select(
    `ORF-RATER name` = orfname,
    `ORF type` = orftype,
    `ORF-RATER score` = orfrating,
    `Associated Gene symbol` = symbol, 
    `Associated Gene Name` = name, 
    `Transcript family` = tfam,
    `Transcript ID` = tid,
    `Start codon` = codon,
    `Length (AAs)` = AAlen,
    `Start Annotated?` = annot_start,
    `Stop Annotated` = annot_stop,
    `Transcript start position` = tcoord,
    `Transcript stop position` =  tstop,
    `Chromosome` = chrom,
    `Strand` = strand.x,
    `Genomic start position` = tcoord,
    `Genomic stop position` =  tstop,
  ) %>%
  arrange(-`ORF-RATER score`)

write_xlsx(list(ORFs = table_s2),
  path = paste(results_dir,"Table S2.xlsx", sep=""),
  format_headers = TRUE)

#save(table_s2, file = paste(results_dir,"results_3_2.RData", sep=""))
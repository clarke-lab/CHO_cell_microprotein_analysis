#!/usr/bin/env Rscript --vanilla
#### Description: Filters the raw ORF-RATER output
####              
#### 
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

package_list <- c("tidyverse","writexl")

invisible(lapply(package_list, require, character.only = TRUE, quietly=TRUE))

root_dir <-c("/mnt/HDD2/colin/ribosome_footprint_profiling/")

results_dir <- paste0(root_dir,"results/section2.2/")
if (!dir.exists(results_dir)) {
  dir.create(results_dir,recursive = TRUE)
}

orftable <- read.csv(file = paste0(root_dir, "/orfrater_analysis/rate_regression.ncbi.csv"), header = T)

orftable_filtered <- orftable %>%
  filter(orfrating >= 0.5 & AAlen >= 5)

print("removal of aa<5 & score<5")
table(orftable_filtered$orftype)
print("####")
sum(table(orftable_filtered$orftype))
print("####")

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
print("removal of truc,internal, LOOF")
table(orftable_filtered$orftype)
print("####")
sum(table(orftable_filtered$orftype))
print("####")
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

print("Chinese hamster ORFs Identified")
table(orftable_filtered$orftype)
print("------")
# add annotation information to identified ORFs
suppressWarnings(
  suppressMessages(
    CriGri_PICRH_1_0_annotation <- read_delim(paste0("reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_feature_table.txt"), 
          "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
          mutate(tfam= `product_accession`)
          )
          )

table_s3 <- left_join(orftable_filtered, CriGri_PICRH_1_0_annotation, by="tfam") 

table_s3 <-table_s3 %>%
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
    `Genomic stop position` =  tstop) %>%
  arrange(-`ORF-RATER score`)

write_xlsx(list(ORFs = table_s3),
  path = paste(results_dir,"Table S3.xlsx", sep=""),
  format_headers = TRUE)

save(table_s3, file = paste(results_dir,"results_2_2.RData", sep=""))


## output ORF lists for different analyses
print("Outputting ORF lists for downstream analysis")

# 1. amino acid frequencies
print("--Amino acid analysis--")
# annotated PCGs
long_orfs <- table_s3 %>%
  filter(`Length (AAs)` > 100) %>%
  filter(`ORF type` == "Annotated") 
write(long_orfs$`ORF-RATER name`, file=paste0("orf_lists/long_orfs_for_AA_freq.txt"))
paste0(length(long_orfs$`ORF-RATER name`), 
  " ORFs written to ", 
  paste0(root_dir,"orf_lists/long_orfs_for_AA_freq.txt"))


# non-coding RNA ORFs 
new_lncrna_short_orfs <-table_s3 %>% 
  filter(`Length (AAs)` <= 100) %>%
  filter(`ORF type` == "New") %>%
  filter(str_detect(`Transcript ID`, "XR|NR"))
write(new_lncrna_short_orfs$`ORF-RATER name`, file="orf_lists/lncrna_short_orfs_for_AA_freq.txt")
paste0(length(new_lncrna_short_orfs$`ORF-RATER name`), 
  " ORFs written to ", 
  paste0(root_dir,"orf_lists/lncrna_short_orfs_for_AA_freq.txt"))


# upstream short ORFs
upstream_short_orfs <-table_s3 %>% 
  filter(`Length (AAs)` <= 100) %>%
  filter(`ORF type` == "Upstream" | `ORF type` == "Start overlap") 

write(upstream_short_orfs$`ORF-RATER name`, file="orf_lists/upstream_short_orfs_for_AA_freq.txt")
paste0(length(upstream_short_orfs$`ORF-RATER name`), 
  " ORFs written to ", 
  paste0(root_dir,"orf_lists/upstream_short_orfs_for_AA_freq.txt"))

# 2. gene level differential expression analysis
print("--Differential expression--")
# where there are multiple new-transcripts take the longest
# this is for plastid counting requirment for 1 ORF per transcript
# select new ORFs on non-coding RNAs onle
selected_new_for_de <- table_s3 %>% 
  filter(`ORF type` == "New") %>% 
  filter(str_detect(`ORF-RATER name`, "NR|XR")) %>%
  #filter(AAlen >= 5 & AAlen <= 100) %>%
  group_by(`Transcript ID`) %>%
  top_n(`Length (AAs)`, n=1) 

# write the ids to allow filtering of the GTF file
write(selected_new_for_de$`ORF-RATER name`, file=paste0(root_dir,"orf_lists/non_coding_orfs.txt"))
paste0(length(selected_new_for_de$`ORF-RATER name`), 
  " ORFs written to ", 
  paste0(root_dir,"orf_lists/non_coding_orfs_for_de.txt"))


# 3. Proteomics
print("--Proteomics--")
selected_for_proteomics <- table_s3 %>%
  filter(`ORF type` != "Annotated" & `ORF type` != "Isoform" & `ORF type` != "Extension")
write(selected_for_proteomics$`ORF-RATER name`, file="orf_lists/orfs_for_proteomics.txt")
paste0(length(selected_for_proteomics$`ORF-RATER name`), 
  " ORFs written to ", paste0(root_dir,
  "orf_lists/orfs_for_proteomics"))
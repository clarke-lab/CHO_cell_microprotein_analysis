#!/usr/bin/env Rscript --vanilla
#### Description: Filters the raw ORF-RATER output
####              
#### 
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

package_list <- c(
  "tidyverse", "DESeq2", "writexl", "ggpubr",
  "ggthemes", "viridis", "patchwork", "WebGestaltR", "fuzzyjoin",
  "GenomicFeatures", "wiggleplotr", "readr", "heatmaply", "ggvenn", "data.table", "ggrepel"
)

lapply(package_list, require, character.only = TRUE)

root_dir <- c("/mnt/HDD2/colin/ribosome_footprint_profiling/")
results_dir <- paste0(root_dir,"results/section3.5/")
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

te_combined_data <- read.delim(paste0(root_dir, "quantitation/gene/count_output/combined_counts.txt"),
  sep = "\t",
  header = TRUE,
  row.names = "region"
)

te_sample_table <- read.csv(paste0(root_dir, "data/de_translation_design.txt"),
  header = TRUE,
  row.names = "sample_name"
)

# import the annotation
reference_annotation_plastid_merged <- read_delim(paste0(root_dir, "diff_translation_analysis/reference_annotation_merged.txt"),
  "\t",
  escape_double = FALSE, col_names = FALSE,
  trim_ws = TRUE, skip = 27
) %>%
  mutate(X3 = str_remove(X1, "gene_")) %>%
  mutate(`product_accession` = str_remove(X3, "_1"))

# import the ncbi annotation
CriGri_PICRH_1_0_annotation <- read_delim(paste0(root_dir, "reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_feature_table.txt"),
  "\t",
  escape_double = FALSE, trim_ws = TRUE
)

# load the mouse annotation - used to match gene names to complete LOC to symbol conversion
mouse_feature_table <- read_delim(paste0(root_dir, "reference_genome/GCF_000001635.27_GRCm39_feature_table.txt"),
  "\t",
  escape_double = FALSE, trim_ws = TRUE
)

# set the differential expression counts
reference_annotation_plastid_merged <- read_delim(paste0(root_dir, "diff_translation_analysis/reference_annotation_merged.txt"),
  "\t",
  escape_double = FALSE, col_names = FALSE,
  trim_ws = TRUE, skip = 27
) %>%
  mutate(X3 = str_remove(X1, "gene_")) %>%
  mutate(`product_accession` = str_remove(X3, "_1"))

# import the ncbi annotation
CriGri_PICRH_1_0_annotation <- read_delim(paste0(root_dir, "reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_feature_table.txt"),
  "\t",
  escape_double = FALSE, trim_ws = TRUE
)

# load the mouse annotation - used to match gene names to complete LOC to symbol conversion
mouse_feature_table <- read_delim(paste0(root_dir, "reference_genome/GCF_000001635.27_GRCm39_feature_table.txt"),
  "\t",
  escape_double = FALSE, trim_ws = TRUE
)

fold_change_threshold <- 1.5
pval_threshold <- 0.05
base_mean_threshold <- 0
average_counts <- 20

# RNA-seq only differential expression
rnaseq_data <- te_combined_data %>%
  dplyr::select(contains("rnaseq"))

# retain genes with a least 10 counts
detected.rna <- rownames(rnaseq_data)[rowMeans(rnaseq_data) >= average_counts]
rnaseq_data <- rnaseq_data[detected.rna, ]

rnaseq_design <- te_sample_table %>%
  filter(assay == "rnaseq")

# 2) RNASeq DESeq
dds_rnaseq <- DESeqDataSetFromMatrix(
  countData = rnaseq_data,
  colData = rnaseq_design,
  design = ~condition
)

dds_rnaseq$condition <- relevel(dds_rnaseq$condition, ref = "nts")
dds_rnaseq <- estimateSizeFactors(dds_rnaseq)
dds_rnaseq <- estimateDispersions(dds_rnaseq)

# track plotting scaling
bam_coverage_scaling_rnaseq <- 1 / sizeFactors(dds_rnaseq)

dds_rnaseq <- DESeq(dds_rnaseq)

res_rnaseq <- results(dds_rnaseq)

sig_res_rnaseq <- res_rnaseq %>%
  as_tibble(rownames = "geneid") %>%
  dplyr::filter(abs(log2FoldChange) >= log2(fold_change_threshold) &
    padj < pval_threshold & baseMean >= base_mean_threshold)

# annotate protein coding genes
pcg_sig_res_rnaseq <- sig_res_rnaseq %>%
  filter(!str_detect(geneid, "NR|XR")) %>%
  mutate(symbol = gsub("_.*", "", geneid)) %>%
  left_join(CriGri_PICRH_1_0_annotation, by = "symbol") %>%
  filter(`# feature` == "mRNA") %>%
  mutate(name = gsub(",.*", "", name)) %>%
  distinct(geneid, .keep_all = T) %>%
  dplyr::select("symbol", "name", "GeneID", c(colnames(sig_res_rnaseq)))

# gene names causes a problem with the regex in fuzzy match
pcg_sig_res_rnaseq <- pcg_sig_res_rnaseq %>%
  mutate(symbol = replace(symbol, name == "amine oxidase [flavin-containing] A", "Maoa"))

pcg_without_loc_ids <- pcg_sig_res_rnaseq %>%
  filter(!str_detect(symbol, "LOC"))

pcg_with_loc_ids <- pcg_sig_res_rnaseq %>%
  filter(str_detect(symbol, "LOC"))

# match the names against mouse
loc_id_gene_symbol <- mouse_feature_table %>%
  mutate(ncbi_symbol = symbol) %>%
  fuzzy_inner_join(pcg_with_loc_ids, by = "name", match_fun = str_detect) %>%
  distinct(symbol.y, .keep_all = T)

pcg_with_loc_ids <- pcg_with_loc_ids %>%
  mutate("symbol.y" = symbol) %>%
  left_join(loc_id_gene_symbol, by = "symbol.y") %>%
  mutate(ncbi_symbol = coalesce(ncbi_symbol, symbol.y)) %>%
  mutate(
    "geneid" = geneid.x,
    "baseMean" = baseMean.x,
    "log2FoldChange" = log2FoldChange.x,
    "lfcSE" = lfcSE.x, "stat" = stat.x,
    "pvalue" = pvalue.x, "padj" = padj.x,
    "symbol" = ncbi_symbol
  ) %>%
  dplyr::select(symbol, name, GeneID, c(colnames(sig_res_rnaseq)))

nc_sig_res_rnaseq <- sig_res_rnaseq %>%
  filter(str_detect(geneid, "NR|XR")) %>%
  mutate(
    symbol = "New",
    name = "New",
    GeneID = 0
  ) %>%
  dplyr::select(
    "symbol",
    "name",
    "GeneID",
    c(colnames(sig_res_rnaseq))
  )

sig_res_rnaseq <- bind_rows(pcg_without_loc_ids, pcg_with_loc_ids, nc_sig_res_rnaseq) %>%
  arrange(-log2FoldChange) %>%
  dplyr::select(c(
    "geneid",
    "GeneID",
    "symbol",
    "name",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj"
  )) %>%
  dplyr::rename(Plastid_ID = "geneid")


# RPF only differential expression
riboseq_data <- te_combined_data %>%
  dplyr::select(contains("riboseq"))

riboseq_design <- te_sample_table %>%
  filter(assay == "riboseq")

# at least 10 count required to retain a given gene
detected.ribo <- rownames(riboseq_data)[rowMeans(riboseq_data) >= average_counts]
riboseq_data <- riboseq_data[detected.ribo, ]

# DESeq2
dds_riboseq <- DESeqDataSetFromMatrix(
  countData = riboseq_data,
  colData = riboseq_design,
  design = ~condition
)

dds_riboseq$condition <- relevel(dds_riboseq$condition, ref = "nts")

dds_riboseq <- estimateSizeFactors(dds_riboseq)
dds_riboseq <- estimateDispersions(dds_riboseq)

bam_coverage_scaling_riboseq <- 1 / sizeFactors(dds_riboseq)

dds_riboseq <- DESeq(dds_riboseq)
res_riboseq <- results(dds_riboseq)

sig_res_riboseq <- res_riboseq %>%
  as_tibble(rownames = "geneid") %>%
  dplyr::filter(abs(log2FoldChange) > log2(fold_change_threshold) &
    padj < pval_threshold & baseMean >= base_mean_threshold) %>%
  arrange(-log2FoldChange)

# annotate protein coding genes
pcg_sig_res_riboseq <- sig_res_riboseq %>%
  filter(!str_detect(geneid, "NR|XR")) %>%
  mutate(symbol = gsub("_.*", "", geneid)) %>%
  left_join(CriGri_PICRH_1_0_annotation, by = "symbol") %>%
  filter(`# feature` == "mRNA") %>%
  mutate(name = gsub(",.*", "", name)) %>%
  distinct(geneid, .keep_all = T) %>%
  dplyr::select("symbol", "name", "GeneID", c(colnames(sig_res_riboseq)))

# gene names causes a problem with the regex in fuzzy match
pcg_sig_res_riboseq <- pcg_sig_res_riboseq %>%
  mutate(symbol = replace(symbol, name == "amine oxidase [flavin-containing] A", "Maoa"))

pcg_without_loc_ids <- pcg_sig_res_riboseq %>%
  filter(!str_detect(symbol, "LOC"))

pcg_with_loc_ids <- pcg_sig_res_riboseq %>%
  filter(str_detect(symbol, "LOC"))

# match the names against mouse
loc_id_gene_symbol <- mouse_feature_table %>%
  mutate(ncbi_symbol = symbol) %>%
  fuzzy_inner_join(pcg_with_loc_ids, by = "name", match_fun = str_detect) %>%
  distinct(symbol.y, .keep_all = T)

pcg_with_loc_ids <- pcg_with_loc_ids %>%
  mutate("symbol.y" = symbol) %>%
  left_join(loc_id_gene_symbol, by = "symbol.y") %>%
  mutate(ncbi_symbol = coalesce(ncbi_symbol, symbol.y)) %>%
  mutate(
    "geneid" = geneid.x,
    "baseMean" = baseMean.x,
    "log2FoldChange" = log2FoldChange.x,
    "lfcSE" = lfcSE.x, "stat" = stat.x,
    "pvalue" = pvalue.x,
    "padj" = padj.x,
    "symbol" = ncbi_symbol
  ) %>%
  dplyr::select(symbol, name, GeneID, c(colnames(sig_res_riboseq)))

nc_sig_res_riboseq <- sig_res_riboseq %>%
  filter(str_detect(geneid, "NR|XR")) %>%
  mutate(symbol = "New", name = "New", GeneID = 0) %>%
  dplyr::select("symbol", "name", "GeneID", c(colnames(sig_res_riboseq)))

sig_res_riboseq <- bind_rows(pcg_without_loc_ids, pcg_with_loc_ids, nc_sig_res_riboseq) %>%
  arrange(-log2FoldChange) %>%
  dplyr::select(c("geneid", "GeneID", "symbol", "name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")) %>%
  dplyr::rename("Plastid_ID" = geneid)

# Differential translation
detected.both <- intersect(detected.rna, detected.ribo)

length(detected.both)

# DESeq2
dds_te <- DESeqDataSetFromMatrix(
  countData = te_combined_data[detected.both, ],
  colData = te_sample_table,
  design = ~ assay + condition + assay:condition
)

dds_te$condition <- relevel(dds_te$condition, ref = "nts")

dds_te <- estimateSizeFactors(dds_te)
dds_te <- estimateDispersions(dds_te)
bam_coverage_scaling_te <- 1 / sizeFactors(dds_te)
# correct for RNASeq changes
dds_te <- nbinomLRT(dds_te,
  full = ~ assay + condition + assay:condition,
  reduced = ~ assay + condition
)

res_te <- DESeq(dds_te)
res_te <- results(dds_te, name = "condition_ts_vs_nts")

sig_res_te <- res_te %>%
  as_tibble(rownames = "geneid") %>%
  dplyr::filter(abs(log2FoldChange) > log2(fold_change_threshold) & padj < pval_threshold & baseMean >= base_mean_threshold) %>%
  arrange(-log2FoldChange)

pcg_sig_res_te <- sig_res_te %>%
  filter(!str_detect(geneid, "NR|XR")) %>%
  mutate(symbol = gsub("_.*", "", geneid)) %>%
  left_join(CriGri_PICRH_1_0_annotation, by = "symbol") %>%
  filter(`# feature` == "mRNA") %>%
  mutate(name = gsub(",.*", "", name)) %>%
  distinct(geneid, .keep_all = T) %>%
  dplyr::select("symbol", "name", "GeneID", c(colnames(sig_res_te)))

# gene names causes a problem with the regex in fuzzy match
pcg_sig_res_te <- pcg_sig_res_te %>%
  mutate(symbol = replace(symbol, name == "amine oxidase [flavin-containing] A", "Maoa"))

pcg_without_loc_ids <- pcg_sig_res_te %>%
  filter(!str_detect(symbol, "LOC"))

pcg_with_loc_ids <- pcg_sig_res_te %>%
  filter(str_detect(symbol, "LOC"))

# match the names against mouse
loc_id_gene_symbol <- mouse_feature_table %>%
  mutate(ncbi_symbol = symbol) %>%
  fuzzy_inner_join(pcg_with_loc_ids, by = "name", match_fun = str_detect) %>%
  distinct(symbol.y, .keep_all = T)

pcg_with_loc_ids <- pcg_with_loc_ids %>%
  mutate("symbol.y" = symbol) %>%
  left_join(loc_id_gene_symbol, by = "symbol.y") %>%
  mutate(ncbi_symbol = coalesce(ncbi_symbol, symbol.y)) %>%
  mutate("geneid" = geneid.x, "baseMean" = baseMean.x, "log2FoldChange" = log2FoldChange.x, "lfcSE" = lfcSE.x, "stat" = stat.x, "pvalue" = pvalue.x, "padj" = padj.x, "symbol" = ncbi_symbol) %>%
  dplyr::select(symbol, name, GeneID, c(colnames(sig_res_te)))

nc_sig_res_te <- sig_res_te %>%
  filter(str_detect(geneid, "NR|XR")) %>%
  mutate(symbol = "New", name = "New", GeneID = 0) %>%
  dplyr::select("symbol", "name", "GeneID", c(colnames(sig_res_te)))

sig_res_te <- bind_rows(pcg_without_loc_ids, pcg_with_loc_ids, nc_sig_res_te) %>%
  arrange(-log2FoldChange) %>%
  dplyr::select(c("geneid", "GeneID", "symbol", "name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")) %>%
  rename("Plastid_ID" = geneid)


save(res_rnaseq, res_riboseq, res_te,dds_te, file = paste(results_dir, "results.3.5.RData", sep = ""))
library("DESeq2")

te_combined_data <- read.delim("diff_translation_analysis/te_combined_reference_counts.txt",
                                sep="\t",
                                header=TRUE,
                                row.names="region")

te_sample_table <- read.csv("data/de_translation_design.txt",
                               header=TRUE,
                               row.names="sample_name")

dds <- DESeqDataSetFromMatrix(countData = te_combined_data,
                               colData = te_sample_table,
                               design = ~ assay + condition + assay:condition)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

dds <- nbinomLRT(dds,
                  full= ~ assay + condition + assay:condition,
                  reduced= ~ assay + condition )

res <- DESeq(dds)
res <- results(dds)
summary(res)

resOrdered <- res[order(res$padj),]

write.table(as.data.frame(resOrdered),
             sep="\t",quote=FALSE,
             file="diff_translation_analysis/te_change_tempshift.txt")


# sneak a peak
head(resOrdered)
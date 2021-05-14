import pandas as pd
sample_names = ["rnaseq_nts_r1",
                "rnaseq_nts_r2",
                "rnaseq_nts_r3",
                "rnaseq_nts_r4", 
                "rnaseq_ts_r1",
                "rnaseq_ts_r2",
                "rnaseq_ts_r3",
                "rnaseq_ts_r4",
                "riboseq_nts_r1",
                "riboseq_nts_r2",
                "riboseq_nts_r3",
                "riboseq_nts_r4",
                "riboseq_ts_r1",
                "riboseq_ts_r2",
                "riboseq_ts_r3",
                "riboseq_ts_r4"]

# load samples as DataFrames
samples = { K : pd.read_table("diff_translation_analysis/reference_counts/%s.txt" % K,sep="\t",header=0,comment="#",index_col=None) for K in sample_names }

# combine count columns to single DataFrame
combined_df = samples["rnaseq_nts_r1"][["region"]]
for k,v in samples.items():
    combined_df[k] = v["cds_reads"]


# save
combined_df.to_csv("diff_translation_analysis/te_combined_reference_counts.txt",sep="\t",header=True,index=False,
                   columns=["region","rnaseq_nts_r1","rnaseq_nts_r2","rnaseq_nts_r3","rnaseq_nts_r4",
                   "rnaseq_ts_r1","rnaseq_ts_r2","rnaseq_ts_r3","rnaseq_ts_r4", "riboseq_nts_r1","riboseq_nts_r2",
                   "riboseq_nts_r3", "riboseq_nts_r4","riboseq_ts_r1","riboseq_ts_r2","riboseq_ts_r3", "riboseq_ts_r4"])
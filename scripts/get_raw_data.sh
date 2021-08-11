#!/bin/bash
#### Description: Downloads the RIBO-seq and RNA-seq data 
####              from ENA  
#### 
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie 

# 1. Total RNASeq (Single-end)
mkdir -p data/rnaseq_se/raw_data

# 2. CHX riboseq
mkdir -p data/riboseq_chx/raw_data

# 3. Harr riboseq
mkdir -p data/riboseq_harr/raw_data

# 4. No drug riboseq
mkdir -p data/riboseq_nd/raw_data
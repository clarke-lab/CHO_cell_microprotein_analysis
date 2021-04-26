#!/bin/bash

# download the raw data from public repository

# 1. Total RNASeq (Paired-end)
mkdir -p data/rnaseq_pe/raw_data

# 2. Total RNASeq (Single-end)
mkdir -p data/rnaseq_se/raw_data

# 3. CHX riboseq
mkdir -p data/riboseq_chx/raw_data

# 4. Harr riboseq
mkdir -p data/riboseq_harr/raw_data

# 5. No drug riboseq
mkdir -p data/riboseq_nd/raw_data

# ribosome_footprint_profiling

# Dependancies
  * cutadapt 1.18
  * Trimmomatic-0.36
  * STAR-2.7.2d
  * stringtie-2.1.5
  * xtail 1.1.5

## Dowload the raw RiboSeq and RNASeq data
```bash
./scripts/get_raw_data.sh
```
## Download the Chinese hamster genome sequence and annotation
Get the reference genome and create a STAR index for mapping
```bash
./get_reference_genome.sh
```

## Preprocess and map sequencing data

```bash
# Trim adapters from the Cycoheximide, Harringtonine and No drug Ribosome footprint profiling data
./preprocess_rnaseq_pe.sh
./preprocess_riboseq.sh
```

## ORF-RATER 

docker pull 
./scripts/orfrater_analysis.sh
```

# 2. Gene level differences in translation efficiency
```bash
./xtail_analysis_gene.sh
```

# 3. Translation initiation site identification
```bash

```

# 1. Differential Translation Analysis

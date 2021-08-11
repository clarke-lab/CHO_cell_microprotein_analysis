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

## Preprocessing and mapping
```bash
# Trim adapters from the Cycoheximide, Harringtonine and No drug Ribosome footprint profiling data
./preprocess_rnaseq_pe.sh
./preprocess_riboseq.sh
```

## Offset determination
```

```

## Idenfication of novel CHO cell ORFs
```bash
./scripts/orfrater_analysis.sh
```

## Transcript level quantitation

```bash
```

# 2. Gene level quantitation

To enable the identification of differential translation between the NTS and TS conditions the mapped Ribo-seq and RNA-seq CDS counts are determined. For the Riboeq

```bash

```

# 1. Differential Translation Analysis
```bash
```
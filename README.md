# CHO Cell Ribo-seq 

The code contained in this repositority enable the reproduction of the results of

Castro-Rivadeneyra *et. al* 

The publication is freely availiable here:

## Preparation
### Dependencies
  * cutadapt 1.18
  * trimmomatic-0.36
  * STAR-2.7.8a
  * ORFATER
  * Docker
  * Plastid 
  * R
  * DESeq2

   
### Dowload the raw RiboSeq and RNASeq data

To be completed when data is uploaded to SRA and ENA

```bash
./scripts/get_raw_data.sh
```
### Reference genome

Download the PICR-H reference genome from NCBI and create a STAR index for mapping

```bash
mkdir -p reference_genome

# NCBI location
url=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/668/045/

# get the sequence and annotation
wget "$url"GCF_003668045.3_CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna.gz \
-P reference_genome

wget "$url"/GCF_003668045.3_CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf.gz \
-P reference_genome

# extract
gunzip reference_genome/*.gz
```

### Create reference index
A STAR index is created to map the Ribo-seq and RNA-seq data

```bash
# set the path
star_path=../bin/STAR-2.7.8a/bin/Linux_x86_64

#create a directory 
mkdir reference_genome/star_index_ncbi

# build the index
$star_path/STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --sjdbOverhang 124 \
     --genomeChrBinNbits 16 \
     --genomeDir reference_genome/star_index_ncbi \
     --genomeFastaFiles reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna \
     --sjdbGTFfile reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf
```

### Preprocessing and mapping
The raw sequencing data for the Ribo-seq and RNASeq is preprocessed to remove adapter sequence. In the case of the Ribo-seq data reads aligning to rRNA, snoRNA or tRNA are removed and only read lengths (28-31nt) where the expected peridocity is observed for 60% of reads were retained.

```bash
# Trim adapters from the Cycoheximide, Harringtonine and No drug Ribosome footprint profiling data
./preprocess_rnaseq_pe.sh
./preprocess_riboseq.sh
```

### Offset determination
Calculation of the P-site offset for RPFs for the merged and individual samples.  
```bash
./scripts/offset_determination.sh
```

## ORF-RATER analysis

We have built a docker image with ORF-RATER and required packages to ensure future compatability

Download the ORF-RATER docker image 
```bash
docker pull clarkelab/orfrater
```

### Run ORF-RATER 

The merged BAM files for Harringtone, cycloheximide and no-drug Ribo-seq as well as the Chinese hamster 
refercen annotation are inputted to ORF-RATER

```bash
./scripts/orfrater_analysis.sh
```
### ORF-RATER filtering

## Transcript level quantitation

```bash
```

## 2. Gene level quantitation

To enable the identification of differential translation between the NTS and TS conditions the mapped Ribo-seq and RNA-seq CDS counts are determined. For the Riboeq

```bash

```

## 1. Differential Translation Analysis
```bash
```

# Proteomics
Make a fasta file containing ORFRATER predictions and the current release of Uniprot proteins for the Chinese hamster 
```bash
./scripts/make_proteomics_fasta.sh
```
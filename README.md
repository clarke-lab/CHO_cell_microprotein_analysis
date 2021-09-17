# CHO Cell Ribo-seq 

The code contained in this repositority enable the reproduction of the results of

Castro-Rivadeneyra *et. al* 

The publication is freely availiable here:

## 1. Dependencies
  * cutadapt 1.18
  * trimmomatic-0.36
  * STAR-2.7.8a
  * ORFATER
  * Docker
  * Plastid 
  * R
  * DESeq2 
  * samtools


# Section 1: Analysis
## 2. Dowload the raw RiboSeq and RNASeq data

To be completed when data is uploaded to SRA and ENA

```bash
./scripts/get_raw_data.sh
```
## 3. Reference genome

Download the PICR-H reference genome from NCBI and create a STAR index for mapping

```bash
mkdir -p reference_genome

# NCBI location
url=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/668/045/

# get the sequence, annotation and feature table
wget "$url"GCF_003668045.3_CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna.gz \
-P reference_genome

wget "$url"/GCF_003668045.3_CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf.gz \
-P reference_genome

wget "$url"/GCF_003668045.3_CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_feature_table.txt \
-P reference_genome

# extract
gunzip reference_genome/*.gz
```

## 4. Create reference index
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

## 5. RPF and read mapping
This script preprocesses the raw sequencing data. For all data types the adapters are removed as well as low quality bases. 
For Ribo-seq data contaminating RNA species (rRNA, tRNA and snoRNA) are removed following mapping to individual indexes, 
remaining reads are filtered based on length with only those within the expected RPF range (28-31nt) retained. Finally the reads from all replicate 
```bash
./scripts/preprocess_reads.sh

# count the reads removed by filtering as well as the final RPFs
./scripts/
```

## 6. Finding the RPF Offset 
Calculation of the P-site offset and analysis of triplet periodicty for RPFs for the merged and individual samples.  
```bash
./scripts/p-site_identification.sh
```

## 7. ORF identification

### Get the docker image
We have built a docker image with ORF-RATER and required packages to ensure future compatability

Download the ORF-RATER docker image 
```bash
docker pull clarkelab/orfrater
```

### Identify CHO cell ORFs 

The merged BAM files for Harringtone, cycloheximide and no-drug Ribo-seq as well as the Chinese hamster 
refercen annotation are inputted to ORF-RATER

```bash
./scripts/identify_ORFs.sh
```
### Filter ORFs

Filter the ORF-RATER output to remove: 
  * ORFs < 5aa & ORFRATER score < 0.05
  * Truncations and Interal ORFs
  * When other ORFs overlap and have the same stop codon retain the longest

A list of ORFs in non-coding RNAs is created for downstream differential expression analysis

```
Rscript ./scripts/filter_ORFs.R
```

## 8. Amino acid sequences

### sORF amino acid usage
```

```

### Mass spectrometry fasta
Here we extract the amino acid sequences for short ORFs and combine with the Uniprot Chinese hamster proteins. A database can be created for Mass spec based HCP analysis
```
./scripts/create_ms_fasta.sh
```

## 9. Plastid reference 
Here a plastid reference is created to enable the determination of the RPKM of transcripts and gene-level counting. A mask is created to elminate the first 5 and last 15 codons of ORFs >100aa, for ORFs <= 100aa the first and last codons are exlcuded from the counting process
```
mkdir plastid
./scripts/make_plastid_reference.sh
```

## 9. Transcript level quantiation
The RPKM is calculated for each annotated transcript
```bash
./scripts/calculate_rpkm.sh
```

## 10. Gene level quantitation

To enable the identification of differential translation between the NTS and TS conditions the mapped Ribo-seq and RNA-seq CDS counts are determined. For the Riboeq

```bash
./scripts/calculate_gene_cds_counts.sh
```

## 11. Differential Translation Analysis
First we count the RPFs and RNA-seq reads mapping to CDS regions using Plastid
```bash
./scripts/calculate_gene_cds_count.sh
```

Then DESeq2 is using to calculate differential expression from the Ribo-seq and RNA-seq reads separately before differential translation is carried out.
```
Rscript ./scripts/run_deseq2.R
```

# Proteomics
Make a fasta file containing ORFRATER predictions and the current release of Uniprot proteins for the Chinese hamster 
```bash
./scripts/make_proteomics_fasta.sh
```

# Section 2: Reproduction of Tables and Figures

## Make alignment tracks
Here we make alignment tracks from the genome and transcriptome BAMs using individual replicates and the merged data
```
```


## Results section 2.1: Quality control 
Assessment of prerocessing and read length, phasing, metagene profiles for the Ribo-seq data
```
Rscript ./scripts/section_2.1.R
```

## Results section 2.2 ORF identification
Outputs of the ORF-RATER algorithm for the Chinese hamster genome
```
Rscript ./scripts/section_2.2.R
```

## Results section 2.3 Upstream ORF analysis
Analysis of the global effect of uORFs at the transcript level
```
Rscript ./scripts/section_2.3.R
```

## Results section 2.4 Differential Translation analysis of annotated protein coding genes
DESeq2 analysis of the NCBI annotated protein coding genes annotated in the Chinese hamster genome
```
Rscript ./scripts/section_2.4.R
```

## Results section 2.5 Differential Translation analysis of ORFs in non-coding RNA
DESeq2 analysis of ORFs identified in the Chinese hamster genome
```
Rscript ./scripts/section_2.5.R
```

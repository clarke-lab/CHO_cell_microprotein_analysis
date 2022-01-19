# Found in translation: Microproteins are a new class of host cell impurity in mAb drug products  

[![DOI](https://zenodo.org/badge/449655379.svg)](https://zenodo.org/badge/latestdoi/449655379)

The code contained in this repositority enable the reproduction of the results of:

Castro-Rivadeneyra *et. al* 2022. **Found in translation: Microproteins are a new class of host cell impurity in mAb drug products**

The publication is freely availiable here: xxxxxxx&nbsp;

**Abstract:**
<p style='text-align: justify;'>
Mass spectrometry (MS) has emerged as a powerful approach for the detection of Chinese hamster ovary (CHO) cell protein impurities in antibody drug products. The incomplete annotation of the Chinese hamster genome, however, limits the coverage of MS-based host cell protein (HCP) analysis.</p> &nbsp;

<p style='text-align: justify;'>
In this study, we performed ribosome footprint profiling (Ribo-seq) of translation initiation and elongation to refine the Chinese hamster genome annotation. Analysis of these data resulted in the identification of thousands of previously uncharacterised non-canonical proteoforms in CHO cells, such as N-terminally extended proteins and short open reading frames (sORFs) predicted to encode for microproteins. MS-based HCP analysis of adalimumab and trastuzumab with the extended protein sequence database, resulted in the detection of CHO cell microprotein impurities in mAb drug product for the first time. Further analysis revealed that the CHO cell microprotein population is altered over the course of cell culture and, in response to a change in cell culture temperature. The annotation of non-canonical Chinese hamster proteoforms permits a more comprehensive characterisation of HCPs in antibody drug products using MS.
</p>
&nbsp;

## Dependencies  

| Software | R packages      ||
| ------------- | --------------- | --------------- |
| [cutadapt 1.18](https://cutadapt.readthedocs.io/en/stable/)     | [tidyverse](https://tidyr.tidyverse.org) | [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html) |
| [STAR-2.7.8a](https://github.com/alexdobin/STAR) | [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) | [ggForce](https://ggforce.data-imaginist.com) |
| [trimmomatic-0.36](http://www.usadellab.org/cms/?page=trimmomatic) | [patchwork](https://patchwork.data-imaginist.com) | [BioStrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) |
| [Plastid](https://plastid.readthedocs.io/en/latest/) | [writexl](https://github.com/ropensci/writexl) | [readxl](https://readxl.tidyverse.org) |
| [ORF-RATER](https://github.com/alexfields/ORF-RATER) | [ggpp](https://cran.r-project.org/web/packages/ggpp/readme/README.html) | [ggpmisc](https://cran.r-project.org/web/packages/ggpmisc/index.html) |
| [Docker](https://www.docker.com/) | [wiggleplotr](https://bioconductor.org/packages/release/bioc/html/wiggleplotr.html) | [WebGestaltR](https://cran.r-project.org/web/packages/WebGestaltR/index.html) |
| [samtools](http://www.htslib.org/) | [GenomicFeatures](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html) | [heatmaply](https://cran.r-project.org/web/packages/heatmaply/index.html) |
| [Deeptools](https://deeptools.readthedocs.io/en/develop/) | [viridis](https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html) | [ggvenn](https://github.com/yanlinlin82/ggvenn) |
| [agat](https://github.com/NBISweden/AGAT) | [ggpubr](https://rpkgs.datanovia.com/ggpubr/) | [ggrepel](https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html) |
| [Kent Utilities](https://hgdownload.soe.ucsc.edu/admin/exe/) | [cowplot](https://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html) | [proDA](https://www.bioconductor.org/packages/release/bioc/html/proDA.html) |
|  | [scales](https://scales.r-lib.org) | [pheatmap](https://cran.r-project.org/web/packages/pheatmap/index.html) |
|  | [fuzzyjoin](https://cran.r-project.org/web/packages/fuzzyjoin/index.html) |

## Identification of CHO cell ORFs  

### 1. Download the raw Ribo-seq and RNA-seq data  

To be completed when data is uploaded to SRA and ENA  

```bash
# change to working directory
cd ribosome_footprint_profiling

# download data
./scripts/get_raw_data.sh
```  

### 2. Download the Chinese hamster genome

Download the PICR-H reference genome from NCBI and create a STAR index for mapping

```bash
# create a directory
mkdir -p reference_genome

# NCBI url
url=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/668/045/

# get the sequence, feature table and GTF,GFF annotation files
wget "$url"GCF_003668045.3_CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna.gz \
-P reference_genome

wget "$url"/GCF_003668045.3_CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf.gz \
-P reference_genome

wget "$url"/GCF_003668045.3_CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_feature_table.txt.gz \
-P reference_genome 

wget "$url"/GCF_003668045.3_CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gff.gz \
-P reference_genome 

# unzip
gunzip reference_genome/*.gz
```

### 3. Create read mapping index  

A STAR index is created to map the Ribo-seq and RNA-seq data

```bash
#create a directory 
mkdir reference_genome/star_index_ncbi

# set the path to STAR
star_path=../bin/STAR-2.7.8a/bin/Linux_x86_64

# build the index
$star_path/STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --sjdbOverhang 74 \
     --genomeChrBinNbits 16 \
     --genomeDir reference_genome/star_index_ncbi \
     --genomeFastaFiles reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna \
     --sjdbGTFfile reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf
```

### 4. RPF and read mapping  

<p style='text-align: justify;'>

This script preprocesses the raw sequencing data. For all data types the adapters are removed as well as low quality bases. For Ribo-seq data contaminating RNA species (rRNA, tRNA and snoRNA) are removed following mapping to individual indexes, remaining reads are filtered based on length with only those within the expected RPF range (28-31nt) retained. Finally the reads from all replicate.  
</p>

```bash
# preprocess
./scripts/preprocess_reads.sh

# count the reads removed by filtering as well as the final RPFs
./scripts/fastq_read_count.sh
```

### 5. Finding the RPF Offset  

Calculation of the P-site offset and analysis of triplet periodicty for RPFs for the merged and individual samples.  

```bash
./scripts/identify_RPF_psite.sh
```

### 6. ORF identification

We have built a docker image with ORF-RATER and required packages to ensure future compatability. The merged BAM files for Harringtone, cycloheximide and no-drug Ribo-seq as well as the Chinese hamster rare used with ORF-RATER to identify ORFs

```bash
# get docker image
docker pull clarkelab/orfrater:final

# run ORF-RATER
./scripts/identify_ORFs.sh
```  

Filter the ORF-RATER output to remove:

* ORFs < 5aa & ORFRATER score < 0.05

* Truncations and Interal ORFs

* When other ORFs overlap and have the same stop codon retain the longest

A list of ORFs in non-coding RNAs is created for downstream differential expression analysis

```bash
# create a directory to store ids for amino acid analysis and plastid quantitation
mkdir orf_lists 

# mkdir to store results
mkdir results/section_2.2

# filter the ORF-RATER output
Rscript ./scripts/filter_ORFs.R
```

## ORF analysis  

### 1. Amino acid usage

```bash
./scripts/get_ORF_amino_acid_sequences.sh
```

### 2. Preparing for HCP and proteomics analysis  

Here we extract the amino acid sequences for short ORFs and combine with the Uniprot Chinese hamster proteins. A database can be created for Mass spec based HCP analysis.  

```bash
mkdir proteomics_db

# create the protein sequence database
./scripts/create_ms_fasta.sh
```

## Quantitation  

### 1. Plastid reference  

<p style='text-align: justify;'>
Here a plastid reference is created to enable the determination of the RPKM of transcripts and gene-level counting. A mask is created to elminate the first 5 and last 15 codons of ORFs >100aa, for ORFs < 100aa the first and last codons are exlcuded from the counting process.  
</p>

```bash
mkdir plastid_reference

# make the plastid reference
./scripts/make_plastid_reference.sh
```

### 2. Transcript level quantiation

The RPKM is calculated for each annotated transcript

```bash
mkdir -p quantitation/transcript_cds_rpkm

# quantitate
./scripts/calculate_rpkm.sh
```

### 3. Gene level quantitation

To enable the identification of differential translation between the NTS and TS conditions the mapped Ribo-seq and RNA-seq CDS counts are determined. For the Riboeq

```bash
mkdir quantitation/gene_cds_counts

# count
./scripts/calculate_gene_cds_counts.sh
```

### 4. Differential Translation Analysis  

First we count the RPFs and RNA-seq reads mapping to CDS regions using Plastid

```bash
# count
./scripts/calculate_gene_cds_count.sh
```

Then DESeq2 is using to calculate differential expression from the Ribo-seq and RNA-seq reads separately before differential translation is carried out.

```bash
# the mouse annotation is used to replace missing gene CGR gene symbols
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_feature_table.txt.gz \
-P reference_genome

gunzip reference_genome/GCF_000001635.27_GRCm39_feature_table.txt.gz

# determine differential expression and translation
Rscript ./scripts/run_deseq2.R
```

### 5. Make alignment tracks  

Here we make the required alignment tracks for figures from the genome and transcriptome BAMs using individual replicates and the merged data

```bash
./scripts/make_coverage_tracks.sh
```

## Results

The following R notebooks allow the reproduction of figures and tables in the manuscript

### 2.1 Transcriptome wide analysis of CHO cell translation initiation and elongation using Ribo-seq

```bash
results/r_scripts/section_2_1.Rmd
```

### 2.2 Ribo-seq enables the characterisation of novel Chinese hamster proteoforms

Outputs of the ORF-RATER algorithm for the Chinese hamster genome

```bash
results/r_scripts/section_2_2.Rmd
```

### 2.3 The Chinese hamster genome harbours thousands of short open reading frames

Analysis of the global effect of uORFs at the transcript level

```bash
results/r_scripts/section_2_3.Rmd
```  

### 2.4 Detection of host cell microprotein contamination in adalimumab and trastuzamab drug products  

```bash
results/r_scripts/section_2_4.Rmd
```  

### 2.5 The translation efficiency of sORFs found in non-coding RNA genes is altered in response to mild hypothermia in CHO cells

DESeq2 analysis of ORFs identified in the Chinese hamster ncRNA  

```bash
results/r_scripts/section_2_5.Rmd
```  

### 2.6 Microproteins are differentially expressed between the exponential and stationary phases of CHO cell culture

```bash
results/r_scripts/section_2_6.Rmd
```

### Supplementaty Results

```bash
results/r_scripts/Supplementary_Results.Rmd
```

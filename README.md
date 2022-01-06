# Found in translation: Microproteins are a new class of host cell impurity in mAb drug products  

The code contained in this repositority enable the reproduction of the results of

Castro-Rivadeneyra *et. al* Found in translation: Microproteins are a new class of host cell impurity in mAb drug products

The publication is freely availiable here: xxxxxxx

**Abstract:**
Mass spectrometry (MS) has emerged as a powerful approach for the detection of Chinese hamster ovary (CHO) cell protein impurities in antibody drug products. The incomplete annotation of the Chinese hamster genome, however, limits the coverage of MS-based host cell protein (HCP) analysis.  
In this study, we performed ribosome footprint profiling (Ribo-seq) of translation initiation and elongation to refine the Chinese hamster genome annotation. Analysis of these data resulted in the identification of thousands of previously uncharacterised non-canonical proteoforms in CHO cells, such as N-terminally extended proteins and short open reading frames (sORFs) predicted to encode for microproteins. MS-based HCP analysis of adalimumab and trastuzumab with the extended protein sequence database, resulted in the detection of CHO cell microprotein impurities in mAb drug product for the first time. Further analysis revealed that the CHO cell microprotein population is altered over the course of cell culture and, in response to a change in cell culture temperature. The annotation of non-canonical Chinese hamster proteoforms permits a more comprehensive characterisation of HCPs in antibody drug products using MS.

## 1. Dependencies

To be added: viridis, ggpp, proDA, ggvenn

| Software | R packages |
| ------------- | --------------- |
| [cutadapt 1.18](https://cutadapt.readthedocs.io/en/stable/)     | [tidyverse](https://tidyr.tidyverse.org) |
| [STAR-2.7.8a](https://github.com/alexdobin/STAR) | [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) |
| [trimmomatic-0.36](http://www.usadellab.org/cms/?page=trimmomatic) | [patchwork](https://patchwork.data-imaginist.com) |
| [Plastid](https://plastid.readthedocs.io/en/latest/) | [writexl](https://github.com/ropensci/writexl) |
| [ORF-RATER](https://github.com/alexfields/ORF-RATER) | [ggpp](https://cran.r-project.org/web/packages/ggpp/readme/README.html) |
| [Docker](https://www.docker.com/) | [wiggleplotr](https://bioconductor.org/packages/release/bioc/html/wiggleplotr.html) |
| [samtools](http://www.htslib.org/) | [GenomicFeatures](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html) |
| [Deeptools](https://deeptools.readthedocs.io/en/develop/) | [viridis](https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html) |
| [agat](https://github.com/NBISweden/AGAT) | [ggpubr](https://rpkgs.datanovia.com/ggpubr/) |
| [Kent Utilities](https://hgdownload.soe.ucsc.edu/admin/exe/) | [cowplot](https://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html) |
|  | [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html) |
|  | [ggForce](https://ggforce.data-imaginist.com) |
|  | [BioStrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) |
|  | [readxl](https://readxl.tidyverse.org) |
|  | [ggpmisc]](https://cran.r-project.org/web/packages/ggpmisc/index.html) |
|  | [WebGestaltR]](https://cran.r-project.org/web/packages/WebGestaltR/index.html) |
|  | [heatmaply]](https://cran.r-project.org/web/packages/heatmaply/index.html) |
|  | [ggvenn]](https://github.com/yanlinlin82/ggvenn) |
|  | [ggrepel]](https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html) |
|  | [proDA]](https://www.bioconductor.org/packages/release/bioc/html/proDA.html) |
|  | [scales]](https://scales.r-lib.org) |
|  | [pheatmap]](https://cran.r-project.org/web/packages/pheatmap/index.html) |
|  | [fuzzyjoin]](https://cran.r-project.org/web/packages/fuzzyjoin/index.html) |

<br/>  

## Section 1: Analysis  


Make a directory for the analysis and enter it

```bash
mkdir ribosome_footprint_profiling
cd ribosome_footprint_profiling
```

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

wget "$url"/GCF_003668045.3_CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_feature_table.txt.gz \
-P reference_genome 

wget "$url"/GCF_003668045.3_CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gff.gz \
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
     --sjdbOverhang 74 \
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
./scripts/fastq_read_count.sh
```

## 6. Finding the RPF Offset 
Calculation of the P-site offset and analysis of triplet periodicty for RPFs for the merged and individual samples.  
```bash
./scripts/identify_RPF_psite.sh
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
# create a directory to store ids for amino acid analysis and plastid quantitation
mkdir orf_lists 

# mkdir to store results
mkdir results/section_2.2

# filter the ORF-RATER output
Rscript ./scripts/filter_ORFs.R
```

## 8. Amino acid sequences

### sORF amino acid usage
```
./scripts/get_ORF_amino_acid_sequences.sh
```

### Mass spectrometry fasta
Here we extract the amino acid sequences for short ORFs and combine with the Uniprot Chinese hamster proteins. A database can be created for Mass spec based HCP analysis
```
mkdir proteomics_db
./scripts/create_ms_fasta.sh
```

## 9. Plastid reference 
Here a plastid reference is created to enable the determination of the RPKM of transcripts and gene-level counting. A mask is created to elminate the first 5 and last 15 codons of ORFs >100aa, for ORFs <= 100aa the first and last codons are exlcuded from the counting process
```
mkdir plastid_reference
./scripts/make_plastid_reference.sh
```

## 9. Transcript level quantiation
The RPKM is calculated for each annotated transcript
```bash
mkdir -p quantitation/transcript_cds_rpkm
./scripts/calculate_rpkm.sh
```

## 10. Gene level quantitation

To enable the identification of differential translation between the NTS and TS conditions the mapped Ribo-seq and RNA-seq CDS counts are determined. For the Riboeq

```bash
mkdir quantitation/gene_cds_counts
./scripts/calculate_gene_cds_counts.sh
```

## 11. Differential Translation Analysis
First we count the RPFs and RNA-seq reads mapping to CDS regions using Plastid
```bash
./scripts/calculate_gene_cds_count.sh
```

Then DESeq2 is using to calculate differential expression from the Ribo-seq and RNA-seq reads separately before differential translation is carried out.
```
# the mouse annotation is used to replace missing gene CGR gene symbols
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_feature_table.txt.gz \
-P reference_genome

gunzip reference_genome/GCF_000001635.27_GRCm39_feature_table.txt.gz


Rscript ./scripts/run_deseq2.R
```

# Section 2: Reproduction of Tables and Figures

## Make alignment tracks
Here we make the required alignment tracks for figures from the genome and transcriptome BAMs using individual replicates and the merged data
```
./scripts/make_coverage_tracks.sh
```


## Quality control of RNA-seq and Ribo-seq data
Assessment of prerocessing and read length, phasing, metagene profiles for the Ribo-seq data
```
results/r_scripts/section_2_1.Rmd
```

## ORF identification 
Outputs of the ORF-RATER algorithm for the Chinese hamster genome
```
results/r_scripts/section_2_2.Rmd
```

## Upstream ORF analysis
Analysis of the global effect of uORFs at the transcript level
```
results/r_scripts/section_2_3.Rmd
```

## Differential Expression/Translation analysis for annotated protein coding genes
DESeq2 analysis of the NCBI annotated protein coding genes annotated in the Chinese hamster genome
```
results/r_scripts/section_2_4.Rmd
```

## Differential Expression/Translation analysis for ORFs in non-coding RNA
DESeq2 analysis of ORFs identified in the Chinese hamster ncRNA
```
results/r_scripts/section_2_5.Rmd
```

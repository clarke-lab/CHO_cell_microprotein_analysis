
# add start and stop codons
agat_sp_add_start_and_stop.pl \
-gff transcripts.fasta.transdecoder.genome.gff3 \
--fasta Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.fa \
-o transcripts.fasta.transdecoder.genome.stop.gff3

# convert to gtf
agat_convert_sp_gff2gtf.pl \
-gff transcripts.fasta.transdecoder.genome.stop.gff3 \
--gtf_version 3 \
-o transcripts.fasta.transdecoder.genome.gtf


# get the ENSEMBL transcripts found to be protein coding; remove the p* tag from transdecoder
awk '$3 == "transcript" {print $0}' transcripts.fasta.transdecoder.genome.gtf | \
awk -F'\'t '{print $9}' | awk -F';' '{print $2}' | sed s/transcript_id//g | \
sed s/\"//g | sed 's/[[:space:]]//g'| grep ENS | sed s/p..\*//g | sort | uniq | \
sed 's/\.//g' > reference_transcripts_in_transdecoder.txt

# get the stringtie transcripts found to be protein coding; remove the p* tag from transdecoder
awk '$3 == "transcript" {print $0}' transcripts.fasta.transdecoder.genome.gtf | \
awk -F'\'t '{print $9}' | awk -F';' '{print $2}' | sed s/transcript_id//g | \
sed s/\"//g | sed 's/[[:space:]]//g'| grep MSTRG | sed 's/...$//' > stringtie_transcripts_in_transdecoder.txt

# all protein coding transcripts identified by transdecoder
cat reference_transcripts_in_transdecoder.txt  stringtie_transcripts_in_transdecoder.txt \
> transcripts_in_transdecoder.txt

grep -Fv -f transcripts_in_transdecoder.txt cgr_stringtie_assembly_raw.gtf | awk '$3 == "transcript" {print $0}' | grep -v ENS | \
awk -F'\'t '{print $9}' | awk -F';' '{print $1}' | sed s/gene_id//g | \
sed s/\"//g | sed 's/[[:space:]]//g' | sort | uniq | awk '{print ""$1"\";"}' > non_coding_stringtie_genes.txt

# loop over stringtie gene ids and select the stringtie only genes
rm non_coding_stringtie_only.gtf
while read -ra a ; do
  echo ${a[0]}
  grep -w ${a[0]} cgr_stringtie_assembly_raw.gtf | grep -v ENS  | cat  >> non_coding_stringtie_only.gtf
done < non_coding_stringtie_genes.txt

# Select non-coding transcript IDs from the reference
grep -Fv -f reference_transcripts_in_transdecoder.txt Cricetulus_griseus_picr.CriGri-PICR.103.gtf | \
grep -v protein_coding | grep transcript | awk -F'\'t '{print $9}' | awk -F';' '{print $3}' | \
sed s/transcript_id//g | sed s/\"//g | sed 's/[[:space:]]//g' | sort | uniq > non_coding_reference_transcripts.txt

# Select non-coding transcrits annotated in ENSEMBL from stringtie
grep -F -f non_coding_reference_transcripts.txt cgr_stringtie_assembly_raw.gtf | awk -F'\'t '{print $9}' | \
awk -F';' '{print $1}' | sed s/gene_id//g | sed s/\"//g | sed 's/[[:space:]]//g'| \
sort | uniq | awk '{print "gene_id "$1"\";"}' > non_coding_stringtie_genes.txt

grep -F -f non_coding_reference_transcripts.txt cgr_stringtie_assembly_raw.gtf > non_coding_reference_genome.gtf

  agat_sp_merge_annotations.pl \
  -f non_coding_reference_genome.gtf \
  -f non_coding_stringtie_only.gtf \
  -f transcripts.fasta.transdecoder.genome.stop.gff3 \
  -o completed_annotation.gtf

  agat_convert_sp_gff2gtf.pl \
  -gff completed_annotation.gtf.gff \
  --gtf_version 2.5 \
  -o completed_annotation_final.gtf

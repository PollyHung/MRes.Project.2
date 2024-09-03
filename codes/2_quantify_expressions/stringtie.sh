#!/bin/bash
#PBS -l select=1:ncpus=8:mem=100gb:ompthreads=8
#PBS -l walltime=9:00:00          
#PBS -N stringtie

# To assess the support for LR-seq isoforms by short-read sequencing, we 
# performed RNA-seq and quantified isoform expression in 29 of our 30 
# LR-seq profiled breast samples. Briefly, 76–base pair long paired-end 
# RNA-seq libraries were sequenced at an average depth of 46 million reads 
# per sample and mapped to our LR-seq breast cancer transcriptome using hisat2 
# and quantified using StringTie. While 89% of the annotated isoforms (FSM) 
# were detected by RNA-seq [FPKM (fragments per kilobase per million mapped 
# reads) > 0.5], NIC and NNC isoforms have average detection rates of 
# 62 and 41%, respectively (Fig. 1E).

## Directories 
HOME="/rds/general/user/ph323"
REFERENCE="$HOME/home/MRes.project.2/docs/0_preprocessing"
DATA="$HOME/ephemeral/MRes.project.2/raw_data"
SHORT_READ="$DATA/short_reads/RNA-seq/G72_run990"
SQANTI="$DATA/tappAS/sqanti"

## Packages 
HISAT2="$HOME/home/MRes.project.2/codes/0_packages/hisat2"
STRINGTIE="$HOME/home/MRes.project.2/codes/0_packages/stringtie-2.2.3"
module load anaconda3/personal
module load samtools
source activate Renv

## Files 
hg38_fasta="$REFERENCE/hg38_fasta/hg38_customized_noChr/hg38_customized.fasta"
hg38_hisat2="$REFERENCE/hg38_fasta/hg38_hisat2_index"
short_read_id="$HOME/home/MRes.project.2/codes/0_sample_id/short_reads/short_reads_id.txt"
sr_id="$HOME/home/MRes.project.2/codes/0_sample_id/short_reads/sr_id.txt"
batch="$HOME/home/MRes.project.2/codes/0_sample_id/batch/batch_9.txt"

## Code 
# cd "$HISAT2"
# ./hisat2-build "$hg38_fasta" "$hg38_hisat2"

while read -r sr_id; do 
  echo "working on $sr_id"

  ## perform hisat2 ------------------------------------------------------------
  cd "$SHORT_READ/$sr_id"
  read_1=$(ls *R1_001.fastq.gz)
  read_2=$(ls *R2_001.fastq.gz)
  R1_fq="$SHORT_READ/$sr_id/$read_1"
  R2_fq="$SHORT_READ/$sr_id/$read_2"
  output="$SHORT_READ/$sr_id/${sr_id}_aligned.sam"
  # cd "$HISAT2"
  # ./hisat2 --dta --threads 8 -x "$hg38_hisat2" -1 "$R1_fq" -2 "$R2_fq" -S "$output"
  
  ## clean up ------------------------------------------------------------------
  cd "$SHORT_READ/$sr_id"
  bam="$SHORT_READ/$sr_id/${sr_id}_aligned.bam"
  sorted_bam="$SHORT_READ/$sr_id/${sr_id}_sorted.aligned.bam"
  fixmate_bam="$SHORT_READ/$sr_id/${sr_id}_sorted.fixmate.aligned.bam"
  nodup_bam="$SHORT_READ/$sr_id/${sr_id}_sorted.nodup.aligned.bam"
  # samtools view -S -b -@ 7 -o "$bam" "$output"
  # samtools sort -n -@ 7 -o "$sorted_bam" "$bam"
  # samtools fixmate -m -@ 7 "$sorted_bam" "$fixmate_bam"
  # samtools sort -@ 7 -o "$sorted_bam" "$fixmate_bam"
  # samtools markdup -r -@ 7 "$sorted_bam" "$nodup_bam"
  # rm "$output" 
  
  ## stringtie -----------------------------------------------------------------
  cd $STRINGTIE
  sqanti_gtf="$SQANTI/qc_corrected.gtf"
  sqanti_nochr_gtf="$SQANTI/qc_corrected_nochr.gtf"
  output_gtf="$SHORT_READ/$sr_id/${sr_id}_stringtie.gtf"
  gene_abundance="$SHORT_READ/$sr_id/${sr_id}_abundance.tsv"
  ./stringtie -p 8 -e -G "$sqanti_nochr_gtf" -o "$output_gtf" -A "$gene_abundance" "$nodup_bam" # estimate the pacbio transcript expression

done < "$short_read_id"









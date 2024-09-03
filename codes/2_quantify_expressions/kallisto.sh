#!/bin/bash

#PBS -l select=1:ncpus=100:mem=100gb
#PBS -l walltime=8:00:00          
#PBS -N kallisto

## establishing environment ----------------------------------------------------
cd /rds/general/user/ph323/home
source ~/.bashrc
conda activate base


## meta directories ------------------------------------------------------------
HOME="/rds/general/user/ph323/home"
PACKAGE="$HOME/MRes.project.2/codes/0_packages"
ANACONDA="$HOME/anaconda3/bin"
DOCS="$HOME/MRes.project.2/docs/0_preprocessing"
SHORT_READS="$DOCS/short_reads"
DATA_DIR="/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/short_reads/RNA-seq/G72_run990"
tappAS="/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/tappAS"

## functions and packages 
kallisto="$ANACONDA/kallisto" #kallisto 0.50.1


## files -----------------------------------------------------------------------
# references 
#gencode_fasta="$DOCS/gencode.v39.transcripts.fa"
#gencode_index="$DOCS/gencode_v39_transcriptome.idx"
transcript_fa="$DOCS/hg38.fa"
transcript_idx="$DOCS/hg38.idx"
short_read_fofn="$DOCS/pooled_short_reads.fofn"
short_reads_id="$HOME/MRes.project.2/codes/0_sample_id/short_reads_id.txt"


## Running kallisto ------------------------------------------------------------
# getting all the short read ids 
#cd $DATA_DIR
#ls -d * > $short_reads_id

# Building a public short read database 
#$kallisto index -i $gencode_index $gencode_fasta
$kallisto index -i $transcript_idx $transcript_fa

# Quantify Transcript Abundance for each sample 
while read -r sample_id; do

  start_time=$(date)
  TEMP_DIR="$DATA_DIR/$sample_id"
  short_read_R1=$(find "$TEMP_DIR" -type f -name "*R1_001.fastq")
  short_read_R2=$(find "$TEMP_DIR" -type f -name "*R2_001.fastq")
  cd $TEMP_DIR
  
  $kallisto quant \
    -i $gencode_index \
    -o "$transcription_expression/$sample_id" \
    -b 100 $short_read_R1 $short_read_R2 \
    -t 100
  end_time=$(date)
  echo "Finish processing $sample_id between $start_time and $end_time"
  
done < $short_reads_id










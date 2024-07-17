#!/bin/bash

#PBS -l select=1:ncpus=256:mem=300gb
#PBS -l walltime=4:00:00          
#PBS -N deeploc2

## meta directories ------------------------------------------------------------
HOME="/rds/general/user/ph323/home"
PACKAGE="$HOME/MRes.project.2/codes/0_packages"
ANACONDA="$HOME/anaconda3/"
DOCS="$HOME/MRes.project.2/docs/0_preprocessing"
deeploc="$ANACONDA/envs/Renv/bin/deeploc2"

## Uniprot ---------------------------------------------------------------------
uniprot_fasta="$DOCS/uniprot/human_fasta"
chunk_fasta="$uniprot_fasta/human_chunk_28.fasta"
$deeploc -f $chunk_fasta -o "$DOCS/uniprot/deeploc_split" -m "Accurate"


## loop over each sample -------------------------------------------------------
while read -r sample_id; do
  sample_id="T16-088_FT2"
  start_time=$(date)
  DATA_DIR="/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/raw/$sample_id"
  PROCESSED_DATA="$DATA_DIR/processed_data"
  SQANTI="$DATA_DIR/sqanti_output"
  TRANSDECODER="$DATA_DIR/transcoder"
  DEEPLOC="$DATA_DIR/deeploc"
  if [ ! -d "$DEEPLOC" ]; then 
    mkdir -p "$DEEPLOC"
  fi
  cd $DEEPLOC
  old="$TRANSDECODER/qc_corrected.fasta.transdecoder.pep"
  Amino_acids="$TRANSDECODER/qc_corrected_transdecoder.fasta"
  if [ ! -f "$Amino_acids" ]; then
    cp "$old" "$Amino_acids"
  fi
  cd "$PACKAGE/deeploc2"
  deeploc2 -f $Amino_acids -o "$DEEPLOC" -m "Fast"
  end_time=$(date)
  echo "Performed DeepLoc2.0 for $sample_id from $start_time to $end_time"
done < "/rds/general/user/ph323/home/MRes.project.2/codes/0_sample_id/sample_list.txt"









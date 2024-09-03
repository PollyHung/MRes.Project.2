#!/bin/bash

#PBS -l select=1:ncpus=70:mem=100gb
#PBS -l walltime=4:00:00          
#PBS -N flair

## directories 
CODE="/rds/general/user/ph323/home/MRes.project.2/codes"
PACKAGE="/rds/general/user/ph323/home/MRes.project.2/codes/0_packages"
DOCS="/rds/general/user/ph323/home/MRes.project.2/docs/0_preprocessing"
DATA_DIR="/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/raw"
tappAS="/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/tappAS"

## install packages 
#cd "$PACKAGE"
#conda create -n flair -c conda-forge -c bioconda flair

## files 
manifest="$CODE/0_sample_id/flair_manifest.txt"
fasta="$tappAS/sqanti3/filter_rule.filtered.fasta"

## activate flair and quantify 
cd "$PACKAGE"
source ~/.bashrc
conda activate flair
flair quantify \
  --reads_manifest "$manifest" \
  --isoforms "$fasta" \
  --output "$tappAS/flair" \
  --threads 64 \
  --temp_dir "$tappAS/flair/temp" \
  --sample_id_only \
  --trust_ends 
  
 












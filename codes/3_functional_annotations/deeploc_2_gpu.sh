#!/bin/bash
#PBS -l select=1:ncpus=4:mem=300gb:ngpus=1
#PBS -l walltime=3:00:00          
#PBS -N chunk_6

eval "$(~/miniconda3/bin/conda shell.bash hook)"
conda  activate deeploc2

start=$(date)
tappAS="/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/tappAS"
CHUNKS="$tappAS/transcoder/chunks"
  
chunk_fasta="$CHUNKS/chunk_6.fasta"
deeploc2 -f $chunk_fasta -o "$tappAS/deeploc_gpu" -m "Fast"
end=$(date)
echo "deeploc2 for $chunk_fasta done between $start and $end"


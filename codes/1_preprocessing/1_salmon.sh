#!/bin/bash

#PBS -l select=1:ncpus=256:mem=100gb
#PBS -l walltime=8:00:00          
#PBS -N salmon

# install package 
#conda install bioconda::salmon

echo "load conda environment"
module load anaconda3/personal
source activate Renv

echo "set up directories and functions"
# directory 
PACKAGE="/rds/general/user/ph323/home/MRes.project.2/codes/0_packages"
ANACONDA="/rds/general/user/ph323/home/anaconda3/envs/Renv/bin"
tappAS="/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/tappAS"
salmon_dir="$tappAS/salmon"

# functions 
salmon="$ANACONDA/salmon"
fasta="$tappAS/sqanti3/qc_corrected.fasta"
short_read_list="/rds/general/user/ph323/home/MRes.project.2/codes/0_sample_id/short_reads_path.txt"

# Star Indexing and Alignment --------------------------------------------------
#$salmon index -t "$fasta" -i "$salmon_dir/transcript_index"

while IFS= read -r line
do
  # Split the line into R1 and R2
  set -- $line
  R1=$1
  R2=$2

  # Extract sample and pair names from the file path
  shortRead_id=$(basename $R1 | cut -d'_' -f1)
  echo "processing salmon for $shortRead_id"
  
  # Output directory for this sample pair
  output_dir="$salmon_dir/output/$shortRead_id"
  mkdir -p $output_dir
  
  # Run Salmon quant
  $salmon quant -i "$salmon_dir/transcript_index" \
    --libType A \
    -1 $R1 \
    -2 $R2 \
    --validateMappings \
    --threads 8 \
    -o $output_dir \
    --useEM \
    --numBootstraps 10 
  
done < "$short_read_list" 









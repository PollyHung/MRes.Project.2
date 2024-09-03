#!/bin/bash

#PBS -l select=1:ncpus=10:mem=100gb
#PBS -l walltime=8:00:00          
#PBS -N isoAnnotLite

## start time?
start_time=$(date)
echo "Start performing IsoAnnotLite for $sample_id at $start_time"

module load anaconda3/personal

## directories =================================================================
PACKAGE="/rds/general/user/ph323/home/MRes.project.2/codes/0_packages"
DOCS="/rds/general/user/ph323/home/MRes.project.2/docs/0_preprocessing"
DATA_DIR="/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/raw/$sample_id"
PROCESSED_DATA="$DATA_DIR/processed_data"
SQANTI_OUTPUT="$DATA_DIR/sqanti_output"
if [ ! -d "$SQANTI_OUTPUT" ]; then 
  mkdir -p "$SQANTI_OUTPUT"
fi

## packages ====================================================================
SQANTI="$PACKAGE/SQANTI3-5.2.1"

## functions ===================================================================
IsoAnnotLite="$SQANTI/utilities/IsoAnnotLite_SQ3.py"

## Loading environment =========================================================
# activate SQANTI3 
cd "$SQANTI"
source ~/.bashrc
conda activate SQANTI3.env
# change to data directory 
cd "$SQANTI_OUTPUT"

## files ==============================
corrected_gtf="$SQANTI_OUTPUT/qc_corrected.gtf"
classification="$SQANTI_OUTPUT/qc_classification.txt"
junctions="$SQANTI_OUTPUT/qc_junctions.txt"
ref_gff="$DOCS/human_tappas_gencode_annotation_file.gff3"

python $IsoAnnotLite "$corrected_gtf" "$classification" "$junctions" \
  -gff3 "$ref_gff" \
  -o "isoAnnot" \
  -novel \
  -stdout "stats"


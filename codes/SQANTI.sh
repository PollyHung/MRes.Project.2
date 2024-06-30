#!/bin/bash

#PBS -l select=1:ncpus=10:mem=100gb
#PBS -l walltime=8:00:00          
#PBS -N preprocess

## start time?
start_time=$(date)
echo "Start performing SQANTI for $sample_id at $start_time"

module load anaconda3/personal

sample_id="090061A"

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
quality_control="$SQANTI/sqanti3_qc.py"
filter="$SQANTI/sqanti3_filter.py"
rescue="$SQANTI/sqanti3_rescue.py"

## Loading environment =========================================================
# activate SQANTI3 
cd "$SQANTI"
source ~/.bashrc
conda activate SQANTI3.env
# change to data directory 
cd "$SQANTI_OUTPUT"

## files =======================================================================
# sequencing files 
isoforms="$PROCESSED_DATA/collapsed.gff"
abundance="$PROCESSED_DATA/collapsed.abundance.txt"
classification="$SQANTI_OUTPUT/qc_classification.txt"
rule_filtered_classification="$SQANTI_OUTPUT/filter_rule_RulesFilter_result_classification.txt"
filtered_corrected_classification="$SQANTI_OUTPUT/filter_rule_RulesFilter_result_classification_corrected.txt"
qc_corrected="$SQANTI_OUTPUT/qc_corrected.gtf"
qc_faa="$SQANTI_OUTPUT/qc_corrected.faa"
isoform_corrected_fasta="$SQANTI_OUTPUT/qc_corrected.fasta"
filtered_gtf="$SQANTI_OUTPUT/filter_rule.filtered.gtf"

# reference files 
annotation="$DOCS/gencode.v46.chr_patch_hapl_scaff.annotation.sorted.gtf"
genome="$DOCS/hg38.fa"
polyAmotifs="$DOCS/mouse_and_human.polyA_motif.txt"
polyApeaks="$DOCS/atlas.clusters.2.0.GRCh38.96.bed"
cage="$DOCS/human.refTSS_v3.1.hg38.bed"
isoAnnotGFF="$DOCS/gencode.v46.chr_patch_hapl_scaff.annotation.gff3"
intropolis="$DOCS/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified.gz"
rules="$DOCS/filter_custom.json"

## Start doing SQANTI quality control ==========================================
python $quality_control "$isoforms" "$annotation" "$genome" \
  --dir "$SQANTI_OUTPUT" \
  --fl_count "$abundance" \
  --polyA_motif_list "$polyAmotifs" \
  --polyA_peak "$polyApeaks" \
  --CAGE_peak "$cage" \
  --aligner_choice "minimap2" \
  --cpus 4 \
  --report both \
  --output "qc"
  
## Start doing SQANTI filter ===================================================
# section 1: rules filtering, using the customized rule adapted from Veiga et al., (2022), 
# Note that minimum coverage was not included in exclusion parameter because short-read
# was not provided in sqanti_qc.py step. 
python $filter rules "$classification" \
  --output "filter_rule" \
  --dir "$SQANTI_OUTPUT" \
  --json_filter "$rules" \
  --gtf "$qc_corrected" \
  --faa "$qc_faa"
  
## section 2: run ML filtering: 
#Error in `dplyr::mutate()`:
#ℹ In argument: `variable = factor(variable) %>% forcats::fct_reorder(importance)`.
#Caused by error:
#! object 'variable' not found
#https://github.com/ConesaLab/SQANTI3/issues/290

awk -F'\t' 'BEGIN {OFS="\t"} {if ($49 == "Isoform") $49 = "isoform"; print}' "$rule_filtered_classification" > temp_file && mv temp_file "$filtered_corrected_classification"


## Start doing SQANTI rescue ===================================================
python $rescue rules "$filtered_corrected_classification" \
  --isoforms "$isoform_corrected_fasta" \
  --gtf "$filtered_gtf" \
  -g "$annotation" \
  -f "$genome" \
  -k "$classification" \
  -d "$SQANTI_OUTPUT" \
  -j "$rules" \
  -o "rescue" \
  --mode automatic \
  -e all
  

end_time=$(date)
echo "Finishing performing SQANTI for $sample_id at $end_time"
echo -e "\a"

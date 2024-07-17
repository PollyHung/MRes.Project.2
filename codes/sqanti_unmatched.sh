#!/bin/bash

#PBS -l select=1:ncpus=256:mem=400gb
#PBS -l walltime=8:00:00          
#PBS -N T16-035-FN2

module load anaconda3/personal

#while read -r sample_id; do
  
  sample_id="T16-035-FN2"
  ## start time?
  start_time=$(date)
  
  ## directories =================================================================
  PACKAGE="/rds/general/user/ph323/home/MRes.project.2/codes/0_packages"
  DOCS="/rds/general/user/ph323/home/MRes.project.2/docs/0_preprocessing"
  DATA_DIR="/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/raw/$sample_id"
  ISOSEQ3="$DATA_DIR/isoseq3"
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
  cd "$SQANTI_OUTPUT"
  
  ## files =======================================================================
  # IsoSeq3 ----------------------------------------------------------------------
  isoforms="$ISOSEQ3/collapsed.gff"
  abundance="$ISOSEQ3/collapsed.abundance.txt"
  
  # Short Reads ------------------------------------------------------------------
  short_read_fofn="$DOCS/pooled_short_reads.fofn"
  
  # SQANTI3 quality control ------------------------------------------------------
  classification="$SQANTI_OUTPUT/qc_classification.txt"
  qc_corrected="$SQANTI_OUTPUT/qc_corrected.gtf"
  qc_faa="$SQANTI_OUTPUT/qc_corrected.faa"
  isoform_corrected_fasta="$SQANTI_OUTPUT/qc_corrected.fasta"
  
  # SQANTI3 filtering ------------------------------------------------------------
  rule_filtered_classification="$SQANTI_OUTPUT/filter_rule_RulesFilter_result_classification.txt"
  filtered_corrected_classification="$SQANTI_OUTPUT/filter_rule_RulesFilter_result_classification_corrected.txt"
  filtered_gtf="$SQANTI_OUTPUT/filter_rule.filtered.gtf"
  filtered_fasta="$SQANTI_OUTPUT/filter_rule.filtered.fasta"
  
  # reference files 
  annotation="$DOCS/gencode.v39.chr_patch_hapl_scaff.annotation.sorted.gtf"
  genome="$DOCS/hg38.fa"
  polyAmotifs="$DOCS/mouse_and_human.polyA_motif.txt"
  polyApeaks="$DOCS/atlas.clusters.2.0.GRCh38.96.bed"
  cage="$DOCS/human.refTSS_v3.1.hg38.bed"
  isoAnnotGFF="$DOCS/human_tappas_gencode_annotation_file.gff3"
  intropolis="$DOCS/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified.gz"
  rules="$DOCS/filter_custom.json"

  ## Start doing SQANTI quality control ==========================================
  python $quality_control "$isoforms" "$annotation" "$genome" \
    --dir "$SQANTI_OUTPUT" \
    --fl_count "$abundance" \
    --short_reads "$short_read_fofn" \
    --polyA_motif_list "$polyAmotifs" \
    --polyA_peak "$polyApeaks" \
    --CAGE_peak "$cage" \
    --isoAnnotLite \
    --gff3 "$isoAnnotGFF" \
    --aligner_choice "minimap2" \
    --cpus 256 \
    --report both \
    --output "qc"
    
  ## Start doing SQANTI filter ===================================================
  python $filter rules "$classification" \
    --output "filter_rule" \
    --dir "$SQANTI_OUTPUT" \
    --json_filter "$rules" \
    --gtf "$qc_corrected" \
    --faa "$qc_faa" \
    --isoforms "$isoform_corrected_fasta"
  
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

#done < "/rds/general/user/ph323/home/MRes.project.2/docs/0_preprocessing/unmatched_samples.txt"




#!/bin/bash

#PBS -l select=1:ncpus=64:mem=400gb
#PBS -l walltime=3:00:00          
#PBS -N preprocess

# Read Set_6.txt line by line
while read -r sample_id; do

  ## start time?
  start_time=$(date)

  ## directories =================================================================
  PACKAGE="/rds/general/user/ph323/home/MRes.project.2/codes/0_packages"
  DOCS="/rds/general/user/ph323/home/MRes.project.2/docs/0_preprocessing"
  DATA_DIR="/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/raw/$sample_id"
  PROCESSED_DATA="$DATA_DIR/processed_data/"
  REPORT_FILES="$DATA_DIR/report_files/"
  
  ## packages ====================================================================
  SMRT_LINK="$PACKAGE/smrtlink_v13/smrtcmds/bin"
  scISA="$PACKAGE/scISA-Tools/bin"
  BGI="$PACKAGE/bgi_commands/bgi_pipeline/bin"
  ANACONDA="/rds/general/user/ph323/home/anaconda3/bin/"
  module load blast+/2.11.1
  
  ## functions ===================================================================
  ccs="$SMRT_LINK/ccs"
  isoseq="$SMRT_LINK/isoseq"
  pbmm2="$SMRT_LINK/pbmm2"
  pigeon="$SMRT_LINK/pigeon"
  samtools="$SMRT_LINK/samtools"
  classify_by_primer="$BGI/classify_by_primer.pl"
  
  ## files =======================================================================
  # reference files 
  primer="/rds/general/user/ph323/home/MRes.project.2/codes/1_preprocessing/primer.fasta"
  ref_fa="$DOCS/hg38.fa"
  ref_gtf="$DOCS/gencode.v46.annotation.gtf"
  sorted_gtf="$DOCS/gencode.v46.annotation.sorted.gtf"
  
  # sequencing files 
  ccs_bam="$PROCESSED_DATA/ccs.bam"
  ccs_sam="$PROCESSED_DATA/ccs.sam"
  ccs_fasta="$PROCESSED_DATA/ccs.fa"
  mapped_m7="$PROCESSED_DATA/mapped.m7"
  flnc_bam="$PROCESSED_DATA/isoseq_flnc.bam"
  clustered_bam="$PROCESSED_DATA/clustered.bam"
  mapped_bam="$PROCESSED_DATA/mapped.bam"
  collapsed_gff="$PROCESSED_DATA/collapsed.gff"
  collapsed_sorted_gff="$PROCESSED_DATA/collapsed.sorted.gff"
  flnc_count="$PROCESSED_DATA/collapsed.flnc_count.txt"
  classification="$PROCESSED_DATA/collapsed_classification.txt"
  classification_lite="$PROCESSED_DATA/collapsed_classification.filtered_lite_classification.txt"
  saturation="$PROCESSED_DATA/saturation.txt"
  
  # reports 
  ccs_report="$REPORT_FILES/ccs_report.csv"
  ccs_log="$REPORT_FILES/ccs_log.txt"
  isoseq2_cluster2_log="$REPORT_FILES/isoseq2_cluster2_log.txt"  
  echo "All file path pointed, start processing sample $sample_id"
  
  ## running pipeline ============================================================
  # step 1: from subreads.bam to ccs.bam  
  if [ ! -f "$ccs_bam" ]; then
    echo "ccs.bam not found, start making ccs.bam"
    for subreads in "$DATA_DIR"/*.bam; do                                           
      $ccs "$subreads" "$ccs_bam" \
        --min-passes 0 \
        --min-length 50 \
        --max-length 21000 \
        --min-rq 0.75 \
        --num-threads 60 \
        --log-level INFO \
        --reportFile "$ccs_report" \
        --logFile "$ccs_log"
    done
  else
    echo "ccs.bam already exists! Skipping to next step."
  fi
  
  # step 2: classify CCS by primer blast
  echo "step 2: classify ccs by primer blast!"
  $samtools view "$ccs_bam" > "$ccs_sam"
  $samtools view "$ccs_bam" | awk '{print ">"$1"\n"$10}' > "$ccs_fasta"
  makeblastdb -in "$primer" -dbtype nucl
  blastn -query "$ccs_fasta" -db "$primer" -outfmt 7 -word_size 5 > "$mapped_m7"
  chmod +x "$mapped_m7"
  chmod +x "$ccs_fasta"
  perl "$classify_by_primer" \
    -blastm7 "$mapped_m7" \
    -ccsfa "$ccs_fasta" \
    -outdir "$PROCESSED_DATA"\
    -umilen 8 \
    -min_primerlen 16 \
    -min_isolen 200
  
  # step 3: clustering
  echo "step 3: clustering and collapsing!"
  $isoseq cluster2 "$flnc_bam" "$clustered_bam"
  $pbmm2 align --preset ISOSEQ --sort "$clustered_bam" "$ref_fa" "$mapped_bam"
  $isoseq collapse --do-not-collapse-extra-5exons "$mapped_bam" "$collapsed_gff"
  
  # step 4: prepare refernce files for pigeon 
  echo "step 4: pigeon to produce gff!"
  $pigeon prepare "$ref_gtf"
  $pigeon prepare "$collapsed_gff"
  $pigeon classify "$collapsed_sorted_gff" "$sorted_gtf" "$ref_fa" --fl "$flnc_count"
  $pigeon filter "$classification" --isoforms "$collapsed_sorted_gff"
  $pigeon report "$classification_lite" "$saturation"

  ## finish time
  end_time=$(date)
  echo "Finish performing ccs_bam conversion for $sample_id at $end_time"
  
done < "/rds/general/user/ph323/home/MRes.project.2/codes/1_preprocessing/sample_list.txt"






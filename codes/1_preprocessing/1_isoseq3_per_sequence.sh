#!/bin/bash

#PBS -l select=1:ncpus=10:mem=100gb
#PBS -l walltime=8:00:00          
#PBS -N isoseq3

module load blast+/2.11.1

while read -r sample_id; do

  ## start time?
  start_time=$(date)

  ## directories =================================================================
  PACKAGE="/rds/general/user/ph323/home/MRes.project.2/codes/0_packages"
  DOCS="/rds/general/user/ph323/home/MRes.project.2/docs/0_preprocessing"
  DATA_DIR="/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/raw/$sample_id"
  ISOSEQ3="$DATA_DIR/isoseq3"
  REPORT_FILES="$DATA_DIR/report_files"
  
  ## packages ====================================================================
  SMRT_LINK="$PACKAGE/smrtlink_v13/smrtcmds/bin"
  scISA="$PACKAGE/scISA-Tools/bin"
  BGI="$PACKAGE/bgi_commands/bgi_pipeline/bin"
  ANACONDA="/rds/general/user/ph323/home/anaconda3/bin"
  
  ## functions ===================================================================
  ccs="$SMRT_LINK/ccs"
  isoseq="$SMRT_LINK/isoseq"
  pbmm2="$SMRT_LINK/pbmm2"
  pigeon="$SMRT_LINK/pigeon"
  samtools="$SMRT_LINK/samtools"
  classify_by_primer="$BGI/classify_by_primer.pl"
  flnc2sam="$BGI/flnc2sam.pl"
  
  ## files =======================================================================
  # reference files 
  primer="/rds/general/user/ph323/home/MRes.project.2/codes/1_preprocessing/primer.fasta"
  ref_fa="$DOCS/hg38.fa"
  ref_gtf="$DOCS/gencode.v39.chr_patch_hapl_scaff.annotation.gtf"
  sorted_gtf="$DOCS/gencode.v39.chr_patch_hapl_scaff.annotation.sorted.gtf"
  
  # sequencing files 
  ccs_bam="$ISOSEQ3/ccs.bam"
  ccs_sam="$ISOSEQ3/ccs.sam"
  ccs_fasta="$ISOSEQ3/ccs.fa"
  mapped_m7="$ISOSEQ3/mapped.m7"
  flnc_fasta="$ISOSEQ3/isoseq_flnc.fasta"
  flnc_sam="$ISOSEQ3/isoseq_flnc.sam"
  flnc_bam="$ISOSEQ3/isoseq_flnc.bam"
  clustered_bam="$ISOSEQ3/clustered.bam"
  mapped_bam="$ISOSEQ3/mapped.bam"
  collapsed_gff="$ISOSEQ3/collapsed.gff"
  collapsed_sorted_gff="$ISOSEQ3/collapsed.sorted.gff"
  flnc_count="$ISOSEQ3/collapsed.flnc_count.txt"
  classification="$ISOSEQ3/collapsed_classification.txt"
  classification_lite="$ISOSEQ3/collapsed_classification.filtered_lite_classification.txt"
  saturation="$ISOSEQ3/saturation.txt"
  
  # reports 
  ccs_report="$REPORT_FILES/ccs_report.csv"
  ccs_log="$REPORT_FILES/ccs_log.txt"
  isoseq2_cluster2_log="$REPORT_FILES/isoseq2_cluster2_log.txt"
  
  cd $ISOSEQ3
  echo "Currently at: $ISOSEQ3"

  ## running pipeline ============================================================
  # step 1: from subreads.bam to ccs.bam  
  if [ ! -f "$collapsed_gff" ]; then
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
      -outdir "$ISOSEQ3"\
      -umilen 8 \
      -min_primerlen 16 \
      -min_isolen 200
    
    # step 3: clustering
    echo "step 3: clustering and collapsing!"
    $flnc2sam "$ccs_sam" "$flnc_fasta" > "$flnc_sam"
    samtools view -bS "$flnc_sam" > "$flnc_bam"
  
    else
    echo "$collapsed_gff already exists! Skipping to next step."
  fi
  
done < "/rds/general/user/ph323/home/MRes.project.2/codes/1_preprocessing/sample_list.txt"



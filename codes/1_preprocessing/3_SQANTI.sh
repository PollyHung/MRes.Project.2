#!/bin/bash

#PBS -l select=1:ncpus=32:mem=100gb
#PBS -l walltime=24:00:00          
#PBS -N SQANTI3

module load anaconda3/personal

## Global Directories 
PACKAGE="/rds/general/user/ph323/home/MRes.project.2/codes/0_packages"
DOCS="/rds/general/user/ph323/home/MRes.project.2/docs/0_preprocessing"
DATA_DIR="/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/raw"
tappAS="/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/tappAS"

## Global Controls 
SQANTI_QC=TRUE
SQANTI_FILTER=TRUE

# packages 
SQANTI="$PACKAGE/SQANTI3-5.2.1"

# functions 
quality_control="$SQANTI/sqanti3_qc.py"
filter="$SQANTI/sqanti3_filter.py"

# reference files 
primer="/rds/general/user/ph323/home/MRes.project.2/codes/1_preprocessing/primer.fasta"
ref_fa="$DOCS/hg38.fa"
ref_gtf="$DOCS/gencode.v39.chr_patch_hapl_scaff.annotation.gtf"
sorted_gtf="$DOCS/gencode.v39.chr_patch_hapl_scaff.annotation.sorted.gtf"
sample_list="/rds/general/user/ph323/home/MRes.project.2/codes/1_preprocessing/sample_list.txt"
polyAmotifs="$DOCS/mouse_and_human.polyA_motif.txt"
polyApeaks="$DOCS/atlas.clusters.2.0.GRCh38.96.bed"
cage="$DOCS/human.refTSS_v3.1.hg38.bed"
isoAnnotGFF="$DOCS/human_tappas_gencode_annotation_file.gff3"
intropolis="$DOCS/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified.gz"
rules="$DOCS/filter_custom.json"
short_read_fofn="$DOCS/pooled_short_reads.fofn"

## tappAS files 
collapsed_gff="$tappAS/isoseq3/collapsed.gff"
abundance="$tappAS/isoseq3/collapsed.abundance.txt"
classification="$tappAS/sqanti3/qc_classification.txt"
qc_corrected_gtf="$tappAS/sqanti3/qc_corrected.gtf"
qc_corrected_fasta="$tappAS/sqanti3/qc_corrected.fasta"
qc_corrected_faa="$tappAS/sqanti3/qc_corrected.faa"

## Step 5 - SQANTI3 Quality Control --------------------------------------------
if [ "$SQANTI_QC" = TRUE ]; then
    echo "Performing SQANTI3 quality control on the merged sample collapse.gff..."
    # activate environment 
    cd "$SQANTI"
    source ~/.bashrc
    conda activate SQANTI3.env
    # change to tappAS directory 
    cd "$tappAS/sqanti3"
    # perform 
    python $quality_control "$collapsed_gff" "$sorted_gtf" "$ref_fa" \
      --dir "$tappAS/sqanti3" \
      --fl_count "$abundance" \
      --coverage "$tappAS/sqanti3/STAR_mapping/" \
      --SR_bam "$tappAS/sqanti3/STAR_mapping/" \
      --polyA_motif_list "$polyAmotifs" \
      --polyA_peak "$polyApeaks" \
      --CAGE_peak "$cage" \
      --isoAnnotLite \
      --gff3 "$isoAnnotGFF" \
      --aligner_choice "minimap2" \
      --cpus 250 \
      --report both \
      --output "qc"
else
    echo "SQANTI3 quality control on the merged sample collapse.gff done, proceed to filtering"
fi

## Step 6 SQANTI3 filtering ----------------------------------------------------
if [ "$SQANTI_FILTER" = TRUE ]; then
    echo "Performing SQANTI3 quality control on the merged sample collapse.gff..."
    # activate environment 
    cd "$SQANTI"
    source ~/.bashrc
    conda activate SQANTI3.env
    # change to tappAS directory 
    cd "$tappAS/sqanti3"
    # perform 
    python $filter rules "$classification" \
      --output "filter_rule" \
      --dir "$tappAS/sqanti3" \
      --json_filter "$rules" \
      --gtf "$qc_corrected_gtf" \
      --faa "$qc_corrected_faa" \
      --isoforms "$qc_corrected_fasta"
    
else
    echo "SQANTI3 quality control on the merged sample collapse.gff done, proceed to filtering"
fi





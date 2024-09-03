#!/bin/bash

#PBS -l select=1:ncpus=256:mem=100gb
#PBS -l walltime=8:00:00          
#PBS -N isoseq3

## Global Directories 
PACKAGE="/rds/general/user/ph323/home/MRes.project.2/codes/0_packages"
DOCS="/rds/general/user/ph323/home/MRes.project.2/docs/0_preprocessing"
DATA_DIR="/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/raw"
tappAS="/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/tappAS"

## Global Controls 
ISOSEQ_LIMA=FALSE
UPDATE_BAM_HEADER=FALSE
MERGE_SMRT_CELLS=FALSE
ISOSEQ_CLUSTER_CLASSIFICATION=FALSE
SQANTI_QC=FALSE
SQANTI_FILTER=TRUE

# packages 
SMRT_LINK="$PACKAGE/smrtlink_v13/smrtcmds/bin"
scISA="$PACKAGE/scISA-Tools/bin"
BGI="$PACKAGE/bgi_commands/bgi_pipeline/bin"
ANACONDA="/rds/general/user/ph323/home/anaconda3/bin"
SQANTI="$PACKAGE/SQANTI3-5.2.1"

# functions 
ccs="$SMRT_LINK/ccs"
isoseq="$SMRT_LINK/isoseq"
lima="$SMRT_LINK/lima"
pbmerge="$SMRT_LINK/pbmerge"
pbmm2="$SMRT_LINK/pbmm2"
quality_control="$SQANTI/sqanti3_qc.py"
filter="$SQANTI/sqanti3_filter.py"
rescue="$SQANTI/sqanti3_rescue.py"

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
output_fofn="$tappAS/isoseq3/all_samples.fofn"
clustered_bam="$tappAS/isoseq3/clustered.bam"
mapped_bam="$tappAS/isoseq3/mapped.bam"
collapsed_gff="$tappAS/isoseq3/collapsed.gff"
abundance="$tappAS/isoseq3/collapsed.abundance.txt"
classification="$tappAS/sqanti3/qc_classification.txt"
qc_corrected_gtf="$tappAS/sqanti3/qc_corrected.gtf"
qc_corrected_fasta="$tappAS/sqanti3/qc_corrected.fasta"
qc_corrected_faa="$tappAS/sqanti3/qc_corrected.faa"


## Perform IsoSeq using the lima pipeline to obtain read group informations ---------------------
if [ "$ISOSEQ_LIMA" = TRUE ]; then
    echo "performing isoseq lima..."
    while read -r sample_id; do
      ## Directories
        start_time=$(date)
        DATA="$DATA_DIR/$sample_id"
        ISOSEQ3="$DATA/isoseq3"
        REPORT_FILES="$DATA/report_files"
      ## Files
        ccs_bam="$ISOSEQ3/ccs.bam"
        fl_bam="$ISOSEQ3/new_fl.bam"
        fl_bam_5to3="$ISOSEQ3/new_fl.primer_5p--primer_3p.bam"
        flnc_bam="$ISOSEQ3/new_flnc.bam"
        ccs_report="$REPORT_FILES/ccs_report.csv"
        ccs_log="$REPORT_FILES/ccs_log.txt"
        isoseq2_cluster2_log="$REPORT_FILES/isoseq2_cluster2_log.txt"
      ## Classify CCS by primer blast
        cd $ISOSEQ3
        echo "Currently at: $ISOSEQ3, now classifing ccs by primer blast!"
        $lima "$ccs_bam" $primer $fl_bam --isoseq  
        $isoseq refine "$fl_bam_5to3" $primer $flnc_bam
    done < "$sample_list"
else
    echo "isoseq lima is performed, code block skipped"
fi

## Step 3a - update bam header ===================================================
if [ "$UPDATE_BAM_HEADER" = TRUE ]; then
    echo "updating bam headers..."
    while read -r sample_id; do
        ## point toward directory 
          DATA="$DATA_DIR/$sample_id"
        ## files 
          isoseq_bam="$DATA/isoseq3/isoseq_flnc.bam"
          new_bam="$DATA/isoseq3/new_flnc.bam"
          isoseq_header="$DATA/isoseq3/isoseq_flnc_header.sam"
          new_header="$DATA/isoseq3/new_flnc_header.sam"
          updated_bam="$DATA/isoseq3/updated_isoseq_flnc.bam"
        ## updating...
          samtools view -H "$isoseq_bam" > "$isoseq_header" # create header 
          samtools view -H "$new_bam" > "$new_header"
          grep -E '^@HD|^@RG' "$new_header" > "${new_header}_hdrg" # Extract @HD and @RG
          awk 'NR==FNR { if (/^@HD/ || /^@RG/) {print; next} } 1' "${new_header}_hdrg" "$isoseq_header" > "${isoseq_header}_new" # Combine headers
          samtools reheader "${isoseq_header}_new" "$isoseq_bam" > "$updated_bam" # updating 
          samtools view -H "$updated_bam"
    done < "$sample_list"
else
    echo "bam headers are updated, code block skipped"
fi

## Step 3b - Merge SMRT Cells --------------------------------------------------
if [ "$MERGE_SMRT_CELLS" = TRUE ]; then
    echo "creating a fofn file recording all updated_isoseq_flnc.bam..."
    > "$output_fofn" # empty the file 
    while read -r sample_id; do
      flnc_bam="$DATA_DIR/$sample_id/isoseq3/updated_isoseq_flnc.bam"
      echo "$flnc_bam" >> "$output_fofn"
    done < "/rds/general/user/ph323/home/MRes.project.2/codes/0_sample_id/sample_list.txt"
else
    echo "fofn file recording all updated_isoseq_flnc.bam created, proceed to clustering"
fi

## Step 4 - Clustering and classification --------------------------------------
if [ "$ISOSEQ_CLUSTER_CLASSIFICATION" = TRUE ]; then
    echo "Clustering and classifying merged all sample bam file..."
    $isoseq cluster2 "$output_fofn" "$clustered_bam"
    $pbmm2 align --preset ISOSEQ --sort "$clustered_bam" "$ref_fa" "$mapped_bam"
    $isoseq collapse --do-not-collapse-extra-5exons "$mapped_bam" "$collapsed_gff"
else
    echo "Clustering and classifying of merged all sample bam file performed, proceed to SQANTI3"
fi

















  
#!/bin/bash
#PBS -l select=1:ncpus=8:mem=100gb:ompthreads=8
#PBS -l walltime=9:00:00          
#PBS -N X603

# Load Modules: 
module load picard/2.6.0
module load samtools/1.2
module load anaconda3/personal
source activate Renv

# Global Control: 
ADD_READ_GROUP=TRUE
MUTECT2=TRUE
GET_PILEUP_SUMMARY=TRUE
CALCULATE_CONTAMINATION=TRUE
FILTER_MUTECT2=TRUE
FUNCTIONAL_ANNOTATION=TRUE

# Directories: 
HOME="/rds/general/user/ph323"
REFERENCE="$HOME/home/MRes.project.2/docs/0_preprocessing/hg38_fasta/"
WES="$HOME/ephemeral/MRes.project.2/raw_data/WES"
TMP="$HOME/ephemeral/temp"

# Batch? 
paired_samples="$HOME/home/MRes.project.2/codes/0_sample_id/batch/temp.txt"

# Reference Files: 
reference_fasta="$REFERENCE/hg38_customized_noChr/hg38_customized.fasta"
reference_dict="$REFERENCE/hg38_customized_noChr/hg38_customized.dict"
reference_fasta_chr="$REFERENCE/hg38_customized_Chr/hg38_customized.fasta"
reference_dict_chr="$REFERENCE/hg38_customized_Chr/hg38_customized.dict"
interval_list_chr="$REFERENCE/interval_list/target_chr.interval_list"
interval_list="$REFERENCE/interval_list/target.interval_list"
germline="$REFERENCE/gnomad/af-only-gnomad.hg38.updated.vcf.gz"
PoN="$REFERENCE/panel_of_normal/somatic-hg38_1000g_pon.hg38.updated.vcf"
small_exac="$REFERENCE/small_exac/somatic-hg38_small_exac_common_3.hg38.updated.vcf.gz"
data_source="$REFERENCE/funcotator_dataSources_v1.7_hg38_somatic" 
conversion_table_rev="$REFERENCE/chr_conversion/chr_conv_map2.txt"

sample_id="X603"

# Run Mutect2 from GATK 
#while read -r sample_id; do

  mkdir -p "$WES/output/$sample_id"
  cd "$WES/output/$sample_id"   # go into output sub-directories 
  
  # if [ "$ADD_READ_GROUP" = TRUE ]; then 
  #   tumour_bam="$WES/Alignment_ov/${sample_id}T_sorted_nodup.bam"
  #   tumour_bam_rg="$WES/Alignment_ov/${sample_id}T_sorted_nodup_rg.bam"
  #   picard AddOrReplaceReadGroups I="$tumour_bam" O="$tumour_bam_rg" \
  #     RGID="$sample_id" RGLB=S07604514 RGPL=illumina RGPU=unit1 RGSM="${sample_id}T"
  #   samtools index "$tumour_bam_rg"
  # fi 
  
  # the raw bam files 
  tumour_bam="$WES/Alignment_ov/${sample_id}T_sorted_nodup_rg.bam"
  normal_bam="$WES/Alignment_ov/${sample_id}N_sorted_nodup_rg.bam"
  
  # the outputs 
  bamout="${sample_id}.bamout.bam"
  mutect2_vcf="${sample_id}.somatic.vcf.gz"
  mutect2_filtered_vcf="${sample_id}.filtered.somatic.vcf.gz"
  mutect2_filtered_chr_vcf="${sample_id}.filtered.somatic.chrRenamed.vcf.gz"
  mutect2_funcAnnot_vcf="${sample_id}.funcAnnot.somatic.vcf.gz"
  table_T="${sample_id}T_pileupsummaries.table"
  table_N="${sample_id}N_pileupsummaries.table"
  contamination="${sample_id}_contaminations.table"
  stats="${sample_id}.somatic.vcf.gz.stats"
  
  #samtools view -H "$tumour_bam" | grep "@SQ" | cut -f2 | sed 's/SN://g'
  
  if [ "$MUTECT2" = TRUE ]; then 
    echo "Performing Mutect2, estimate 80 minutes needed......"
    gatk Mutect2 \
      -R "$reference_fasta" \
      -I "$tumour_bam" \
      -I "$normal_bam" \
      -O "$mutect2_vcf" \
      --tmp-dir "$TMP" \
      --germline-resource "$germline" \
      --panel-of-normals "$PoN" \
      --normal-sample "${sample_id}N" \
      --bam-output "$bamout" \
      --native-pair-hmm-threads 8
      # --intervals "$interval_list" \
  fi 
  
  if [ "$GET_PILEUP_SUMMARY" = TRUE ]; then 
    echo "Getting pile up summaries, estimate 20 minutes needed......"
    gatk GetPileupSummaries \
      --input "$tumour_bam" \
      --intervals "$small_exac" \
      --variant "$small_exac" \
      --tmp-dir "$TMP" \
      --output "$table_T"
    gatk GetPileupSummaries \
      --input "$normal_bam" \
      --intervals "$small_exac" \
      --variant "$small_exac" \
      --tmp-dir "$TMP" \
      --output "$table_N"
  fi
  
  if [ "$CALCULATE_CONTAMINATION" = TRUE ]; then 
    echo "calculating contamination, estimated "
    gatk CalculateContamination \
      --input "$table_T" \
      --tmp-dir "$TMP" \
      --output "$contamination" \
      --matched-normal "$table_N"
    # gatk CalculateContamination \
    #   --input "$table_T" \
    #   --tmp-dir "$TMP" \
    #   --output "$contamination"
  fi
    
    
  if [ "$FILTER_MUTECT2" = TRUE ]; then 
    echo "filtering mutect2 vcf based on contamination"
    gatk FilterMutectCalls \
      --reference "$reference_fasta" \
      --variant "$mutect2_vcf" \
      --output "$mutect2_filtered_vcf" \
      --tmp-dir "$TMP" \
      --contamination-table "$contamination" \
      --stats "$stats" 
  fi

  # bcftools annotate --rename-chrs "$conversion_table_rev" "$mutect2_filtered_vcf" -Ov -o "$mutect2_filtered_chr_vcf"
  # gatk IndexFeatureFile -I "$mutect2_filtered_chr_vcf"
  # 
  # if [ "$FUNCTIONAL_ANNOTATION" = TRUE ]; then
  #   echo "performing functional annotation on filtered vcf file..."
  #   gatk Funcotator \
  #     --data-sources-path "$data_source" \
  #     --output "$mutect2_funcAnnot_vcf" \
  #     --output-file-format VCF \
  #     --reference "$reference_fasta_chr" \
  #     --ref-version hg38 \
  #     --tmp-dir "$TMP" \
  #     --intervals "$interval_list_chr" \
  #     --variant "$mutect2_filtered_chr_vcf"
  # fi

#done < "$paired_samples"














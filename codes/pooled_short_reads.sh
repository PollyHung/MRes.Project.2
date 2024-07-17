#!/bin/bash

output_fofn="/rds/general/user/ph323/home/MRes.project.2/docs/0_preprocessing/pooled_short_reads.fofn"
> "$output_fofn"
sample_list="/rds/general/user/ph323/home/MRes.project.2/docs/0_preprocessing/matched_samples.txt"
mapping_file="/rds/general/user/ph323/home/MRes.project.2/docs/0_preprocessing/mapping_table.csv"

# Loop through each sample_id in the matched_samples file
while read -r sample_id; do
  short_read_id=$(awk -F, -v sample="\"$sample_id\"" 'NR > 1 && $2 == sample {print $4}' "$mapping_file" | tr -d '"')
  SHORT_READS="/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/short_reads/RNA-seq/G72_run990/$short_read_id"

  short_read_R1=$(find "$SHORT_READS" -type f -name "*R1_001.fastq")
  short_read_R2=$(find "$SHORT_READS" -type f -name "*R2_001.fastq")
  echo $short_read_R1
  echo $short_read_R2

  echo "$short_read_R1 $short_read_R2" >> "$output_fofn"
done < $sample_list



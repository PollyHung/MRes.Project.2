#!/bin/bash

#PBS -l select=1:ncpus=1:mem=100gb
#PBS -l walltime=3:00:00          
#PBS -N bam2fastq


## convert ccs.bam to ccs.fastq for FLAIR 
module load bedtools/2.25
sample_list="/rds/general/user/ph323/home/MRes.project.2/codes/0_sample_id/sample_list.txt"
while read -r sample_id; do 
  cd "/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/raw/$sample_id/isoseq3"
  bedtools bamtofastq -i ccs.bam -fq ccs.fq
done < "$sample_list"

## check all short reads in sub-folders in directory
## if they are not in fastq and not fastq.gz format, gzip it. 
directory="/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/short_reads/RNA-seq/G72_run990"
short_read="/rds/general/user/ph323/home/MRes.project.2/codes/0_sample_id/short_reads_id.txt"

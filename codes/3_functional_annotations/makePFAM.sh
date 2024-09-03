#!/bin/bash
#PBS -l select=1:ncpus=10:mem=100gb
#PBS -l walltime=8:00:00          
#PBS -N makePFAM


cd /rds/general/user/ph323/home/MRes.project.2/docs/0_preprocessing/PFAM/

## download 
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz

## hmmpress
hmmpress Pfam-A.hmm

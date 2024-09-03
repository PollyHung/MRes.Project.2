#!/bin/bash

#PBS -l select=1:ncpus=127:mem=200gb
#PBS -l walltime=24:00:00          
#PBS -N transDecoder

## establishing environment 
cd /rds/general/user/ph323/home
source ~/.bashrc
conda activate base

## meta directories 
HOME="/rds/general/user/ph323/home/"
PACKAGE="$HOME/MRes.project.2/codes/0_packages"
ANACONDA="$HOME/anaconda3/"
DOCS="$HOME/MRes.project.2/docs/0_preprocessing"

## packages and functions 
module load blast+/2.11.1
module load hmmer/3.1 
LongOrfs="$ANACONDA/bin/TransDecoder.LongOrfs"
predict="$ANACONDA/bin/TransDecoder.Predict"
tmhmm="/rds/general/user/ph323/home/MRes.project.2/codes/0_packages/tmhmm-2.0c/bin/tmhmm"
makeDB="/rds/general/user/ph323/home/MRes.project.2/codes/1_preprocessing/makeDB.sh"
makePFAM="/rds/general/user/ph323/home/MRes.project.2/codes/1_preprocessing/makePFAM.sh"

start_time=$(date)
  
# directory ------------------------------------------------------------------
DATA_DIR="/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/tappAS"
SQANTI="$DATA_DIR/sqanti"
TRANSCODER="$DATA_DIR/transcoder"

cd $TRANSCODER
  
# files ----------------------------------------------------------------------
SQANTI_FA="$SQANTI/qc_corrected.fasta"
BLAST_DB="$DOCS/uniprot/human.fasta"
PFAM_DB="$DOCS/PFAM/Pfam-A.hmm"
PEP="$TRANSCODER/qc_corrected.fasta.transdecoder_dir/longest_orfs.pep"
RETAIN_PEP="$TRANSCODER/qc_corrected.fasta.transdecoder.pep"
pfam_domtblout="$DOCS/PFAM/pfam.domtblout"
blast_out="$TRANSCODER/blastp.out"
  
# do TransDecoder! -----------------------------------------------------------
#$LongOrfs -t $SQANTI_FA
  
#blastp -query $PEP  \
#       -db $BLAST_DB \
#       -max_target_seqs 1 \
#       -outfmt 6 \
#       -evalue 1e-5 \
#       -num_threads 256 > blastp.out

hmmscan --cpu 122 --domtblout pfam.domtblout $PFAM_DB $PEP
hmmscan_time=$(date)
echo "hmmscan done at $hmmscan_time"

$predict \
    -t $SQANTI_FA \
    --retain_pfam_hits pfam.domtblout \
    --retain_blastp_hits blastp.out \
    --single_best_only
predict_time=$(date)  
echo "transdecoder predict done at $predict_time"
  
cat $RETAIN_PEP | $tmhmm -workdir test > out_tmhmm.txt
grep -v "#" out_tmhmm.txt > tmhmm_regions.txt
tmhmm_time=$(date)
echo "transdecoder predict done at $tmhmm_time"
  
# time stamp -----------------------------------------------------------------
end_time=$(date)
echo "TransDecoder for $sample_id from $start_time to $end_time"


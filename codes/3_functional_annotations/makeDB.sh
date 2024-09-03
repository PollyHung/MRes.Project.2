#!/bin/bash
#PBS -l select=1:ncpus=10:mem=100gb
#PBS -l walltime=8:00:00          
#PBS -N makeDB

module load blast+/2.11.1
seqkit="/rds/general/user/ph323/home/anaconda3/bin/seqkit"

cd /rds/general/user/ph323/home/MRes.project.2/docs/0_preprocessing/uniprot

## download 
wget -c "ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.fasta.gz"
wget -c "ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz"
wget -c "ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz"

## unzip it 
gunzip uniprot_sprot.fasta.gz
gunzip uniprot_trembl.fasta.gz
gunzip uniprot_sprot_varsplic.fasta.gz

## combine the SwissProt, TrEMBL, and VarSplice datasets into one file
cat uniprot_sprot.fasta uniprot_trembl.fasta uniprot_sprot_varsplic.fasta > allProteins_human.fasta

# Filter for unique sequences
$seqkit rmdup -s -i allProteins_human.fasta -o allProteins_unique_human.fasta
## [INFO] 23008755 duplicated records removed

## count how many protein sequences is in this combined fasta 
grep -c "^>" allProteins_human.fasta
# there are 245523696 sequences in total
grep -c "^>" allProteins_unique_human.fasta
# after removing duplicates, there are still 222514941 sequences 

## Okay turns out that mostly are not human, we'll filter for human
grep -c "OS=Homo sapiens" allProteins_unique_human.fasta
# after filtering for HUMAN only sequences, there are 196955 in total, which is 
# still 100000 more than the 95915 mentioned in paper but much better than 222514941

## let's filter 
input_file="allProteins_unique_human.fasta"
output_file="human.fasta"

# Filter for sequences where OS=Homo sapiens
awk '
    BEGIN {print_seq = 0} 
    /^>/ {print_seq = ($0 ~ /OS=Homo sapiens/)}
    {if (print_seq) print}
' $input_file > $output_file

# count again...
grep -c "^>" human.fasta

# now remove the old files 
rm allProteins*
rm uniprot*

## makeblastdb
makeblastdb -in human.fasta -dbtype prot -parse_seqids


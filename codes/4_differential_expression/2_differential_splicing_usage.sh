## This script performs differential splicing usage analysis 
## Aim: identify changes in the relative usage of different transcripts (isoforms) 
##      of the same gene between normal and tumour
## Methods: 
##        SUPPA –– 


## install SUPPA ---------------------------------------------------------------
pip3 install SUPPA==2.3 
pip3 check SUPPA
pip3 show SUPPA
#Name: SUPPA
#Version: 2.3
#Summary: A tool to study splicing across multiple conditions at high speed and accuracy.
#Home-page: https://github.com/comprna/SUPPA
#Author: GP Alamancos
#Author-email: eduardo.eyras@upf.edu
#License: MIT
#Location: /Library/Frameworks/Python.framework/Versions/3.9/lib/python3.9/site-packages
#Requires: numpy, pandas, scikit-learn, scipy, statsmodels
#Required-by: 

## establish connections 
tappAS="/Users/pollyhung/Desktop/tappAS"
suppa="$tappAS/packages/SUPPA/suppa.py"
split_files="$tappAS/packages/SUPPA/scripts/split_file.R"
split_psi="$tappAS/packages/SUPPA/scripts/split_psi.R"
sqanti_gtf="$tappAS/docs/sqanti/filter_rule.filtered.gtf"
expression_matrix="$tappAS/docs/processed/tpm_normalized.txt"

## generate alternative splicing events ----------------------------------------
cd "$tappAS/docs/suppa"
# run the command 
python3 "$suppa" generateEvents -i "$sqanti_gtf" -o "sqanti_filter.events" -e SE SS MX RI FL -f ioe
# put all the ioe events in the same file 
awk '
    FNR==1 && NR!=1 { while (/^<header>/) getline; }
    1 {print}
' *.ioe > sqanti_filter.events.ioe
# put all gtf files in the same file
cat *.gtf > sqanti_filter.events.gtf


## obtain PSI values for splicing events ---------------------------------------
# calculate PSI for local events 
python3 "$suppa" psiPerEvent -i sqanti_filter.events.ioe -e "$expression_matrix" -o "psi_local" -m INFO
# calculate PSI per transcript Isoform 
python3 "$suppa" psiPerIsoform -g "$sqanti_gtf" -e "$expression_matrix" -o "psi" -m INFO


## differential splicing analysis ----------------------------------------------
# split the PSI and TPM files between two conditions 
Rscript "$split_files" "$expression_matrix" "HH12000124FN2,HH12000157FN4,T14_042FN3,T15_036_FN2,T16_002FN1,T16_035_FN2" "0700055A,0700181A,0800062A,090061A,10_149A,1100074A,HH15000080_FT1,T14_049_FT2,T15_022_FT2,T15_051_FT4,T15_162_FT2,T15_171_FT4,T16_081_FT2,T16_088_FT2,T16_106_FT1,T17_046_FT3,T17_205_FT3,T18_119_FT13,T18_158_FT1" normal.tpm cancer.tpm -i 
Rscript "$split_psi" psi_local.psi "HH12000124FN2,HH12000157FN4,T14_042FN3,T15_036_FN2,T16_002FN1,T16_035_FN2" "0700055A,0700181A,0800062A,090061A,10_149A,1100074A,HH15000080_FT1,T14_049_FT2,T15_022_FT2,T15_051_FT4,T15_162_FT2,T15_171_FT4,T16_081_FT2,T16_088_FT2,T16_106_FT1,T17_046_FT3,T17_205_FT3,T18_119_FT13,T18_158_FT1" normal.psi cancer.psi -e
# perform differential splicing analysis 
python3 "$suppa" diffSplice -m empirical -i sqanti_filter.events.ioe -p cancer.psi normal.psi -e cancer.tpm normal.tpm -o "diffSplice" --save_tpm_events




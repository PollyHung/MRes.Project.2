## packages 
HOME="/Volumes/Wild_Flower/OV_SpliceVariants"
ANNOVAR="$HOME/packages/annovar"
MUTECT2="$HOME/data/4_mutect2/output"

## functions 
CONVERT2ANNOVAR="$ANNOVAR/convert2annovar.pl"
TABLE_ANNOVAR="$ANNOVAR/table_annovar.pl"
ANNOTATE_VARIATION="$ANNOVAR/annotate_variation.pl"
RETRIEVE_SEQ_FROM_FASTQ="$ANNOVAR/retrieve_seq_from_fasta.pl"

## files 
sample_id="$HOME/code/0_sample_id/WES/WES_paired_patient_list.txt"
humandb="$ANNOVAR/humandb"

## Download databases 
#$ANNOTATE_VARIATION -downdb -buildver hg38 -webfrom annovar cytoBand humandb/
#$ANNOTATE_VARIATION -downdb -buildver hg38 -webfrom annovar ensGene humandb/
#$ANNOTATE_VARIATION -downdb -buildver hg38 -webfrom annovar refGene humandb/
#$ANNOTATE_VARIATION -downdb -buildver hg38 -webfrom annovar dbnsfp30a humandb/
#$ANNOTATE_VARIATION -downdb -buildver hg38 -webfrom annovar exac03 humandb/
#$ANNOTATE_VARIATION -downdb -buildver hg38 -webfrom annovar avsnp147 humandb/
#$ANNOTATE_VARIATION -downdb -buildver hg38 -downdb cytoBand humandb/

## turn vcf into a txt 
#while read -r wes_id; do 
  wes_id="X603"
  cd "$MUTECT2/$wes_id"
  
  vcf_gz="$MUTECT2/$wes_id/$wes_id.filtered.somatic.vcf.gz"
  vcf="$MUTECT2/$wes_id/$wes_id.filtered.somatic.vcf"
  avinput="$MUTECT2/$wes_id/annovar_filt.${wes_id}T.avinput"
  
  gunzip $vcf_gz
  
  $CONVERT2ANNOVAR \
    --format vcf4 \
    --coverage 25 \
    --fraction 0.05 \
    --allsample \
    --outfile "annovar_filt" \
    $vcf
  
  $TABLE_ANNOVAR "$avinput" "$humandb" \
    --protocol ensGene,cytoBand,exac03,avsnp147,dbnsfp30a \
    --operation gx,r,f,f,f  \
    --outfile "annovar_annot" \
    --buildver hg38 \
    --remove \
    --polish \
    --nastring . 
  
#done < "$sample_id"



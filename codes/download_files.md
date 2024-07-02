# File Download Details    
This is a script showing how and were packages and reference files are acquired, most actions took place in command line, but if otherwise it will be stated. 

## Download Packages 
### SMRTLINK v13 
Website: https://www.pacb.com/support/software-downloads/   
```
wget https://downloads.pacbcloud.com/public/software/installers/smrtlink-release-sequel2_13.1.0.221970.zip
unzip smrtlink-release-sequel2_13.1.0.221970.zip
chmod +x smrtlink-release-sequel2_13.1.0.221970_linux_x86-64_libc-2.17_anydistro.run
./smrtlink-release-sequel2_13.1.0.221970_linux_x86-64_libc-2.17_anydistro.run --rootdir smrtlink --smrttools-only
```

### BGI functions
Website: https://github.com/shizhuoxing/BGI-Full-Length-RNA-Analysis-Pipeline/tree/master for bulk sequencing    
Website: https://github.com/shizhuoxing/scISA-Tools/tree/master for single cell sequencing    
```
git clone https://github.com/shizhuoxing/BGI-Full-Length-RNA-Analysis-Pipeline.git
mv BGI-Full-Length-RNA-Analysis-Pipeline bgi_commands
ls bgi_commands/bin
## here is where the commands are located
## classify_by_primer.fullpa.pl  classify_by_primer.pl  flnc2sam.pl  PolymeraseReads.stat.pl  SubReads.stat.pl
```

### SQUANTI3
Website: https://github.com/ConesaLab/SQANTI3    
```
wget https://github.com/ConesaLab/SQANTI3/archive/refs/tags/v5.2.1.tar.gz
tar -xf v5.2.1.tar.gz
cd SQANTI3-5.2.1
conda env create -f SQANTI3.conda_env.yml
conda activate SQANTI3.env
## once activated the virtual environemnt, you should see prompt changing to:
(SQANTI3.env)$
```

### TransDecoder 
Website: https://github.com/TransDecoder/TransDecoder/wiki
```
## install from source package
wget https://github.com/TransDecoder/TransDecoder/archive/refs/tags/TransDecoder-v5.7.1.tar.gz
tar -xzvf TransDecoder-v5.7.1.tar.gz
## or use conda install
conda install bioconda::transdecoder
```


## Download Reference Files 
### primer.fasta
Website: https://github.com/shizhuoxing/BGI-Full-Length-RNA-Analysis-Pipeline/tree/master
```
nano primer.fa
```
input: 
```
>primer_F
AAGCAGTGGTATCAACGCAGAGTACATGGGGGGGG
>primer_S
GTACTCTGCGTTGATACCACTGCTTACTAGT
```
save file 

### hg38.fasta
Website: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa
```

### annotation.gtf and annotation.gff
Website: https://www.gencodegenes.org/human/
```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz
gunzip gencode.v46.annotation.gtf.gz
chmod +x gencode.v46.annotation.gtf

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gff3.gz
gunzip gencode.v46.annotation.gff3.gz
chmod +x gencode.v46.annotation.gff3

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz
gunzip gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz
chmod +x gencode.v46.chr_patch_hapl_scaff.annotation.gtf

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.chr_patch_hapl_scaff.annotation.gff3.gz
gunzip gencode.v46.chr_patch_hapl_scaff.annotation.gff3.gz
chmod +x gencode.v46.chr_patch_hapl_scaff.annotation.gff3
```

### polyA data    
There are two polyA files required for SQANTI quality control:    
1. PolyA motif data    
2. PolyA site data
```
wget https://github.com/ConesaLab/SQANTI3/blob/master/data/polyA_motifs/mouse_and_human.polyA_motif.txt
wget https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz
gunzip atlas.clusters.2.0.GRCh38.96.bed.gz
## for the second atlas.cluster.bed file we'll need to modify the missing chr prefix by:
sed -i 's/^/chr/' atlas.clusters.2.0.GRCh38.96.bed
```
and the polyA motif file look something like this: 
```
aataaa
attaaa
agtaaa
tataaa
cataaa
gataaa
aatata
aataca
aataga
aaaaag
actaaa
aagaaa
aatgaa
tttaaa
aaaaca
ggggct
```

### CAGE peak data  
Obtained public CAGE peak data for human and mice from the refTSS database supplied in SQANTI3's `data` folder.     
```
wget https://github.com/ConesaLab/SQANTI3/blob/master/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed
```

### Intropolis junction bed file     
Provided by SQANTI github page in previous version      
```
wget https://github.com/Magdoll/images_public/blob/master/SQANTI2_support_data/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified.gz
```

### uniprot_sprot.fasta
```
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
```
























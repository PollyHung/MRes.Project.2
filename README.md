## Pipeline     

### Section 1 - Preprocessing subreads.bam with BGI pipeline     
Raw PacBio Iso-Seq data in `subreads.bam` format were processed using a combination of the [BGI-Full-Length-RNA-Analysis-Pipeline](https://github.com/shizhuoxing/BGI-Full-Length-RNA-Analysis-Pipeline) and the official [PacBio IsoSeq3 pipeline](https://github.com/PacificBiosciences/IsoSeq). Initially, each sequencing run was processed with `ccs` to generate HiFi reads (predicted accuracy ≥ Q20) using the parameters: min-passes 0, min-length 50, max-length 21000, and min-rq 0.75. This ensured that the resulting HiFi reads had at least 1 full-length subreads, lengths between 50 and 21000 bp, and a minimum read accuracy of 0.75.     

Instead of using PacBio IsoSeq3 pipeline for subsequent processing, we jumpted to BGI pipeline from this step because of [Issue 1](/logs/Issue_1). Primers were removed using in-house scripts from the BGI pipeline, with forward primer (5') sequence `AAGCAGTGGTATCAACGCAGAGTACATGGGGGGGG` and the reverse primer (3') sequence `GTACTCTGCGTTGATACCACTGCTTACTAGT`. These primers were mapped to the CCS reads using NCBI BLAST (version 2.11.1). The BLAST results were utilized by the `classify_by_primer.py` utility to identify full-length transcripts based on the BGI patented multi-transcripts approach. FLNC reads were generated with parameters: minimum UMI length of 8, minimum aligned primer length is 16, and minimum output transcript length without polyA tail is 200. The `classify_by_primer.py` script then performed the following tasks: (1) parsing 5’ and 3’ primers in CCS reads to obtain FLNC reads oriented from 5’ to 3’; (2) trimming the 5’ and 3’ primer sequences, including 28 bp following the 3’ primer; and (3) trimming the 3’ polyA tail using a sliding window algorithm. Subsequently, `flnc.sam` was created from `ccs.sam` and `flnc.fasta` using the `flnc2sam` function and converted to `flnc.bam` using `samtools`. This completed the steps in the BGI pipeline, after which processing continued with the Iso-Seq3 pipeline.   

### Section 2 - Preprocessing isoseq_flnc.bam with IsoSeq3 pipeline      
#### Approach 1 – processing long read sequences separately      
Script: [preprocess.sh](/codes/preprocess.sh)      

For each individual long read sequence, `isoseq3 cluster2` function performs a clustering step beginning by clustering using hierarchical n•log(n) alignemnt and iterative cluater merging, followed by polishing POA sequence generation using a QV-guided consensus approach. The resulting output contains isoforms with at least 2 FLNC reads, all should have a predicted accuracy of ≥0.99, as HiFi reads were used as input.     

#### Approach 2 – processing long read sequences together     
Script: [preprocess_pooled.sh](/codes/preprocess_pooled.sh)    

Due to [Issue 2](/logs/Issue_2), we decide to perform additional step to generate one transcriptome for all sequenced samples because we want to generate a single expression matrix for subsequent differential exon expression analysis. This is done by adding an additional step after `isoseq_flnc.bam` and before `isoseq3 cluster2` to merge the SMRT cells. Since we processed our `isoseq_flnc.bam` using customised BGI script, here we need to add read group information to these `isoseq_flnc.bam`. @RG and @HD are obtained by running the origianl PacBio IsoSeq3 `lima` on `ccs.bam` and `refine` on the `fl.bam` to produce a `new_flnc.bam` containing @RG, @HD, and @PG. Sample corresponded @RG and @HD are stripped from the `new_flnc.bam` using samtools and pasted to our BGI produced `isoseq_flnc.bam`. The paths of new `updated_isoseq_flnc.bam` from all samples are stored in a `output.fofn` file with each path separated by a line. This `output.fofn` file is then used in the `isoseq3 cluster2` function to call all long read samples and produce a single `clustered.bam` file.      

From `clustered.bam` onwards the pipeline remain the same, in short we employed Pigeon to classify full-length transcript isoforms against a reference annotation. Initially, the `clustered.ba`m was mapped to hg38 reference genome using `pbmm2`. Redundant transcripts were then collapsed into unique isoforms based on exonic structures using `isoseq3 collapse`, ensuring not to collapse extra 5' exons, as bulk sequencing was performed. The genome annotation [gencode.v39.chr_patch_hapl_scaff.annotation.gtf](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.chr_patch_hapl_scaff.annotation.gtf.gz), [hg38.fasta](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz), and the `collapsed.gff` from isoseq collapse were then sorted and indexed. The output `collapsed.gff` is then passed to SQANTI3 pipeline for quality control and isoform annotation.            

### Section 3 - Quality Control, Curation and Annotation of Long-read Transcripts using SQANTI3      
Scripts: [sqanti_matched.sh](/codes/sqanti_matched.sh), [pooled_short_reads.sh](/codes/pooled_short_reads.sh) and [sqanti_unmatched.sh](/codes/sqanti_unmatched.sh)     

SQANTI3 is the newest version integrating functionalities from both SQANTI and SQANTI2. The core module, the quality control, evaluates the de novo transcriptome against a reference annotation and outputs a classification file. The `sqanti3_qc.py` script takes in three compulsory files: (1) `collapsed.gff` from IsoSeq3, representing our sample transcript, (2) sorted and indexed `gencode.v39.chr_patch_hapl_scaff.annotation.sorted.gtf`, and (3) `hg38.fasta` file used in previous pipeline. Optional files were supplied when possible, public reference files including [polyA motifs](/docs/mouse_and_human.polyA_motifs.txt), [polyA peaks](/docs/atlas.clusters.2.0.GRCh38.96.bed.zip), [public CAGE peak data](/docs/human.refTSS_v3.1.hg38.bed.zip), [human TAPPAS Gencode v39 annotation](https://app.tappas.org/resources/downloads/gffs/Homo_sapiens_Gencode_v39.zip) and [intropolis junctions](https://github.com/Magdoll/images_public/blob/master/SQANTI2_support_data/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified.gz) were downloaded from SQANTI3 wiki page. Sample specific reference files including FLNC counts from IsoSeq3 pipeline and matched/pooled short reads fastq (R1/R2) were also provided. For processing long read sequences separately, 18/25 samples have [matched short reads](docs/short_read_to_long_read_map.csv) and the other 7/25 used a [pooled short read](pooled_short_read.fofn) as reference. The pooled short read were also used as reference for the compiled long read sequence. Alignment was carried out using `minimap2`.   

The output classification table consists of the following main categories of splicing events: full splice match (FSM), incomplete splice match (ISM), novel in category (NIC), novel not in catalog (NNC) and others (antisense, fusion, genic, genic intron, and intergenic). FSM is defined when reference and query isoform have the same exon count with same internal junction aligning to the reference. ISM is defined when query isoform have fewer outer exons compared to reference but still have all internal junctions aligned to reference. NIC is characterized by new combination of known, catalogued splice junctions (pairs of donor-acceptor sites) while NNC is very similar but have at least one splice site (donor or acceptor) entirely novel and uncatalogued before.    

We employed a [customized filtering rule](filter_custom.json) inspired by Veiga et al., (2022). Briefly, FSM: all retained; ISM: filtered out unreliable 3’end (if percentage of genomic A’s in the downstream 20bp window perc_A_dowmstreamTTS > 0.80); NIC: filtered out unreliable 3’end, intron-retention (in subcategory column), minimum coverage below 10; NNC: filtered out unreliable 3’end, intron-retention, junction labelled as RT-switch, minimum coverage below 10, and non-canonical splice sites; the rest of the splicing events are filtered out for intron-retention, junctions labelled as RT-switch, minimum coverage below 10, non-canonical splicing sites.      

### Section 4 - Functional Annotations of Isoforms     
Scripts: [makeDB.sh](/codes/makeDB.sh), [makePFAM.sh](/codes/makePFAM.sh), [transDecoder.sh](/codes/transDecoder.sh), and [deeploc2](/codes/deeploc2/sh)       

To understand the potential functional consequences of isoforms from our transcriptome at the protein level, we extracted open reading frames (ORFs; i.e., coding sequences) using TransDecoder and predicted domains using PFAM domain and BLAST, transmembrane regions using TMHMM2.0 and hmmer, and subcellular localization using DeepLoc2.      

Human proteome reference is assembled using both canonical ([SwissProt](ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz) + [TrEmbl](ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.fasta.gz)) and spliced isoforms ([VarSplice](ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz)) from Uniprot, filtering for `homo sapiens` as organism and removing duplicated sequences, resulting in 196955 sequences in total. Database is built using the `makeblastdb` function from blast+ (v2.11.1) package with `dbtype` set to `prot` and resulting file is `human.fasta`. Then, possible coding sequences (open reading frames) from isoforms were predicted using TransDecoder using `TransDecoder.LongOrfs` script providing the `qc_corrected.fasta` from SQANTI3 `quality_control.py`. The predicted amino acid sequences from TransDecoders are stored in `longest_orfs.pep` and is feeded as input for local alignment using blastp against the reference proteome `human.fasta` we built previously using SwissProt+TrEmbl+VarSplice with max_target_seqs = 1 and e-value = 10−5 to identify homologs in UniProt.     

We downloaded the compressed Pfam-A HMM database file ([Pfam-A.hmm.gz](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz)) from the Pfam FTP server and ran `hmmpress` from `HMMER` (v.3.1) on `Pfam-A.hmm` to prepare the HMM database for use by indexing it and creating binary files that can be read quickly by HMMER tools. PFAM domains for extracted ORFs were predicted using `hmmscan` with default parameters, and ouput of search is stored in `pfam.domtblout` containing a tabular form of domain annotation results for the query sequences. Then `TransDecoder.Predict` is use to select a single best ORF for each transcript based on domain conservation (pfam.domtblout) and significant sequence homology (blastp.out) to human proteins.       

Prediction and annotation of the subcellular localization and transmembrane helices of the predicted coding sequence (ORF) were carried our by TMHMM2.0 and DeepLoc2. Briefly, the TMHMM2.0 program takes in protein sequences in fasta format and recognizes 20 amino acids and B, Z, and X, all other unknown characters are treated as X. The output is a 5 column table with column 1 representing transcript ID, column 3 indicating where the peptide sequence is located (inside, outside or TMhelix), and column 4 and 5 recording the corresponding amino acid sequences of that segment of peptide. These gives 5 possible categories: peptide is outside, inside, single pass with part of the peptide located out/inside membrane or multi-pass peptide. The DeepLoc2 program takes in the same protein sequence fasta file and outputs a table with column 1 representing transcript ID, column 2 recording the two most possible subcellular location of the peptide, column 3 recording the signal detected and column 4 to 13 recording the possibility of this peptide being at cytoplasm, nucleus, extracellular, cell membrane, mitochondria, plastid, endoplasmic reticulum, lysosome/vacuole, golgi apparatus or peroxisome. Both progams were also applied to the protein reference `human.fasta` to create a reference transmembrane location and subcellular localization table.   

### Section 5 - Differential Transcript Analysis     

























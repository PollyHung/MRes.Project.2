## Merged process or separate process samples?      

#### problem:      
We initially processed our long reads samples separately, however this causes a problem that the pacbio IDs from each samples do not match. Because pacbio IDs were assigned randonly and floats around between different samples. Meaning that PB.100.1 in sample 1 does not label the same transcript or even the same gene for sample 2.      
Separate processing our samples therefore leads to a problem of us unable to merge the isoform expressions for differential exon expression analysis as the PacBio ID do not match and therefore novel transcripts would be lost as they do not have a consensus ensemble transcript ID we could use as a reference.     
To mitigate this, I first thought of creating a new mapping ID for different samples so we could merge everything together, but that didn't work out.      
Then, I found this section in tappAS FAQ        
```
I processed my long read samples separately, but tappAS requires one transcriptome for all samples. How do I generate one?
The best way to use tappAS is to generate one transcriptome for all sequenced samples. To do this, users should pre-process all SMRT cells together using IsoSeq3, and then run SQANTI3 on the output of this joint run. IsoSeq3 documentation contains more info on how to merge SMRT cells -typically all users will need to do is merge the output of the refine command, and then repeat the clustering step. See the IsoSeq3 documentation IsoSeq3 documentation for details.

WARNING: note that the QC report produced by SQANTI3 can help you make informed decisions about isoforms that might be false positives or low quality, and we strongly advise you to remove them from your transcriptome before you continue your analysis. We recommend reading the SQANTI paper to get a better idea of how to produce a high-quality, curated transcriptome.
```
So I went back to the IsoSeq3 step to redo the whole process. Except the only problem is that because we used a [customed script]([url](https://github.com/shizhuoxing/BGI-Full-Length-RNA-Analysis-Pipeline.git)) to do the preprocessing step to obtain isoseq_flnc.bam, therefore when processing each samples separately we could sneak by it but when processing all samples together using the .fofn format `pbmerge` cannot process the non-pacbio pipeline produced isoseq_flnc.bam files as they lack the necessary read group informations. To obtain the corresponding read group information for each bam files, we preprocessed all the samples using the traditional     
```
$lima "$ccs_bam" $primer $fl_bam --isoseq  
$isoseq refine "$fl_bam_5to3" $primer $flnc_bam
```
After that extract the @RG and @HD from the new_flnc.bam and added them to our custom-script-produced isoseq_flnc.bam file. After that we performed    
```
$isoseq cluster2 "$output_fofn" "$clustered_bam"
$pbmm2 align --preset ISOSEQ --sort "$clustered_bam" "$ref_fa" "$mapped_bam"
$isoseq collapse --do-not-collapse-extra-5exons "$mapped_bam" "$collapsed_gff"
```
and successfully obtained the files needed for next step processing.    


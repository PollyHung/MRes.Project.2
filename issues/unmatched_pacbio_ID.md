## Issue 1: Unmatched PacBio IDs if processing Long Read samples separately using IsoSeq3-SQANTI3-IsoAnnotLite pipeline    

#### Problem Description:   
Initially we processed our long read samples separately following the IsoSeq3-SQANTI3-IsoAnnotLite pipeline. However this causes a problem during the differential exon expression analysis step. Specifically, the PacBio IDs from each samples cannot be used universially and is specfic to each sample, that is to say PB.100.1 in sample 1 does not point toward the same gene as PB.100.1 in sample 2. Therefore although we have acquired the isoform expression at transcript level for each sample we cannot merge them into a transcription matrix due to the issue of unmatched PacBio ID.    




To mitigate this problem, we first thought that merging the qc_corrected.gtf, qc_classification.txt, and qc_junctions.txt (from SQANTI3 quality_control.py) from each individual long read samples. Specifically merging the qc_corrected.gtf first, then identify transcript boxes 
<img width="800" alt="Screenshot 2024-07-17 at 12 39 21" src="https://github.com/user-attachments/assets/e81400af-1c14-4920-9ce7-679f8ea5e7b4">     
with each box beginning with a line of transcript followed by several lines of exons. Then removing duplicated boxes and re-order the boxes by their transcript start site, and lastly assigning each box a new PacBio ID. However, soon we found out that PacBio ID assignment was not completely random as we though they would be, they actually reference to reference gene start and end sites meaning that we cannot assign each box a random ID but have to guess if different boxes could potentially be different transcripts of the same gene. And I'm stuck at this step.     

Our second plan is to go back to the IsoSeq3 pipeline and merge the flnc.bam files after the step of IsoSeq3 refine and perform IsoSeq3 and SQANTI3 pipeline on the merged files. This method is inspired by a [FAQ]([url](https://app.tappas.org/faqs/)) seen on tappAS website, specifically:    
<img width="800" alt="Screenshot 2024-07-17 at 12 44 39" src="https://github.com/user-attachments/assets/743d2c9d-c26c-499f-af7c-28b190ac53d4">   
This method requires us to provide a .fofn file recording the pathway toward each flnc.bam. However, because we didn't use the PacBio IsoSeq3 pipeline from beginning to end but rather acquired our flnc.bam through a [customed script]([url](https://github.com/shizhuoxing/BGI-Full-Length-RNA-Analysis-Pipeline.git)) provided by BGI. Therefore when we tried to perform     
```
$isoseq cluster2 "$output_fofn" "$clustered_bam"
```
we keep running into the issue of     
```
[pbbam] version string parsing ERROR: failed to parse:
  version: 
  reason: [pbbam] version string parsing ERROR: empty string
```
Looking into the pbbam source code we located the problem to be coming from this [Version.cpp](https://github.com/PacificBiosciences/pbbam/blob/develop/src/Version.cpp):     
```
...
// string must be "<major>.<minor>.<version>"
Version::Version(const std::string& v) : major_{0}, minor_{0}, revision_{0}
{
    // parse string
    try {
        const auto fields = Split(v, '.');
        const auto numFields = fields.size();
        if (numFields == 0) {
            throw std::runtime_error{"[pbbam] version string parsing ERROR: empty string"};
        }
        major_ = std::stoi(fields.at(0));
        if (numFields > 1) {
            minor_ = std::stoi(fields.at(1));
            if (numFields > 2) {
                revision_ = std::stoi(fields.at(2));
            }
        }
    } catch (std::exception& e) {
        std::ostringstream msg;
        msg << "[pbbam] version string parsing ERROR: failed to parse:\n"
            << "  version: " << v << '\n'
            << "  reason: " << e.what();
        throw std::runtime_error{msg.str()};
    }

    // ensure valid numbers
    Check();
}
...
```
Which essentially asked for a version string from the .bam file passed to it, probably a version string found in read group of the bam file. Therefore we suspected that maybe pacbio IsoSeq3 pipeline produces a read group attached to bam file that is not found in our customized script. As a result, we run a test sample using the IsoSeq3 pipeline and compared the headers between two bam files.    
```
samtools view -H isoseq_flnc.bam > isoseq_flnc_header.sam
samtools view -H new_flnc.bam > new_flnc_header.sam
diff isoseq_flnc_header.sam new_flnc_header.sam
```
Here is what we found:    
```
1,2c1,6
< @PG   ID:samtools     PN:samtools     VN:1.14 CL:samtools view -bS /rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/raw/0700055A/processed_data//isoseq_flnc.sam
< @PG   ID:samtools.1   PN:samtools     PP:samtools     VN:1.18 CL:samtools view -H isoseq_flnc.bam
---
> @HD   VN:1.6  SO:unknown      pb:5.0.0
> @RG   ID:b8eae6f2/0--1        PL:PACBIO       DS:READTYPE=CCS;BINDINGKIT=101-894-200;SEQUENCINGKIT=101-826-100;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100.000000;BarcodeFile=/rds/general/user/ph323/home/MRes.project.2/codes/1_preprocessing/primer.fasta;BarcodeHash=0a8631e04f02d91fc7b338fcd4067299;BarcodeCount=2;BarcodeMode=None;BarcodeQuality=Score     LB:WHPBCDNAPEP00000249  PU:m64048_240607_071033 SM:WHPBCDNAPEP00000249   PM:SEQUELII     BC:AAGCAGTGGTATCAACGCAGAGTACATGGGGGGGG-GTACTCTGCGTTGATACCACTGCTTACTAGT  CM:S/P5-C2/5.0-8M
> @PG   ID:ccs  PN:ccs  VN:8.0.1 (commit v8.0.1)        DS:Generate circular consensus sequences (ccs) from subreads.   CL:/rds/general/user/ph323/home/MRes.project.2/codes/0_packages/smrtlink/install/smrtlink-release_13.1.0.221970/bundles/smrttools/install/smrttools-release_13.1.0.221970/private/pacbio/unanimity/binwrap/../../../../private/pacbio/unanimity/bin/ccs /rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/raw/0700055A/8522212004059.bam /rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/raw/0700055A/processed_data//ccs.bam --min-passes 0 --min-length 50 --max-length 21000 --min-rq 0.75 --num-threads 60 --log-level INFO --reportFile /rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/raw/0700055A/report_files//ccs_report.csv --logFile /rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/raw/0700055A/report_files//ccs_log.txt
> @PG   ID:lima VN:2.10.0 (commit v2.10.0)      CL:/rds/general/user/ph323/home/MRes.project.2/codes/0_packages/smrtlink_v13/install/smrtlink-release_13.1.0.221970/bundles/smrttools/install/smrttools-release_13.1.0.221970/private/pacbio/barcoding/binwrap/../../../../private/pacbio/barcoding/bin/lima /rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/raw/0700055A/isoseq3/ccs.bam /rds/general/user/ph323/home/MRes.project.2/codes/1_preprocessing/primer.fasta /rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/raw/0700055A/isoseq3/new_fl.bam --isoseq
> @PG   ID:refine       VN:4.1.2 (commit v4.1.2)        CL:refine /rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/raw/0700055A/isoseq3/new_fl.primer_5p--primer_3p.bam /rds/general/user/ph323/home/MRes.project.2/codes/1_preprocessing/primer.fasta /rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/raw/0700055A/isoseq3/new_flnc.bam
> @PG   ID:samtools     PN:samtools     PP:refine       VN:1.18 CL:samtools view -H new_flnc.bam
```
Okay now we can confirm that there is two additional tags @HD and @RG in the pacbio IsoSeq3 pipeline produced bam file compared to the bam file produced by BGI custom pipeline. We therefore ran the PacBio IsoSeq3 pipeline on all the samples and transferred the corresponding @HD and @RG from each new_flnc.bam to our previously BGI created isoseq_flnc.bam. After that, when we perform this code again:      
```
$isoseq cluster2 "$output_fofn" "$clustered_bam" 
``` 
We succesfully created a clustered bam file.   


## Issue 1: PacBio IsoSeq3 Lima flagged high number of undesired 3p--3p pairs (~70-80%) in x.lima.summary

### Context:    
After generating hifi_reads.bam (or ccs.bam) from subreads.bam using `ccs`, we proceeded our analysis with `isoseq3 lima`. The quality of ccs.bam looks good as on average around 98% of ccs.bam reads passed Q30. 

### Problem:    
However, upon checking the quality we found that there is a substantial amount of reads not passing quality control (70~90%), mostly due to the reason of undesired 3p--3p pairs, an example report provided below:    
```
ZMWs input                (A) : 138055
ZMWs above all thresholds (B) : 25221 (18.27%)
ZMWs below any threshold  (C) : 112834 (81.73%)

ZMW marginals for (C):
Below min length              : 22 (0.02%)
Below min score               : 0 (0.00%)
Below min end score           : 32874 (29.13%)
Below min passes              : 9503 (8.42%)
Below min score lead          : 0 (0.00%)
Below min ref span            : 31899 (28.27%)
Without SMRTbell adapter      : 3 (0.00%)
Undesired 5p--5p pairs        : 8765 (7.77%)
Undesired 3p--3p pairs        : 85895 (76.13%)
Undesired no hit              : 3 (0.00%)

ZMWs for (B):
With different pair           : 25221 (100.00%)
Coefficient of correlation    : 0.00%

ZMWs for (A):
Allow diff pair               : 128552 (93.12%)
Allow same pair               : 128552 (93.12%)

Reads for (B):
Above length                  : 25221 (100.00%)
Below length                  : 0 (0.00%)
```
When looking through trouble shooting documents by PacBio we found that this issue was described in early SMRT link trouble shooting files (version 6, 7, and 8 only). Specically, they described it as `Iso-Seq reads flagged as 3’--3’ reads in classification step`, with symptoms described matching our case:    
  Symptoms:     
  • The majority of Iso-Seq reads are failing in the classification step.    
  • Lower number of transcripts from Isoseq3.    
  • The total number of CCS reads look good, but the resulting number of reads with 5' and 3' primer is much lower.    
Following their trouble shooting guide, we inspected ccs.bam > ccs.fasta for a 3'ATGGG overhang of the 5'primer. Typical CCS reads of a proper IsoSeq library should contain: `5' primer → ATGGGG overhang → cDNA sequence → polyA tail → 3' primer` and have 5-6 bp of ATGGGG overhang missing from the 3' end of the 5' primer sequence.     
Instead, what we observed is that ATGGGG overhang are missing here and there and have poly A/T tracks on both side of a "transcript", indicating the presence of random priming. Moreover, different from the problem described in their documentation, we also cannot find our 3'end primer in ccs reads, we counted the number of 5'end primer `AAGCAGTGGTATCAACGCAGAGTACATGGGG` and found them in 276110 sequences where as the 3'end primer `GTACTCTGCGTTGATACCACTGCTTACTAGT` was only found in at the end of 6 sequences in an example ccs.bam. We do see a random insertion of inverse 5'end primer with ATGGGG removed `GTACTCTGCGTTGATACCACTGCTT` (partial 3' primer) randomly inserted across the ccs reads (not locating at the sequence end as they should be).     
Trouble shooting documentation suggests that this may be a problem in sample preparation, specifically something might've gone wrong during the library preparation step due to a degraded Template Switching Oligonucleotides (TSO) batch or lack of TSO. 

### Solution    
The proposed solution by trouble shooting guide 













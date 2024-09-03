#!/bin/bash

## Meta directories ------------------------------------------------------------
cd "~/MRes.project.2/codes/0_packages/deeploc2/"
tappAS="/rds/general/user/ph323/ephemeral/MRes.project.2/raw_data/tappAS"
CHUNKS="$tappAS/transcoder/chunks"

## Deeploc2 
chunk_fasta="$CHUNKS/chunk_0.fasta"
deeploc2 -f $chunk_fasta -o "$tappAS/deeploc" -m "Fast"


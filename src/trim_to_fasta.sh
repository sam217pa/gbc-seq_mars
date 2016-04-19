#!/usr/bin/env bash

seqret -sformat fastq -osformat fasta -auto -stdout \
       -sequence data/trimmed.fastq > data/trimmed.fasta

fastx_reverse_complement \
    -i raw/ref/reference1073bis-1392.fasta \
    >> data/trimmed.fasta

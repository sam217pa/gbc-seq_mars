#!/bin/bash

# transforme en fastq
cd data/ws

touch ws_untrimmed.fastq

for file in spectro/*.ab1
do
    seqret \
        -sformat abi \
        -osformat fastq \
        -auto \
        -stdout \
        -sequence $file \
        >> ws_untrimmed.fastq
done

cd ../sw

# TEMPORAIRE, élimine le fichier .ab1 tendancieux pour l'instant
rm spectro/psw76-1073bis.ab1

touch sw_untrimmed.fastq

for file in spectro/*.ab1
do
    seqret \
        -sformat abi \
        -osformat fastq \
        -auto \
        -stdout \
        -sequence $file \
        >> sw_untrimmed.fastq
done

#!/bin/bash

cd data/ws

mkdir fasta seq

find . -name "*.fas"  -maxdepth 1 -exec mv {} ./fasta/ \;
find . -name "*.ab1"  -maxdepth 1 -exec mv {} ../spectro/ \;
find . -name "*.seq"  -maxdepth 1 -exec mv {} ./seq/ \;
find . -name "*.csv"  -maxdepth 1 -exec mv {} ../csv/ \;

cd ../sw

mkdir fasta seq

find . -name "*.fas"  -maxdepth 1 -exec mv {} ./fasta/ \;
find . -name "*.ab1"  -maxdepth 1 -exec mv {} ../spectro/ \;
find . -name "*.seq"  -maxdepth 1 -exec mv {} ./seq/ \;
find . -name "*.csv"  -maxdepth 1 -exec mv {} ../csv/ \;

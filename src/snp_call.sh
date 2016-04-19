   #!/usr/bin/env bash

   # Variant calling using ssaha2 and ssahaSNP

   if [ ! -e data/snp_call ]
   then
       mkdir data/snp_call
   fi

   cd data/snp_call

   if [[ ! -e trimmed.fastq && ! -e ref.fasta ]]
   then
       ln -s    ../trimmed.fastq .
       ln -s ../../raw/ref/reference1073bis-1392.fasta ref.fasta
   fi

   ## alignement à la séquence de référence
   #' -output psl :             format de sortie psl
   #' reference_reverse.fasta : séquence de référence
   #' trimmed.fastq :           séquence à aligner
   #' output.psl :              fichier de sortie
   ../../bin/ssahaSNP/ssaha2 -output psl ref.fasta trimmed.fastq > output.psl

   ## column annotation based on
   ## ftp://ftp.sanger.ac.uk/pub/resources/software/ssahasnp/readme,
   ## part (6) some further information
   # la première ligne du fichier .dat, afin d'être lu dans R
   cat \
       <( echo "match subject_name index_of_subject read_name s_base q_base s_qual q_qual offset_on_subject offset_on_read length_of_snp start_match_of_read end_match_of_read match_direction length_of_subject" ) \
       <( ../../bin/ssahaSNP/ssahaSNP ref.fasta trimmed.fastq | \
                egrep ssaha:SNP | \
                awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}') \
       > snp_calling.dat

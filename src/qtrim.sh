  #!/bin/bash

  #' -qtrim=rl : quality trim right and left
  #' -trimq=28 : trim if quality < 28 (sanger encoding, illumina 1.9)
  #' -minlen=620 : keep only seq with length > 620, after trimming.
  #' -Xmx1g : tells bbduk / java to use 1G of RAM

  FASTQ="data/sw/sw_untrimmed.fastq
  data/ws/ws_untrimmed.fastq
  "

  for f in $FASTQ
  do
      out=$f.trim
      bin/bbmap/bbduk.sh -Xmx1g \
                         -in=$f \
                         -out=$f.qtrim \
                         -qtrim=rl \
                         -trimq=28 \
                         -minlen=600

      ## convertit les bases d'une qualité inférieure à 20 en N.
      seqtk seq \
            -q20 \
            -nN \
            $f.qtrim > $out

      ## convertit le fastq en fasta
      seqret \
          -sformat fastq \
          -osformat fasta \
          -auto -stdout \
          -sequence $out > $f.fasta

      rm $f.qtrim
  done

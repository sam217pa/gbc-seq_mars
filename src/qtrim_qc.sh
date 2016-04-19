  #!/usr/bin/env bash

  cd data

  TRIM="ws/ws_untrimmed.fastq.trim sw/sw_untrimmed.fastq.trim"

  mkdir qtrim_qc

  for f in $TRIM
  do
      fastqc $f -o qtrim_qc
  done

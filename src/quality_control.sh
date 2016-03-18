#!/bin/bash

# contrôle qualité des séquences.
cd data

# crée un dossier qc de quality control pour ws et sw
mkdir -p qc/{ws,sw}

# fastqc sortie dasn le dossier approprié
fastqc ws/ws_untrimmed.fastq -o qc/ws
fastqc sw/sw_untrimmed.fastq -o qc/sw

# transfert du fichier .png qui nous intéresse dans le dir adapté
# et suppression des dossiers inutiles
#' -qq : very quiet
#' -o  : overwrite without warning
#' -d  : out dir
unzip -qq -o \
      qc/sw/sw_untrimmed_fastqc.zip \
      sw_untrimmed_fastqc/Images/per_base_quality.png \
      -d ../anl/ \
    && mv ../anl/sw_untrimmed_fastqc/Images/per_base_quality.png \
          ../anl/sw_per_base_qual.png \
    && rm -r ../anl/sw_untrimmed_fastqc

# transfert du fichier .png qui nous intéresse dans le dir adapté
# et suppression des dossiers inutiles
unzip -qq -o \
      qc/ws/ws_untrimmed_fastqc.zip \
      ws_untrimmed_fastqc/Images/per_base_quality.png \
      -d ../anl/ \
    && mv ../anl/ws_untrimmed_fastqc/Images/per_base_quality.png \
          ../anl/ws_per_base_qual.png \
    && rm -r ../anl/ws_untrimmed_fastqc

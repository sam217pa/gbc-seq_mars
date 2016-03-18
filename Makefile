# range les fichiers classés par extension

.PHONY: clean

anl/sw_per_base_qual.png: data/sw/sw_untrimmed.fastq
		bash src/quality_control.sh

data/sw/sw_untrimmed.fastq:
		bash src/ab1_to_fastq.sh

# met en place la structure de dossier
raw/dir: clean
		mkdir -p data/{ws,sw,csv} anl
		unzip -qq -o raw/1582203.zip -d data/ws/
		unzip -qq -o raw/1582443.zip -d data/sw/
		bash src/sort_into_dir.sh
		find data/ -name "*.pdf" -maxdepth 2 -exec mv {} anl/ \;
		touch raw/dir

# supprime tout sauf le README.org
clean:
		rm -rf data anl

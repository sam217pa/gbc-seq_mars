# range les fichiers classés par extension

.PHONY: clean



# met en place la structure de dossier
raw/dir: clean
		mkdir data anl data/ws data/sw data/csv
		unzip -qq -o raw/1582203.zip -d data/ws/
		unzip -qq -o raw/1582443.zip -d data/sw/
		bash src/sort_into_dir.sh
		find data/ -name "*.pdf" -maxdepth 2 -exec mv {} anl/ \;
		touch raw/dir

# supprime tout sauf le README.org
clean:
		rm -rf data anl

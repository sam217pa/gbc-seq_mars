# range les fichiers classés par extension
vpath %.sh src
vpath %.py src
vpath %.png anl

## ANALYSES

.PHONY: gatc_snpcall
gatc_snpcall: README.org | data/ws
	R --vanilla -e "rmarkdown::render('anl/snp_call/snp_call.r')"

# analyse la qualité via fastqc
.PHONY: check_qual
check_qual: sw_per_base_qual.png

anl/sw_per_base_qual.png: quality_control.sh data/sw/sw_untrimmed.fastq
	bash $<

# converti les spectrogrammes en fastq
data/sw/sw_untrimmed.fastq: ab1_to_fastq.sh | data/ws
	bash $<

# met en place la structure de dossier
data/ws: sort_into_dir.sh | clean
	mkdir -p data/{ws,sw,csv} anl
	unzip -qq -o raw/1582203.zip -d data/ws/
	unzip -qq -o raw/1582443.zip -d data/sw/
	bash $<
	find data/ -name "*.pdf" -maxdepth 2 -exec mv {} anl/ \;

src/*.sh: README.org
	orgtangle $<

README.org:
	orgtangle $@

.PHONY: commit
commit:
	git stage README.* Makefile
	git commit -m "sauvegarde"
	git push

# supprime tout sauf le README.org
.PHONY: clean
clean:
	rm -rf data anl

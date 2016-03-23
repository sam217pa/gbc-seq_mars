# range les fichiers classés par extension
vpath %.sh src test
vpath %.py src
vpath %.png anl

## ANALYSES

.PHONY: all
all: check_qual gatc_snpcall

## snp calling via ssahaSNP
.PHONY: snp_call
snp_call: qtrim


## quality trim
.PHONY: qtrim
qtrim: data/sw/sw_untrimmed.fastq.trim

data/sw/sw_untrimmed.fastq.trim: qtrim.sh data/sw/sw_untrimmed.fastq
	bash $<


## SNP calling de gatc
.PHONY: gatc_snpcall
gatc_snpcall: tangle | data
	R --vanilla -e "rmarkdown::render('anl/gatc_snpcall/snp_call.r')"
	git stage anl/gatc_snpcall/snp_call.html

# analyse la qualité via fastqc
.PHONY: check_qual
check_qual: sw_per_base_qual.png

sw_per_base_qual.png: quality_control.sh data/sw/sw_untrimmed.fastq
	bash $<

# converti les spectrogrammes en fastq
data/sw/sw_untrimmed.fastq: ab1_to_fastq.sh data/sw/seq
	bash $<

# met en place la structure de dossier
data/sw/seq: sort_into_dir.sh | data
	bash $<

.PHONY: tangle
tangle: README.org | data
	orgtangle $<

data:
	mkdir -p data/{ws,sw,csv} anl/{gatc_snpcall,snp_call}
	unzip -qq -o raw/1582203.zip -d data/ws/
	unzip -qq -o raw/1582443.zip -d data/sw/
	find data/ -name "*.pdf" -maxdepth 2 -exec mv {} anl/ \;

%.sh: tangle
%.r:  tangle

.PHONY: test
test: test.sh
	bash $<

.PHONY: commit
commit:
	git stage README.* Makefile
	git commit -m "sauvegarde"
	git push


# supprime tout sauf le README.org
.PHONY: clean cleanall
clean:
	rm -rf data

cleanall:
	rm -rf data anl

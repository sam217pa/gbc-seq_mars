# range les fichiers classés par extension
vpath %.sh src test
vpath %.py src
vpath %.png anl

## ANALYSES

.PHONY: all
all: phruscle_call snp_call check_qual gatc_snpcall

## snp calling via phruscle
.PHONY: phruscle_call
phruscle_call: anl/phruscle_call/phruscle_call.html

anl/phruscle_call/phruscle_call.html: anl/phruscle_call/phruscle_call.R data/phruscle_snpcall.csv
	R --vanilla -e "rmarkdown::render('$<')"
	git stage $< $@

data/phruscle_snpcall.csv: src/phruscler phruscle.py src/make_ref data/sw/seq
	mkdir data/ws/aln
	mkdir data/sw/aln
	mkdir data/s_w/aln
	src/make_ref
	$<

.PHONY: clean_phruscle
phruscle_clean:
	rm -r data/sw/aln data/ws/aln data/s_w/aln data/phruscle_snpcall.csv
## </phruscle>

## snp calling via ssahaSNP
# .PHONY: snp_call
# snp_call: qtrim
.PHONY: snp_call
snp_call: anl/snp_call/snp_call.r data/snp_call/snp_calling.dat data/id_table.dat
	R --vanilla -e "rmarkdown::render('$<')"
	git stage $<
	git stage anl/snp_call/snp_call.html

data/id_table.dat: make_id_table.py data/sw/seq
	python $< > $@

data/snp_call/snp_calling.dat: snp_call.sh data/trimmed.fastq
	bash $<

data/trimmed.fasta: trim_to_fasta.sh data/trimmed.fastq
	bash $<

data/trimmed.fastq: pool_trim.sh data/sw/sw_untrimmed.fastq.trim
	bash $<

data/qtrim_qc: qtrim_qc.sh data/sw/sw_untrimmed.fastq.trim
	bash $<
	git stage -f data/qtrim_qc/*.html

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
data/sw/seq: | sort_into_dir.sh data
	bash src/sort_into_dir.sh

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

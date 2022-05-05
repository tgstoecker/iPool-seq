# Makefile, Copyright 2017 Florian G. Pflug
# 
# This file is part of the iPool-Seq Analysis Pipeline
#
# The iPool-Seq Analysis Pipeline is free software: you can redistribute it
# and/or modify it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# The iPool-Seq Analysis Pipeline is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the iPool-Seq Analysis Pipeline.  If not, see
# <http://www.gnu.org/licenses/

# Options
TRIMMOMATIC_OPTS ?= -threads 4
NGM_OPTS ?= -t 4
TRIM_CORES ?= 4
UMICOUNTS_MAPQ ?= 20
SCRATCHDIR ?= /tmp

# Software
NGM ?= ngm
TRIMMOMATIC ?= trimmomatic
SAMTOOLS ?= samtools
GZIP ?= gzip
PICARD ?= picard
RSCRIPT ?= Rscript
PYTHON ?= python

# UTILS used by various recipes
define run-utils
	set -o pipefail; \
	set -e; \
	die () { \
		echo "$$*" >&2; \
		exit 1; \
	}; \
	outputExists () { \
		for f in $$*; do \
			if [ -e "$$f" ]; then \
				if [ "$$REMAKE" != "y" ]; then \
					die "ERROR: $$f already exists. Set RMAKE=y to remake"; \
				else \
					log "*** Will remake $$f"; \
					rm "$$f"; \
				fi; \
			fi; \
		done \
	}; \
	log() { \
		echo "$$*" >&2; \
	}
endef

# NGM
define run-ngm
	$(run-utils); \
	NGM_LOG="$(patsubst %.bam,%.log,$@)"; \
	NGM_RAW="$(SCRATCHDIR)/tmp-ngm-raw-$(subst /,:,$@)"; \
	NGM_CLEAN="$(SCRATCHDIR)/tmp-ngm-cleaned-$(subst /,:,$@)"; \
	NGM_INVAL="$(patsubst %.bam,%.invalid.txt,$@)"; \
	outputExists "$@"; \
	( \
		if [ "$(word 4,$^)" != "" ]; then \
			log "*** Running NGM on pair-ended inputs $(word 3,$^) and $(word 4,$^) with parameters from $(word 1,$^) "; \
			log "*** Parameters: $$(cat "$(word 1,$^)")"; \
			$(NGM) -r "$(word 2,$^)" -1 "$(word 3,$^)" -2 "$(word 4,$^)" --bam -o "$${NGM_RAW}"\
				$$(cat "$(word 1,$^)") $(NGM_OPTS); \
		else\
			log "*** Running NGM on unpaired input $(word 3,$^)"; \
			log "*** Parameters: $$(cat "$(word 1,$^)")"; \
			$(NGM) -r "$(word 2,$^)" -q "$(word 3,$^)" --bam -o "$${NGM_RAW}" \
				$$(cat "$(word 1,$^)") $(NGM_OPTS); \
		fi;\
		log "*** Scanning $${NGM_RAW} for invalid reads"; \
		$(PICARD) ValidateSamFile INPUT="$${NGM_RAW}" \
			MAX_OUTPUT=1000000000 IGNORE=RECORD_MISSING_READ_GROUP IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND \
			TMP_DIR="$(SCRATCHDIR)" \
			| sed -n 's/^.*Read name \([^,]*\),.*$$/\1/p' \
			| sort \
			| uniq \
			> "$${NGM_INVAL}" \
			|| true; \
		log "*** Found the following $$(wc -l < "$${NGM_INVAL}") invalid reads"; \
		cat "$${NGM_INVAL}"; \
		if [ $$(wc -l < "$${NGM_INVAL}") -gt 0 ]; then \
			log '*** Removing invalid reads, writing $@ sorted by read name'; \
			$(PICARD) FilterSamReads INPUT="$${NGM_RAW}" OUTPUT="$${NGM_CLEAN}" READ_LIST_FILE="$${NGM_INVAL}" \
				FILTER=excludeReadList VALIDATION_STRINGENCY=SILENT \
				WRITE_READS_FILES=false SORT_ORDER=queryname \
				TMP_DIR="$(SCRATCHDIR)"; \
		else \
			log "*** Sorting $@ by read name"; \
			$(PICARD) SortSam INPUT="$${NGM_RAW}" OUTPUT="$${NGM_CLEAN}" SORT_ORDER=queryname \
				TMP_DIR="$(SCRATCHDIR)"; \
		fi; \
		rm "$${NGM_RAW}"; \
		log "*** Fixing mate information in $@"; \
		$(PICARD) FixMateInformation INPUT="$${NGM_CLEAN}" OUTPUT="$@" \
			TMP_DIR="$(SCRATCHDIR)"; \
		rm "$${NGM_CLEAN}"; \
	) 2>&1 | tee "$${NGM_LOG}"
endef

# Don't delete intermediate files
IDS=$(patsubst data/%/,%,$(filter data/%/,$(wildcard data/*/)))
.SECONDARY:	$(foreach id,$(IDS),\
		data/$(id)/raw.1.fq.gz \
		data/$(id)/raw.2.fq.gz \
		data/$(id)/tom.1.fq.gz \
		data/$(id)/tom.2.fq.gz \
		data/$(id)/trim.1.fq.gz \
		data/$(id)/trim.2.fq.gz \
		data/$(id)/ngm.bam \
		data/$(id)/ngm.sorted.bam \
		data/$(id)/ngm.assigned.bam \
		data/$(id)/ngm.counts.tab.gz \
		data/$(id)/ngm.results.rda \
		data/$(id)/ngm.stats.rda)

# Disable built-in rules
.SUFFIXES:

# Print variables
print-%  : ; @echo $* = $($*)

%.1.fq.gz %.2.fq.gz: %.bam
	@$(run-utils); outputExists "$(*).1.fq.gz" "$(*).2.fq.gz" "$(*).1.fq" "$(*).2.fq"; \
	READ1_FQ="$(SCRATCHDIR)/tmp-split-$(subst /,:,$(*)).1.fq"; \
	READ2_FQ="$(SCRATCHDIR)/tmp-split-$(subst /,:,$(*)).2.fq"; \
	echo "*** Splitting $(*).bam into $${READ1_FQ} and $${READ2_FQ}"; \
	$(SAMTOOLS) fastq -1 "$${READ1_FQ}" -2 "$${READ2_FQ}" "$(*).bam" || exit 1; \
	echo "Compressing $${READ1_FQ} and $${READ2_FQ} to $(*).1.fq and $(*).2.fq"; \
	$(GZIP) --fast --stdout "$${READ1_FQ}" > "$(*).1.fq.gz" & READ1_PID=$$!; \
	$(GZIP) --fast --stdout "$${READ2_FQ}" > "$(*).2.fq.gz" & READ2_PID=$$!; \
	ERROR=0; \
	wait "$${READ1_PID}" || ERROR=$$?; \
	wait "$${READ2_PID}" || ERROR=$$?; \
	rm -f "$${READ1_FQ}" "$${READ2_FQ}"; \
	exit $$ERROR;

%/tom.1.fq.gz %/tom.2.fq.gz %/tom.1up.fq.gz %/tom.2up.fq.gz %/tom.log: %/raw.1.fq.gz %/raw.2.fq.gz %/tom.cfg
	@$(run-utils); outputExists "$(*)/tom.1.fq.gz" "$(*)/tom.2.fq.gz"; \
	SCRATCH="$(SCRATCHDIR)/tmp-tom-$(subst /,:,$(*))"; \
	echo "*** Running Trimm-O-Matic on $(*)/raw.1.fq.gz and $(*)/raw.2.fq.gz"; \
	. "$(*)/tom.cfg"; \
	echo "*** Trim Commands: $${CMD}"; \
	$(TRIMMOMATIC) PE $(TRIMMOMATIC_OPTS) \
		"$(*)/raw.1.fq.gz" "$(*)/raw.2.fq.gz" \
		"$${SCRATCH}.1.fq.gz" "$${SCRATCH}.1up.fq.gz" "$${SCRATCH}.2.fq.gz" "$${SCRATCH}.2up.fq.gz" \
		$${CMD} \
		2>&1 | tee "$(*)/tom.log" \
		|| exit 1; \
	echo "*** Synthesiying dummy mates for unpaired output"; \
	(zcat "$${SCRATCH}.1.fq.gz"; \
	 zcat "$${SCRATCH}.1up.fq.gz"; \
	 zcat "$${SCRATCH}.2up.fq.gz" | sed -n 's|^@\(.*\)/2$$|@\1/1\nN\n+\n!|p' \
	) | $(GZIP) -c > "$(*)/tom.1.fq.gz" & READ1_PID=$$!; \
	(zcat "$${SCRATCH}.2.fq.gz"; \
	 zcat "$${SCRATCH}.1up.fq.gz" | sed -n 's|^@\(.*\)/1$$|@\1/2\nN\n+\n!|p'; \
	 zcat "$${SCRATCH}.2up.fq.gz"; \
	) | $(GZIP) -c > "$(*)/tom.2.fq.gz" & READ2_PID=$$!; \
	wait "$${READ1_PID}" || exit 1; \
	wait "$${READ2_PID}" || exit 1; \
	if [ $$(zcat "$(*)/tom.1.fq.gz" | wc -l) -ne $$(zcat "$(*)/tom.2.fq.gz" | wc -l) ]; then \
	  echo "Output not properly paired" >&2; \
	  exit 1; \
	fi

%/trim.1.fq.gz %/trim.2.fq.gz: %/tom.1.fq.gz %/tom.2.fq.gz
	@$(run-utils); outputExists "$(*)/trim.1.fq.gz" "$(*)/trim.2.fq.gz"; \
	echo "*** Running trim.tag.py on $(*)/tom.1.fq.gz and $(*)/tom.2.fq.gz"; \
	$(PYTHON) ./trim.tag.py $(TRIM_CORES) \
		"$(*)/tom.1.fq.gz" "$(*)/tom.2.fq.gz" \
		"$(*)/trim.1.fq.gz" "$(*)/trim.2.fq.gz" \
		2>&1 | tee "$(*)/trim.log" \
		|| exit 1;

%/ngm.bam: %/ngm.cfg %/ref.fasta %/trim.1.fq.gz %/trim.2.fq.gz
	@$(run-ngm)

%.sorted.bam %.sorted.bai: %.bam
	@$(run-utils); outputExists "$(*).sorted.bam"; \
	$(PICARD) SortSam INPUT="$(*).bam" OUTPUT="$(*).sorted.bam" SORT_ORDER=coordinate \
		TMP_DIR="$(SCRATCHDIR)" COMPRESSION_LEVEL=1 \
		|| exit 1; \
	$(PICARD) BuildBamIndex INPUT="$(*).sorted.bam" || exit 1; \

%.assigned.bam %.assigned.bai: %.sorted.bam
	@$(run-utils); outputExists "$(*).assigned.bam"; \
	FEATURES="$$(dirname "$(word 1,$^)")/features.gff"; \
	$(PYTHON) ./assign_to_features.py "$${FEATURES}" "$(*).sorted.bam" "$(*).assigned.bam"
	$(PICARD) BuildBamIndex INPUT="$(*).assigned.bam" || exit 1; \

%.counts.tab.gz: %.assigned.bam
	@$(run-utils); outputExists "$(*).counts.tab.gz"; \
	$(PYTHON) ./umicounts.tag.py $(UMICOUNTS_MAPQ) "$(*).assigned.bam" "$(*).counts.tab.gz" 2>&1 | tee "$(*).counts.log"

%.results.rda: %.counts.tab.gz %.results.cfg
	@$(run-utils); outputExists "$(*).results.rda"; \
	$(RSCRIPT) ./counts2results.R "$(*).counts.tab.gz" "$(*).results.cfg" "$(*).results.rda" 2>&1 | tee "$(*).results.log"

%/ngm.stats.rda: %/raw.bam %/trim.1.fq.gz %/ngm.bam %/ngm.assigned.bam %/ngm.counts.tab.gz %/ngm.results.cfg %/ngm.results.rda
	@$(run-utils); outputExists "$(*)/ngm.stats.rda"; \
	$(RSCRIPT) ./collect_stats.R "$(*)" 2>&1 | tee "$(*)/ngm.stats.log"

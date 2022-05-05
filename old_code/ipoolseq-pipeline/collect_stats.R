#!/usr/bin/env Rscript
# collect_stats.R, Copyright 2017 Florian G. Pflug
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
# You should have received a copy of the GNU Affero General Public License
# along with the iPool-Seq Analysis Pipeline.  If not, see
# <http://www.gnu.org/licenses/

# *****************************************************************************
# ./collect_stats.R <directory>
# The directory must contain the follwing files (bracket states the program
# that creates each file, and a brief description)
#
#   raw.bam                    (raw library reads as published on ENA)
#   ngm.results.cfg            (config file for umicounts.tag)
#   trim.1.fq.gz, trim.2.fq.gz (trim.tag.py reads sans UMIs & adapters)
#   ngm.bam                    (NGM:, mapped reads)
#   ngm.assigned.bam           (assign_to_features: reads in RG per eff.flank)
#   ngm.counts.tab.gz          (umicounts.tag: UMIs per eff.flank)
# *****************************************************************************

library(data.table)

ARGV <- commandArgs(trailingOnly = TRUE)
sampledir <- ARGV[1]

readcount.bam <- function(file, flags_set=0, flags_cleared=0) {
  as.integer(system(paste0('samtools view -f ', flags_set, ' -F ', flags_cleared, ' -c ', file),
                    intern=TRUE))
}

fragcount.bam <- function(file, flags_set=0, flags_cleared=0) {
  as.integer(system(paste0('samtools view -f ', flags_set, ' -F ', flags_cleared, ' ', file,
                           ' | cut -f 1 | sort --buffer-size=2G | uniq | wc -l | cut -d " " -f 1'),
                    intern=TRUE))
}

fragcount.assigned.bam <- function(file, flags_set=0, flags_cleared=0) {
  as.integer(system(paste0('samtools view -f ', flags_set, ' -F ', flags_cleared, ' ', file,
                           ' | egrep -v "RG:Z:(unmatched|ambiguous)"',
                           ' | cut -f 1 | sort --buffer-size=2G | uniq | wc -l | cut -d " " -f 1'),
                    intern=TRUE))
}

fragcount.ambiguous.bam <- function(file, flags_set=0, flags_cleared=0) {
  as.integer(system(paste0('samtools view -f ', flags_set, ' -F ', flags_cleared, ' ', file,
                           ' | egrep "RG:Z:ambiguous"',
                           ' | cut -f 1 | sort --buffer-size=2G | uniq | wc -l | cut -d " " -f 1'),
                    intern=TRUE))
}

fragcount.unmatched.bam <- function(file, flags_set=0, flags_cleared=0) {
  as.integer(system(paste0('samtools view -f ', flags_set, ' -F ', flags_cleared, ' ', file,
                           ' | egrep "RG:Z:unmatched"',
                           ' | cut -f 1 | sort --buffer-size=2G | uniq | wc -l | cut -d " " -f 1'),
                    intern=TRUE))
}

readcount.fq <- function(file) {
  as.integer(system(paste0('zcat ', file, ' | wc -l  | cut -d " " -f 1'), intern=TRUE)) / 4
}

STATS <- list()

message("Counting total number of sequenced fragments")
STATS$frags.total <- fragcount.bam(paste0(sampledir, '/raw.bam'))
message("Counting total number of valid fragments (containing a molecular barcode)")
STATS$frags.trim <- readcount.fq(paste0(sampledir, '/trim.1.fq.gz'))
message("Counting total number of mapped fragments")
STATS$frags.mapped <- fragcount.bam(paste0(sampledir, '/ngm.bam'), 0, 4)
message("Counting total number of unassigned fragments (mapped but not assigned to a knockout)")
STATS$frags.unassigned <- fragcount.unmatched.bam(paste0(sampledir, '/ngm.assigned.bam'))
message("Counting total number of ambigioulsy assigned fragments (mapped and assigne to multiple knockouts)")
STATS$frags.ambiguous <- fragcount.ambiguous.bam(paste0(sampledir, '/ngm.assigned.bam'))
message("Counting total number of uniquely assigned fragments (mapped and assigne to multiple knockouts)")
STATS$frags.assigned <- fragcount.assigned.bam(paste0(sampledir, '/ngm.assigned.bam'))

# Read per-sample config
message("Loading per-sample config file ", paste0(sampledir, "/ngm.results.cfg"))
source(paste0(sampledir, "/ngm.results.cfg"))
CFG <- list(threshold=threshold, molecules=molecules)

message("Loading UMI clustering results ", paste0(sampledir, '/ngm.counts.tab.gz'))
umis <- data.table(read.table(paste0(sampledir, '/ngm.counts.tab.gz'), header=TRUE))
umis.assigned <- umis[!(rg %in% c("unmatched", "ambiguous")), ]
stopifnot(sum(umis.assigned$count) != STATS$frags.assigned)

message("Counting raw UMIs, UMIs after clustering, and UMIs after filtering")
STATS$umis.raw <- sum(umis.assigned$rawumis)
STATS$umis.clustered <- nrow(umis.assigned[count > 0,])
STATS$umis.filtered <- nrow(umis.assigned[count >= threshold,])

message("Counting total number of fragments belonging to a true (i.e. post-filtering) UMI")
STATS$frags.trueumis <- umis.assigned[count >= threshold, sum(count)]

message("Saving stats to ", paste0(sampledir, "/ngm.stats.rda"))
save(STATS, file=paste0(sampledir, "/ngm.stats.rda"))

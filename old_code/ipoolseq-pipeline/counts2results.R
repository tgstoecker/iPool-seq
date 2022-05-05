#!/usr/bin/env Rscript
# counts2results.R, Copyright 2017 Florian G. Pflug
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
# Implements "Correcting for artifacts and lost UMIs"
# ./counts2results.R <UMI table> <Config file> <Output>
# Output is written as an R data file (use extension .rda)
# *****************************************************************************
library(data.table)
library(gwpcR)
library(parallel)

CORES <- 4
ARGV <- commandArgs(trailingOnly = TRUE)
f.in <- ARGV[1]
f.cfg <- ARGV[2]
f.out <- ARGV[3]

# Load data. We use the readgroup to split into individual features
message("Loading ", f.in)
t.in <- data.table(read.table(f.in, header=TRUE))
t.in[, feature := rg ]

# Read per-sample config
source(f.cfg)
CFG <- list(threshold=threshold, molecules=molecules)
message("Using threshold of ", threshold, " reads/UMI, assuming ", molecules, " initial molecules")

# Fit global model
gm <- gwpcrpois.mom(mean(t.in[count >= threshold, count], na.rm=TRUE),
                    var(t.in[count >= threshold, count], na.rm=TRUE),
                    threshold=threshold, molecules=molecules,
                    nonconvergence.is.error=FALSE)
message("Overall efficiency ", round(100*gm$efficiency), "%, depth=", round(gm$lambda0, digits=3), " reads/UMI")

# Correct for unobserved threshold, store corrected counts in umis.tot
# XXX: Setting var.est.distfree=FALSE, use.nonconv.groupest=TRUE, obs.min.ingroup=2
# restores the behaviour of older gwpcR versions -- for new datasets, this isn't
# necessarily the best choice!
message("Fitting gene-wise model")
uc <- gwpcrpois.mom.groupwise(count ~ feature, t.in[count >= threshold,],
                              threshold=threshold, molecules=molecules,
                              ctrl=list(cores=CORES, verbose=TRUE, var.est.distfree=FALSE,
                                        obs.min.ingroup=1, use.nonconv.groupest=TRUE,
                                        use.nonconv.globalest=TRUE))

# Raw UMI frequencies
uf <- list()
dummy <- t.in[, {
  uf[[as.character(feature)]] <<- as.integer(count * 1.0) # Force copy!
}, by="feature"]

# Result
message("Saving results to ", f.out)
DATA <- list(umicounts=uc, umifreqs=uf, global=gm)
save(DATA, CFG, file=f.out)

# model.R, Copyright 2017 Florian G. Pflug
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

## WARNING: What is called "f" here, is 1/f in the model description
#           an model.pdf, and the same holds for "g". Should change
#           this here as well!

dmodel <- function(a, scale, disp, b, f, g, log=FALSE) {
#  print(cbind(a=a, scale=scale, disp=disp, b=b, f=f, g=g))
  dnbinom(a, mu=scale*b*g/f, size=b/(1+disp*b), log=log)
}

pmodel <- function(a, scale, disp, b, f, g, log.p=FALSE) {
  pnbinom(a, mu=scale*b*g/f, size=b/(1+disp*b), log.p=log.p)
}

qmodel <- function(p, scale, disp, b, f, g, log.p=FALSE) {
  qnbinom(p, mu=scale*b*g/f, size=b/(1+disp*b), log.p=log.p)
}

model.fit <- function(a, b, f, g) {
  if ((length(a) != length(b)) || (length(b) != length(f)) || (length(f) != length(g)))
    stop("a, b, f and g must have the same length")
  
  
  p.init <- c(scale=mean(((a * f) / (b * g))[b > 0], na.rm=TRUE),
              disp=1)
  ll.fun <- function(p) {
    if (all(p > 0))
      -sum(dmodel(a, p['scale'], p['disp'], b, f, g, log=TRUE))
    else NA
  }
  r <- optim(par=p.init, fn=ll.fun)
  if (r$convergence == 0)
    as.list(r$par)
  else
    stop("numerical optimazation failed to converge")
}

if (FALSE) {
  b.true <- outer(c(1,3), 10^c(0,0,1,1,2,2,2,3,3,3,4,4,4,5,5,5))
  g <- c(2.5, 1.5, 2, 1.2)[1 + 0:(length(b.true)-1) %% 4]
  b <- mapply(function (bp, gp) rpois(1, bp/gp), bp=b.true, gp=g)
  
  d <- c(0.6, 1.2, 0.8, 1.4, 1)[1 + 0:(length(b.true)-1) %% 4]
  
  scale <- 10
  a.true <- scale * b.true * d
  f <- c(0.6, 1.2, 1.4, 0.8)[1 + 0:(length(b.true)-1) %% 4]
  a <- mapply(function (ap, fp) rpois(1, ap/fp), ap=a.true, fp=f)
  
  print(model.fit(a[b > 0], b[b > 0], f[b > 0], g[b > 0]))
}

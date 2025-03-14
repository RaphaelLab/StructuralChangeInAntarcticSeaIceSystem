# This is the code for analysing the dynamic ARMA fits
# in Nature manuscript 2023-07-12845A
#
# "A Twenty-First Century Structural Change in Antarctica’s Sea ice System"
#  by Marilyn N. Raphael, Thomas J. Maierhofer, Ryan L. Fogt, William R. Hobbs, and Mark S. Handcock
#
# These standard R packages need to be present (and can be installed via 'install.packages("scales")', etc.
library(scales)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(grid)
library(gridExtra, warn.conflicts = FALSE)

library(rjags)
library(runjags, warn.conflicts = FALSE)

library(bayesplot)
library(loo)
#

nsim <- 2000 # Used in actual paper run
nsim <- 2 # To demonstrate the code
rs <- sample(1:2000, size=nsim, replace=FALSE)

get.posterior <- function(pars="ar1[", coda_samples, chain=1:4) {
    sin = NULL
    sel <- grep(pars,dimnames(coda_samples[[chain[1]]])[[2]],fixed=TRUE)
    for(ch in chain){
      sin = rbind(sin, coda_samples[[ch]][,sel])
    }
    return(sin)
}

library("future")
library("doRNG")
library("doFuture")
library("foreach")
doFuture::registerDoFuture()
future::plan(multicore, workers=30)
doRNG::registerDoRNG(1)

#source("Nature-s43247-025-02107-5-dynamic-ARIMA-plot-roots.R")
persistence <- matrix(0, nrow=nsim, ncol=1495)
companion <- matrix(0,ncol=12,nrow=12)
for(k in 1:11){companion[k+1,k] <- 1}
persistence <- foreach::foreach(sim = 1:nsim, .packages = c("runjags"), .export=c("rs"), .combine = 'rbind') %dopar% {
  runJagsOut <- readRDS(file=paste0("../Data/Nature-s43247-025-02107-5-dynamic-ARIMA.",sim,".rds"))
  coda_samples <- as.mcmc.list(runJagsOut$mcmc)

  ar1 <- get.posterior(pars="ar1[", coda_samples, chain=1:3)
  ar2 <- as.vector(t(get.posterior(pars="ar2", coda_samples, chain=1:3)))
  ar3 <- as.vector(t(get.posterior(pars="ar3", coda_samples, chain=1:3)))

  ar1 <- apply(ar1,2,mean)
  ar2 <- mean(ar2)
  ar3 <- mean(ar3)
  persistence.e <- rep(0, 1495)
  persistence.s <- rep(0, 1495)
    for(j in 1:1495){
      arlag <- c(ar1[j],0,0,ar2,0,0,0,0,0,0,0,ar3)
      companion[1,] <- arlag
      ar_roots <- eigen(companion)$values
      persistence.e[j] <- max(abs(ar_roots))
      persistence.s[j] <- 1/ (1-sum(arlag))
    }
  c(persistence.e, persistence.s)
}

tdate.m <- seq(1899+0.5/12,2025,by=1/12)[1:1495]
persistence.e <- persistence[,1:1495]
persistence.s <- persistence[,-c(1:1495)]

#save(persistence.e, persistence.s, file="persistence_m.RData")

pdf(paste0("../Figures/Fig3_persistence_brief.pdf",sep=""),paper="USr",width=11,height=9.5)

B=apply(persistence.e,2,quantile,c(0.05,0.95))
plot(x=tdate.m, y=apply(persistence.e,2,mean),type='l',xlab="Time (years)",ylab="persistence",
     ylim=c(min(B[1,]),max(B[2,])),
     xlim=c(1900,2025), main="", lwd=3, xaxt="n")
lines(x=tdate.m, y=B[1,],col=2,lty=2, lwd=3)
lines(x=tdate.m, y=B[2,],col=2,lty=2, lwd=3)
abline(h=1, lty=1, col="gray")
axis(1, at=c(seq(1900, 2020, by=5), 2024))
abline(h=0, lty=1, col="gray")
abline(v=c(1978.9, 1987.9, 2020, 2024.8), lty=1, col="gray")

B=apply(persistence.s,2,quantile,c(0.05,0.95))
plot(x=tdate.m, y=apply(persistence.s,2,mean),type='l',xlab="Time (years)",ylab="persistence",
 ylim=c(-1,1),
     xlim=c(1900,2025), main="", lwd=3, xaxt="n")
lines(x=tdate.m, y=B[1,],col=2,lty=2, lwd=3)
lines(x=tdate.m, y=B[2,],col=2,lty=2, lwd=3)
abline(h=1, lty=1, col="gray")
axis(1, at=c(seq(1900, 2020, by=5), 2024))
abline(h=0, lty=1, col="gray")
abline(v=c(1978.9, 1987.9, 2020, 2024.8), lty=1, col="gray")

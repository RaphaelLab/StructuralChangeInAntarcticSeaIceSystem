# This is the code for recreating the dynamic ARMA fits
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
#
# First we read in the NSIDC sea ice extent data 
# and the reconstructions of Maierhofer (2023) [See cite details in the manuscript]
load(file="../Data/nsidcV4.RData")
tempFile_location<- tempfile()
download.file(url="https://ucla.box.com/shared/static/ip0wlp72fslp40o7vudnwfvfwcbmsmg2",
              destfile=tempFile_location, method="wget", mode="rb")
reconstructions_total <- readRDS(tempFile_location)
file.remove(tempFile_location)
ls()

#=======================================================================
# Gather parameters and vectors for jags run
#=======================================================================

set.seed(2)

start <- 958
end <- dim(reconstructions_total)[2]
end

tdate.m <- seq(1899+0.5/12,2025,by=1/12)[1:end]
year.m <- trunc(tdate.m)
month.m <- rep(1:12,150)[1:end]

year.m[month.m==1] <- year.m[month.m==1] - 1
year.m[month.m==2] <- year.m[month.m==2] - 1

year.m[year.m < 1899] <- 1899

skip_run <- FALSE
nsim <- 200
nsim <- 2
rs <- sample(1:2000, size=nsim, replace=FALSE)

for(sim in 1:nsim){
sie <- reconstructions_total[rs[sim],,1]
tdate <- tdate.m[seq_along(sie)]
month <- month.m[seq_along(sie)]
year <- year.m[seq_along(sie)]

start_satellite <- 1
end_data <- length(sie) # 2023/12 plus 2024/1 to 2024/8

beta_poly_order <- 5
beta_poly <- poly(seq(along=tdate[(start_satellite:end_data)]), beta_poly_order)

jags_data = list(
            "sie" = sie[(start_satellite:end_data)],
            "tdate" = tdate[(start_satellite:end_data)],
            "mu_sigma" = c(log(0.5),rep(0, beta_poly_order)),
            "mu_abeta1" = rep(0, beta_poly_order+1),
            "tau_sigma" = diag(rep(0.5, beta_poly_order+1)),
            "tau_abeta1" = diag(rep(0.5, beta_poly_order+1)),
            "beta_basis" = cbind(rep(1, nrow(beta_poly)), beta_poly)
                )
#=======================================================================
# initialization
#=======================================================================
initfunction <- function(chain) {
  return(switch(chain,
                "1" = list(
                  .RNG.name="base::Super-Duper",
                  .RNG.seed=1),
                "2" = list(
                  .RNG.name="base::Super-Duper",
                  .RNG.seed=2),
                "3" = list(
                  .RNG.name="base::Super-Duper",
                  .RNG.seed=3),
                ))
}

model_file <-
"model {

  # Priors for population-level effects

  int_1 ~ dt(0, 1/(3)^2, 3) 
  sigma_0 ~ dnorm(0, 1) T(0.00001, )
  sigmabeta ~ dmnorm(mu_sigma, tau_sigma)
  abeta1 ~ dmnorm(mu_abeta1, tau_abeta1)
  mbeta1 ~ dmnorm(0, 0.5)
  ma1 ~ dunif(-1, 1)
  ma2 ~ dunif(-1, 1)
  ar2 ~ dunif(-1, 1)
  ar3 ~ dunif(-1, 1)
  ar10 ~ dunif(-1, 1)
  ar20 ~ dunif(-1, 1)

  # Apply autoregression to the residuals
  mu_[1] = int_1
  epsilon_[1] = sie[1]-mu_[1]
  sigma[1] = exp(inprod(sigmabeta, beta_basis[1,])) +sigma_0
  ar1[1] = inprod(abeta1, beta_basis[1,])
#
  mu_[2] = int_1 + ar1[2] * sie[2 - 1] +
                   ma1 * epsilon_[2 - 1]
  epsilon_[2] = sie[2]-mu_[2]
  sigma[2] = exp(inprod(sigmabeta, beta_basis[2,])) +sigma_0
  ar1[2] = inprod(abeta1, beta_basis[2,])
#
  for (i_ in 3:12) {
    mu_[i_] = int_1 + ar1[i_] * sie[i_- 2] +
                      ma1 * epsilon_[i_- 1]
    epsilon_[i_] = sie[i_]-mu_[i_]
    sigma[i_] = exp(inprod(sigmabeta, beta_basis[i_,])) +
     (tdate[i_] >= 1978.9) * (tdate[i_] < 1987.72) * sigma_0 
    ar1[i_] = inprod(abeta1, beta_basis[i_,])
  }
  for (i_ in 13:length(tdate)) {
    sigma[i_] = exp(inprod(sigmabeta, beta_basis[i_,])) +
     (tdate[i_] >= 1978.9) * (tdate[i_] < 1987.72) * sigma_0 
    ar1[i_] = inprod(abeta1, beta_basis[i_,])
    mu_[i_] = int_1 + 
       ar1[i_] * sie[i_ -  1] + 
       (tdate[i_] >= 1978.9) * (tdate[i_] < 1987.72) * ar10 * sie[i_ -  1] +
       ar2 * sie[i_ -  4] + 
       (tdate[i_] >= 1978.9) * (tdate[i_] < 1987.72) * ar20 * sie[i_ -  4] +
       ar3 * sie[i_ -  12] + 
       ma1 * epsilon_[i_ - 1] +
       ma2 * epsilon_[i_ - 12]
    epsilon_[i_] = sie[i_]-mu_[i_]
  }

  # Model and likelihood
  for (i_ in 1:length(tdate)) {
    # Likelihood and log-density for family = gaussian()
    sie[i_] ~ dnorm(mu_[i_], 1 / sigma[i_]^2)  # SD as precision
    loglik_[i_] = logdensity.norm(sie[i_], mu_[i_], 1 / sigma[i_]^2)  # SD as precision
  }
}"
#=======================================================================
# Conditional setup
#=======================================================================

nChains <- 3
# Next for long run
nAdaptSteps <- 50000
nBurninSteps <- 100000
nThinSteps <- 200
nUseSteps <- 150000

# Next for simple run
nAdaptSteps <- 500
nBurninSteps <- 5000
nThinSteps <- 20
nUseSteps <- 500

jags_params <- c("abeta1",
                 "mbeta1",
                 "sigma",
                 "ar1",
                 "ar2",
                 "ar3",
                 "ar10",
                 "ar20",
                 "sigma_0",
                 "ma1",
                 "ma2",
                 "loglik_"
                )
if(skip_run){
  runJagsOut <- readRDS(file=paste0("../Data/Nature-COMMSENV-s43247-023-01836-A-dynamic-ARIMA.",sim,".rds"))
} else {
  runJagsOut <- run.jags(method = "parallel",
                     model = model_file,
                     monitor = jags_params,
                     data = jags_data,
                     n.chains = nChains,
                     adapt = nAdaptSteps,
                     burnin = nBurninSteps,
                     sample = ceiling(nUseSteps/nChains),
                     thin = nThinSteps,
                     summarise = FALSE,
                     plots = FALSE,
                     inits = initfunction)
  saveRDS(runJagsOut, file=paste0("../Data/Nature-COMMSENV-s43247-023-01836-A-dynamic-ARIMA.",sim,".rds"))
}

#=======================================================================
# coda samples - MCMC
#=======================================================================
coda_samples <- as.mcmc.list(runJagsOut)

#=======================================================================
# combine mcmc draws for chains
#=======================================================================
get.posterior <- function(pars="ar1[", coda_samples, jags_data, chain=1:4) {
    sin = NULL
    sel <- grep(pars,dimnames(coda_samples[[chain[1]]])[[2]],fixed=TRUE)
    for(ch in chain){
      sin = rbind(sin, coda_samples[[ch]][,sel])
    }
    return(sin)
}

abeta1 <- get.posterior(pars="abeta1[", coda_samples, jags_data, chain=1:3)

# This code creates some simple plots
pdf(paste0("../Figures/Nature-COMMSENV-s43247-023-01836-A-dynamic-ARIM.pdf",sep=""),paper="USr",width=11,height=9.5)

ar1 <- get.posterior(pars="ar1[", coda_samples, jags_data, chain=1:3)
ar10 <- get.posterior(pars="ar10", coda_samples, jags_data, chain=1:3)
B=apply(ar1,2,quantile,c(0.05,0.95))
dev.off()
pdf("../Figures/Nature-COMMSENV-s43247-023-01836-A-dynamic-ARIMA-details.pdf",paper="USr",width=11,height=9.5)
range(jags_data$tdate)
plot(x=jags_data$tdate, y=sie,type='l',xlab="Time (years)",ylab="SIE",
     xlim=c(1900,2025), main="nsidc combined", lwd=3, xaxt="n")
plot(x=jags_data$tdate, y=jags_data$sie,type='l',xlab="Time (years)",ylab="SIE",
     xlim=c(1900,2025), main="sie standardized", lwd=3, xaxt="n")
vl <- length(apply(ar1,2,mean))
plot(x=jags_data$tdate, y=apply(ar1,2,mean)[ncol(ar1)+1-(vl:1)],type='l',xlab="Time (years)",ylab="AR(1) coefficient",ylim=c(min(B[1,]),max(B[2,])),
     xlim=c(1900,2025), main="", lwd=3, xaxt="n")
lines(x=jags_data$tdate, y=B[1,ncol(ar1)+1-(vl:1)],col=2,lty=2, lwd=3)
lines(x=jags_data$tdate, y=B[2,ncol(ar1)+1-(vl:1)],col=2,lty=2, lwd=3)
abline(h=1, lty=1, col="gray")
axis(1, at=c(seq(1900, 2020, by=5), 2024))
abline(h=0, lty=1, col="gray")
abline(v=c(1978.9, 1987.9, 2020, 2024.8), lty=1, col="gray")
#
ar1a <- apply(ar1,2,mean)[ncol(ar1)+1-(vl:1)]
ar1a[(tdate >= 1978.9) & (tdate < 1987.72) ] <- ar1a[(tdate >= 1978.9) & (tdate < 1987.72) ] + mean(ar10)
plot(x=jags_data$tdate, y=ar1a,type='l',xlab="Time (years)",ylab="AR(1) coefficient",ylim=c(min(B[1,]),max(B[2,])),
     xlim=c(1900,2025), main="", lwd=3, xaxt="n")
lines(x=jags_data$tdate, y=B[1,ncol(ar1)+1-(vl:1)],col=2,lty=2, lwd=3)
lines(x=jags_data$tdate, y=B[2,ncol(ar1)+1-(vl:1)],col=2,lty=2, lwd=3)
abline(h=1, lty=1, col="gray")
axis(1, at=c(seq(1900, 2020, by=5), 2024))
abline(h=0, lty=1, col="gray")
abline(v=c(1978.9, 1987.9, 2020, 2024.8), lty=1, col="gray")
#dev.off()
print(mean(ar10))

sigma <- get.posterior(pars="sigma[", coda_samples, jags_data, chain=1:3)
B=apply(sigma,2,quantile,c(0.05,0.95))
plot(x=jags_data$tdate, y=apply(sigma,2,mean),type='l',xlab="date",ylab="sigma",
     main="Standard deviation",  xlim=c(1978,2025))
lines(x=jags_data$tdate, y=B[1,],col=3,lty=2)
lines(x=jags_data$tdate, y=B[2,],col=3,lty=2)

#=======================================================================
# Plots
#=======================================================================

plot.trace.param.vec <- function(coda_samples, params, labels) {
    npar <- length(params)
    max_par <- 6
    if(npar > max_par){
        if((npar %% max_par) == 1) {
            max_par <- max_par - 1
        }
        par_chunks <- split(params, ceiling(seq_along(params)/max_par))
        lab_chunks <- split(labels, ceiling(seq_along(labels)/max_par))
        for(i in 1:length(par_chunks)){
            plot.trace.param.vec(coda_samples, par_chunks[[i]], lab_chunks[[i]])
        }
        return(NA)
    }
    posterior <- as.matrix(coda_samples[,params])
    plots <- list(NA)
    cur_plot=1
    for(i in 1:npar){
        pars <- c(params[i], params[i])
        plots[[cur_plot]] <- mcmc_trace(posterior[, c(i, i)], pars=pars)
            labs(y=labels[i])
        cur_plot<-cur_plot+1
        plots[[cur_plot]] <- mcmc_areas(posterior[, c(i,i)],
                                        pars=pars,
                                        prob = 0.8) +
            labs(y=labels[i]) + yaxis_text(on=FALSE)
        cur_plot<-cur_plot+1
    }
    if(npar < 6) {
        for(i in (npar+1):6) {
            for(j in 1:2){
                plots[[cur_plot]] <- grid.rect(gp=gpar(col="white"))
                cur_plot<-cur_plot+1
            }
        }
    }
    do.call("grid.arrange", c(plots, ncol=2))
}

plot.intervals.param.vec <- function(coda_samples, params) {
    posterior <- as.matrix(coda_samples[,params]) 
    mcmc_intervals(posterior)
}

plot.cin.errbar <- function(coda_samples, jags_data, chain=1, minmax=FALSE) {
    c_post <- get.cin.posterior(coda_samples, jags_data, chain)
    if(minmax) {
        probs <- c(0, .5, 1)
        title <- "Median with Min and Max values"
        df_q <- as.data.frame(colQuantiles(c_post, probs=probs)) %>%
            setNames(c("Lower", "Median", "Upper"))
    } else {
        sds <- colSds(c_post)
        meds <- colQuantiles(c_post, probs=c(0.5))
        df_q <- as.data.frame(cbind(meds-1.96*sds, meds, meds+1.96*sds)) %>%
            setNames(c("Lower", "Median", "Upper"))
        probs <- c(.05, .5, .95)
        title <- "Median +/- 1.96*std deviation"
    }
    ggplot(df_q, aes(x=row.names(df_q), y=Median)) +
        geom_point(size=2) +
        geom_errorbar(aes(ymin=Lower, ymax=Upper)) +
        scale_x_discrete(limits=c(1:dim(df_q)[1])) + 
        labs(x="Lag (Days)", y="Coefficient",
             title=title)
}
 
#-----------------------------------------------------------------------
## MCMC diagnostics
#-----------------------------------------------------------------------

pdf("../Figures/Nature-COMMSENV-s43247-023-01836-A-dynamic-ARIMA-mcmc_plots.pdf")

# Plot beta posterior
params <- c(sprintf("abeta1[%d]", c(1:beta_poly_order)), sprintf("mbeta1[%d]", c(1:1)), "sigma[527]")
params <- c(sprintf("abeta1[%d]", c(1:beta_poly_order)), "mbeta1", "sigma[527]", "ar2", "ar3", "ar10", "ar20", "ma1", "ma2","sigma_0")
labels <- params
plot.intervals.param.vec(coda_samples, params)
plot.trace.param.vec(coda_samples, params, labels)

dev.off()

pdf("../Figures/Nature-COMMSENV-s43247-023-01836-A-dynamic-ARIMA-loo.pdf")
# Model comparison
ll <- get.posterior(pars="loglik_", coda_samples, jags_data, chain=1:3)
gc <- rep(1:3,sapply(coda_samples, nrow))
r_eff <- loo::relative_eff(x = exp(ll), chain_id=gc)
loo.m <- loo::loo(ll, r_eff = r_eff, save_psis = FALSE)         
plot(loo.m)
loo.m
saveRDS(loo.m, file=paste0("../Data/Nature-COMMSENV-s43247-023-01836-A-dynamic-ARIMA-loo_",sim,".rds"))
} # end of foreach loop

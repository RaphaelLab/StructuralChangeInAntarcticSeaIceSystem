# This is the code for recreating the probability
# computations in Nature manuscript 2023-07-12845A
#
# "A Twenty-First Century Structural Change in Antarctica’s Sea ice System"
#  by Marilyn N. Raphael, Thomas J. Maierhofer, Ryan L. Fogt, William R. Hobbs, and Mark S. Handcock
#
# These standard R packages need to be present (and can be installed via 'install.packages("scales")', etc.
library(scales)
library(ggplot2)
library(dplyr)
library(mgcv)
library(abind)

#
# First we read in the NSIDC sea ice extent data 
# and the reconstructions of Maierhofer (2023) [See cite details in the manuscript]
load(file="../Data/nsidcV4.RData")
tempFile_location<- tempfile()
download.file(url="https://ucla.box.com/shared/static/zjiyg6pgqbx4xoiuyawvbnrniawi2c80.rds",
              destfile=tempFile_location, method="wget", mode="rb")
reconstructions <- readRDS(tempFile_location)
file.remove(tempFile_location)
total_fits = reconstructions[,,1]
rm(reconstructions)
ls()

set.seed(2)
#
pdf("../Figures/SuppFig3-details.pdf")
# ##############################################################################
# how unusual was the 2014-2017 decline in a historic context?
# First plot the recorded (total) SIE
# the largest anomaly was in January 2015 at 1.901878 the smallest in December 2016 at -2.20637
# check how often a jump of 4.108247 or more occured within 2 years
total_fits %>% str()
which(nsidc_data$Year==1979 & nsidc_data$Month == 12)
# Now look at the probability that over a 44 year window (like our satelite era)
# the annual min was less than
# that observed during the satellite era.
# It looks at 1/1899 to 12/1979
# and computes the min across all windows.
ann_min <- NULL
ann_min_when <- NULL
for(y in 1:(1979-1898-44)){
  ann_min = c(ann_min,apply(total_fits[,((12*y)-11):((y+43)*12)], 1, function(x) { min(x[seq(2,43*12+2,by=12)]) } ) )
  ann_min_when = c(ann_min_when,apply(total_fits[,((12*y)-11):((y+43)*12)], 1, function(x) { which.min(x[seq(2,43*12+2,by=12)]) } ) )
}
# year.m   sie.m
#   2017 2 2.238460
#   2022 2 2.222484
#   2023 2 2.044898
#   Year    total
#   2017 2 -0.9260761
#   2022 2 -0.9420516
#   2023 2 -1.1196377

off <- 3.164536
plot(density(ann_min+off), xlab = "annual min", main = "Distribution of annual min during each 44-year window",
     sub="The lines are the observed min for 2017, 2022 and 2023")
abline(v = -0.9260761+off, lty = 2, col=3, lwd=2)
abline(v = -0.9420516+off, lty = 2, col=4, lwd=2)
abline(v = -1.1196377+off, lty = 2, col=1, lwd=2)
legend("topleft", bty="n", lty=2, col=c(3,4,1), lwd=2, legend=c("2017", "2022","2023"))
# The probability of seeing such a drop at any time since 1/1899
mean(ann_min < -0.9260761)
mean(ann_min < -0.9420516)
mean(ann_min < -1.1196377)
#
# Now look at the probability that over a 44 year window (like our satelite era)
# the annual max was greater than
# that observed during the satellite era.
# It looks at 1/1899 to 12/1979
# and computes the max across all windows.
ann_max <- NULL
ann_max_when <- NULL
for(y in 1:(1979-1898-44)){
  ann_max = c(ann_max,apply(total_fits[,((12*y)-11):((y+43)*12)], 1, function(x) { max(x[seq(2,43*12+2,by=12)]) } ) )
  ann_max_when = c(ann_max_when,apply(total_fits[,((12*y)-11):((y+43)*12)], 1, function(x) { which.max(x[seq(2,43*12+2,by=12)]) } ) )
}
str(ann_max)
off <- 18.81309
plot(density(ann_max+off), xlab = "annual max", main = "Distribution of annual max during each 44-year window",
     sub="The lines are the observed max for (September) 2012 and 2014")
abline(v = 0.6450608+off, lty = 2, col=3, lwd=2)
abline(v = 1.2755211+off, lty = 2, col=1, lwd=2)
legend("topleft", bty="n", lty=2, col=c(3,1), lwd=2, legend=c("2012", "2014"))
# The probability of seeing such a drop at any time since 1/1899
mean(ann_max > 0.6450608)
mean(ann_max > 1.2755211)
#
table(ann_max_when)
a = apply(total_fits[,seq(2,81*12+2,by=12)], 1, function(x) { which.min(x) } )
ncol(total_fits)
81*12+2
table((1899:1980)[a])
#
y <- 1978:2023
ax.nophase.m <- rep((1:12),length(y)+1)
ax.nophase.m <- ax.nophase.m[1:length(sie.m)]
# Remove 1987.958 13.684670
sie.m[120]
sie.m[120] <- NA
# Remove 1987.958 13.684670
sie.m[115]
sie.m[115] <- NA

ac.fit.noamplitude.nophase <- gam(sie.m ~ s(ax.nophase.m,bs="cc"),knots=list(ax.nophase.m=c(0,12)),method="REML",subset=tdate.m >= 1979)
ay.noamplitude.nophase <- predict(ac.fit.noamplitude.nophase,newdata=data.frame(ax.nophase.m=ax.nophase.m, tdate.m=tdate.m))
anom=sie.m - ay.noamplitude.nophase
tail(cbind(seq(along=tdate.m),tdate.m, sie.m, anom))
actual_win_sie <- mean(sie.m[546:548])
actual_win_sie
actual_win_min <- mean(anom[546:548])
actual_win_min
# Now look at the prbability that over a 44 year window (like our satelite era)
# the annual winter season average was lower than
# that observed during the satellite era - that is winter 2023.
# It looks at 1/1899 to 12/1979
# and computes the min across all windows.
win_min <- NULL
win_min_when <- NULL
for(y in 1:(1979-1898-44)){
  win_min = c(win_min,apply(total_fits[,((12*y)-11):((y+43)*12)], 1, function(x) {
     xx <- rep(0,44)
     for(i in 1:44){
      xx[i] <-(x[6+12*(i-1)] + x[7+12*(i-1)] + x[8+12*(i-1)])/3
     }
     min(xx)
    }))
}
off <- actual_win_sie-actual_win_min
plot(density(win_min+off), xlab = "Winter average SIE", main = "Distribution of Winter average SIE during each 44-year window",
     sub="The line is the observed Winter SIE for 2023")
abline(v = actual_win_min+off, lty = 2, col=3, lwd=2)
# The probability of seeing such a drop at any time since 1/1899
mean(win_min < actual_win_min)
#
# So the 12/1979 was month 972
# and 1/1979 was month 961
# and 12/1978 was month 960
# Here each row is a (stochastic) reconstruction
# It looks at each 2 year (24 month) window
# from 1/1899 to 12/1979
# and computes the maximum drop across all windows.
max_jump_20th = apply(total_fits[,1:960], 1, function(x) {
  this.max = 0
  for (i in 1:(length(x)-24)) {
    this.window = min(length(x), i + 24) # 36
    if (abs(max(x[i:this.window]) - min(x[i:this.window])) > this.max) {
      this.max = abs(max(x[i:this.window],na.rm=TRUE) - min(x[i:this.window],na.rm=TRUE))
    }
  }
  return(this.max)
})
#  Compute the actual maximum jump
  actual_max_jump = 0
  ncol(total_fits)-8
  colnames(total_fits)[c(961,1456, 1464)]
  x <- total_fits[1,961:1464]
  x[108] <- NA # Dec 1987
  for (i in 1:(length(x)-24)) {
    this.window = min(length(x), i + 24) # 36
    if (abs(max(x[i:this.window],na.rm=TRUE) - min(x[i:this.window],na.rm=TRUE)) > actual_max_jump) {
      actual_max_jump = abs(max(x[i:this.window],na.rm=TRUE) - min(x[i:this.window],na.rm=TRUE))
    }
  }
which.min(x)
x[456]
which.max(x)
x[433]
max(x,na.rm=T)
min(x,na.rm=T)
plot(x=seq(1979+0.5/12, 2020+11.5/12, by=1/12),y=x)
abline(v=c(2015+0.5/12, 2016+11.5/12))
actual_max_jump
plot(density(max_jump_20th), xlim = c(0.5, 5), xlab = "Two-year maximal difference", main = "Distribution of max drops during each reconstruction",
     sub="The line in the observed drop (1/2015 to 12/2016)")
abline(v = actual_max_jump, lty = 2)
# The probability of seeing such a drop at any time since 1/1899
mean(max_jump_20th > actual_max_jump)
# The number of reconstructions seeing such a drop at any time since 1/1899
sum(max_jump_20th > actual_max_jump)

# Now look at the probability that over a 44 year window (like our satelite era)
# That is the proportion of 44 year windows where a 2 year drop was greater than
# that observed during hte satellite era.
# It looks at each 2 year (24 month) drop in a 44 year window from 
# 1/1899 to 12/1979
# and computes the maximum drop across all windows.
max_jump_44 <- NULL
for(y in 1:(1979-1898-44)){
max_jump_44 = c(max_jump_44,apply(total_fits[,((12*y)-11):((y+43)*12)], 1, function(x) {
  this.max = 0
  for (i in 1:(length(x)-24)) {
    this.window = min(length(x), i + 24) # 24
    if (abs(max(x[i:this.window]) - min(x[i:this.window])) > this.max) {
      this.max = abs(max(x[i:this.window]) - min(x[i:this.window]))
    }
  }
  return(this.max)
}))
}
str(max_jump_44)
plot(density(max_jump_44), xlim = c(0.5, 5), xlab = "Two-year maximal difference", main = "Distribution of max drops during each 44 year period
over the reconstruction",
     sub="The line in the observed drop (1/2015 to 12/2016)")
abline(v = actual_max_jump, lty = 2)
# The probability of seeing such a drop at any time since 1/1899
mean(max_jump_44 > actual_max_jump)
# The number of reconstructions seeing such a drop at any time since 1/1899
sum(max_jump_44 > actual_max_jump)

# Now all three mins occuring within a 6-year window
min_three <- NULL
for(y in 1:(1979-1898-44)){
min_three = c(min_three,apply(total_fits[,((12*y)-11):((y+43)*12)], 1, function(x) {
  oall <- rank(x)
  w <- rep(0, length(x))
  for (i in 1:(length(x)-(6*12))) {
    this.window = min(length(x), i + (6*12))
    w[i] <- all(c(1,2,3) %in% oall[i:this.window])
  }
  return(w)
}))
}
# The probability of seeing such a drop at any time since 1/1899
mean(min_three)
# The number of reconstructions seeing such a drop at any time since 1/1899
sum(min_three)
dev.off()
#
pdf(file = "../Figures/SuppFig3.pdf", width = 8, height = 8)
par(mfrow=c(2,3))
plot(density(max_jump_20th), xlim = c(0.5, 5), xlab = "Two-year maximal difference", main = "(a)", #Distribution of max drops during each reconstruction",
     cex.sub=0.8, sub="The line in the max recorded drop (1/2015 to 12/2016)")
abline(v = actual_max_jump, lty = 2)
plot(density(max_jump_44), xlim = c(0.5, 5), xlab = "Two-year maximal difference", main = "(b)",#Distribution of max drops during each 44 year period over the reconstruction",
     cex.sub=0.8, sub="The line in the max recorded drop (1/2015 to 12/2016)")
abline(v = actual_max_jump, lty = 2)
off <- 3.164536
plot(density(ann_min+off), xlab = "Annual min", main = "(c)",#Distribution of annual min during each 44-year window",
     cex.sub=0.7, sub="The lines are the recorded min for 2017, 2022 and 2023")
abline(v = -0.9260761+off, lty = 2, col=3, lwd=2)
abline(v = -0.9420516+off, lty = 2, col=4, lwd=2)
abline(v = -1.1196377+off, lty = 2, col=1, lwd=2)
legend("topleft", bty="n", lty=2, col=c(3,4,1), lwd=2, legend=c("2017", "2022","2023"))
off <- 18.81309
plot(density(ann_max+off), xlab = "Annual max", main = "(d)", #Distribution of annual max during each 44-year window",
     cex.sub=0.7,sub="The lines are the recorded max for (September) 2012 and 2014")
abline(v = 0.6450608+off, lty = 2, col=3, lwd=2)
abline(v = 1.2755211+off, lty = 2, col=1, lwd=2)
legend("topright", bty="n", lty=2, col=c(3,1), lwd=2, legend=c("2012", "2014"))
off <- actual_win_sie-actual_win_min
plot(density(win_min+off), xlab = "Winter average SIE", main = "(e)", #Distribution of Winter average SIE during each 44-year window",
     cex.sub=0.7, sub="The line is the observed Winter SIE for 2023")
abline(v = actual_win_min+off, lty = 2, col=3, lwd=2)
dev.off()

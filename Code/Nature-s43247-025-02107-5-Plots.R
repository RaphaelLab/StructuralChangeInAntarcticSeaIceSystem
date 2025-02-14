# This is the code for recreating the plots in
# Nature manuscript 2023-07-12845A
#
# "A Twenty-First Century Structural Change in Antarctica’s Sea ice System"
#  by Marilyn N. Raphael, Thomas J. Maierhofer, Ryan L. Fogt, William R. Hobbs, and Mark S. Handcock
#
# These standard R packages need to be present (and can be installed via
# 'install.packages("scales")', etc.
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
# These are up to 2020 but include the sectors and total.
# The observed are augmented later
reconstructions <- readRDS(tempFile_location)
dim(reconstructions)
file.remove(tempFile_location)
fits = reconstructions[,,2:6]
total_fits = reconstructions[,,1]
rm(reconstructions)
#
tempFile_location<- tempfile()
download.file(url="https://ucla.box.com/shared/static/az3wdc2celjjmjgr25exix5u76myrw96",
              destfile=tempFile_location, method="wget", mode="rb")
# These are the posterior predictive draws
load(file=tempFile_location)
file.remove(tempFile_location)

set.seed(2)

#
pdf(file = "../Figures/Fig1.pdf", width = 8, height = 8)
par(mfrow=c(3,2))
# c(bottom, left, top, right)
par(mar=c(0,4.2,0,0))
# Total
matplot(nsidc_data$tdate[1:957], t(total_fits[c(T,rep(F,5)),1:957]), type = "l", lty = 1, col = alpha(1, 0.01),
        xlab = "Time (years)", ylab = "",
        ylim = c(-2, 2), axes=F,
        xlim = c(1899, 2021))
lines(nsidc_data$tdate[1:957], colMeans(total_fits)[1:957], col = 2)
lines(nsidc_data$tdate, nsidc_data$total, col = "steelblue")
axis(2)
text(x=2020, y=2.05, labels="Total", col=2)
abline(h=0, col="gray")
#
# King_Hakon
par(mar=c(0,0,0,2))
matplot(nsidc_data$tdate[1:957], t(fits[c(T,rep(F,5)),1:957,1]), type = "l", lty = 1, col = alpha(1, 0.01),
        xlab = "Time (years)", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'),
        ylim = c(-1.5, 1.5), axes=FALSE,
        xlim = c(1899, 2021))
lines(nsidc_data$tdate[1:957], colMeans(fits[,1:957,1]), col = 2)
lines(nsidc_data$tdate, nsidc_data$King_Hakon, col = "steelblue")
axis(4)
text(x=2007, y=1.5, labels="King Haakon VII", col=2)
abline(h=0, col="gray")

# Ross
# c(bottom, left, top, right)
par(mar=c(0,4.2,0,0))
matplot(nsidc_data$tdate[1:957], t(fits[c(T,rep(F,5)),1:957,2]), type = "l", lty = 1, col = alpha(1, 0.01),
        xlab = "Time (years)", ylab = "Sea ice anomaly in millions of sq. km", 
        ylim = c(-1.5, 1.5), axes=FALSE,
        xlim = c(1899, 2021))
lines(nsidc_data$tdate[1:957], colMeans(fits[,1:957,2]), col = 2)
lines(nsidc_data$tdate, nsidc_data$Ross, col = "steelblue")
axis(2)
text(x=2007, y=1.4, labels="Ross-Amundsen Sea", col=2)
abline(h=0, col="gray")

# East_Antarctica
par(mar=c(0,0,0,2))
matplot(nsidc_data$tdate[1:957], t(fits[c(T,rep(F,5)),1:957,3]), type = "l", lty = 1, col = alpha(1, 0.01),
        xlab = "Time (years)", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'),
        ylim = c(-0.75, 0.75), axes=F,
        xlim = c(1899, 2021))
lines(nsidc_data$tdate[1:957], colMeans(fits[,1:957,3]), col = 2)
lines(nsidc_data$tdate, nsidc_data$East_Antarctica, col = "steelblue")
axis(4)
text(x=2007, y=0.7, labels="East Antarctica", col=2)
abline(h=0, col="gray")

# Weddell
par(mar=c(4.2,4.2,0,0))
matplot(nsidc_data$tdate[1:957], t(fits[c(T,rep(F,5)),1:957,4]), type = "l", lty = 1, col = alpha(1, 0.01),
        xlab = "Time (years)", ylab = "",
        ylim = c(-1.5, 1.5), axes=F,
        xlim = c(1899, 2021))
lines(nsidc_data$tdate[1:957], colMeans(fits[,1:957,4]), col = 2)
lines(nsidc_data$tdate, nsidc_data$Weddell, col = "steelblue")
text(x=2007, y=1.5, labels="Weddell", col=2)
axis(1)
axis(2, at=c(-1.5, -1,-0.5,0,0.5,1,1.5), labels=c("-1.5", "-1","-0.5","0","0.5","1",""))
abline(h=0, col="gray")

# Bellingshausen_Amundsen_Sea
par(mar=c(4.2,0,0,2))
matplot(nsidc_data$tdate[1:957], t(fits[c(T,rep(F,5)),1:957,5]), type = "l", lty = 1, col = alpha(1, 0.01),
        xlab = "Time (years)", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'),
        ylim = c(-0.75, 0.75), axes=F,
        xlim = c(1899, 2021))
lines(nsidc_data$tdate[1:957], colMeans(fits[,1:957,5]), col = 2)
lines(nsidc_data$tdate, nsidc_data$Bellingshausen_Amundsen_Sea, col = "steelblue")
axis(1)
axis(4)
text(x=2004.5, y=0.765, labels="Bellinghausen Amundsen Sea", col=2, cex=0.81)
abline(h=0, col="gray")
dev.off()

# ##############################################################################
# analyze the predicted annual cycle over the course of the 20th century
# compute average prediction by sector and month
avg_prediction = data.frame("King_Haakon" = NA, "Ross" = NA,
                            "East_Antarctica" = NA, "Weddell" = NA,
                            "Bellingshausen_Amundsen" = NA)
for (i in 1:dim(fits)[2]) {
  avg_prediction[i, 1:5] = colMeans(fits[,i,])
}
avg_prediction$Total = rowSums(avg_prediction)

# compute sds
for (i in 1:dim(fits)[2]) {
  avg_prediction[i, 7:11] = apply(fits[,i,], 2, sd)
}
avg_prediction[, 12] = apply(total_fits[], 2, sd)
names(avg_prediction)[7:12] = c("King_Haakon_sd", "Ross_sd", "East_Antarctica_sd", "Weddell_sd",
                                "Bellingshausen_Amundsen_sd", "Total_sd")

# get a prediction for the total
avg_prediction$Month = 1:12
avg_prediction$Season = c(1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 1, 1)
avg_prediction$Year = rep(1899:2020, each = 12)
avg_prediction$Decade = as.character(rep(seq(1900, 2020, by = 10), each = 12*10))[1:1464]

# compute monthly average per decade
decadal_avg = avg_prediction %>%
  group_by(Decade, Month) %>%
  summarize_all(.funs = mean)

decadal_avg$King_Haakon_text = paste0(weights::rd(decadal_avg$King_Haakon, digits=2, add=FALSE), "(",
                                      weights::rd(decadal_avg$King_Haakon_sd, digits=2, add=FALSE), ")")
decadal_avg$Ross_text = paste0(weights::rd(decadal_avg$Ross, digits=2, add=FALSE), "(",
                               weights::rd(decadal_avg$Ross_sd, digits=2, add=FALSE), ")")
decadal_avg$East_Antarctica_text = paste0(weights::rd(decadal_avg$East_Antarctica, digits=2, add=FALSE), "(",
                                          weights::rd(decadal_avg$East_Antarctica_sd, digits=2, add=FALSE), ")")
decadal_avg$Weddell_text = paste0(weights::rd(decadal_avg$Weddell, digits=2, add=FALSE), "(",
                                  weights::rd(decadal_avg$Weddell_sd, digits=2, add=FALSE), ")")
decadal_avg$Bellingshausen_Amundsen_text = paste0(weights::rd(decadal_avg$Bellingshausen_Amundsen, digits=2, add=FALSE), "(",
                                                  weights::rd(decadal_avg$Bellingshausen_Amundsen_sd, digits=2, add=FALSE), ")")
decadal_avg$Total_text = paste0(weights::rd(decadal_avg$Total, digits=2, add=FALSE), "(",
                                weights::rd(decadal_avg$Total_sd, digits=2, add=FALSE), ")")

pdf(file = "../Figures/Fig2.pdf", width = 8, height = 8)

fit_total <- gam(Total ~ s(Month, Year), data=avg_prediction, method="ML")
summary(fit_total)

grd <- rep(1899:2020, each = 12/0.2)
grdM <- seq(0,12-0.2,by=0.2)
grd <- data.frame( Month=rep(grdM, length(grd)/length(grdM)), Year=grd )
p_cplot = data.frame("Month"=grd$Month, "Year"=grd$Year,
            "Total"= grd$Month,
            "King_Haakon" = grd$Month,
            "Ross" = grd$Month,
            "East_Antarctica" = grd$Month,
            "Weddell" = grd$Month,
            "Bellingshausen_Amundsen" = grd$Month )

p_cplot[,"Total"]=predict(gam(Total ~ te(Month, Year, bs = c("cc", "cs"), k=c(12,20)), data=avg_prediction, method="ML"), newdata=grd)
p_cplot[,"King_Haakon"]=predict(gam(King_Haakon ~ te(Month, Year, bs = c("cc", "cs"), k=c(12,20)), data=avg_prediction, method="ML"), newdata=grd)
p_cplot[,"Ross"]=predict(gam(Ross ~ te(Month, Year, bs = c("cc", "cs"), k=c(12,20)), data=avg_prediction, method="ML"),
newdata=grd)
p_cplot[,"East_Antarctica"]=predict(gam(East_Antarctica ~ te(Month, Year, bs = c("cc", "cs"), k=c(12,20)), data=avg_prediction, method="ML"), newdata=grd)
p_cplot[,"Weddell"]=predict(gam(Weddell ~ te(Month, Year, bs = c("cc", "cs"), k=c(12,20)), data=avg_prediction, method="ML"), newdata=grd)
p_cplot[,"Bellingshausen_Amundsen"]=predict(gam(Bellingshausen_Amundsen ~ te(Month, Year, bs = c("cc", "cs"), k=c(12,20)), data=avg_prediction, method="ML"), newdata=grd)

p_stack <- data.frame(Month=rep(p_cplot$Month, 6), Year=rep(p_cplot$Year,6),
  Sector=rep(c("Total","King Haakon VII","Ross-Amundsen","East Antarctica","Weddell","Bellingshausen Amundsen"),each=length(p_cplot$Total)),
  Anomaly=c(p_cplot[,"Total"],p_cplot[,"King_Haakon"],p_cplot[,"Ross"],p_cplot[,"East_Antarctica"],
            p_cplot[,"Weddell"],p_cplot[,"Bellingshausen_Amundsen"]) )

# Heat map colors
mybreaks <- seq(min(p_stack$Anomaly), max(p_stack$Anomaly), length.out = 13)

mycolors<- function(x) {
   colors<-colorRampPalette(c("darkblue", "yellow"))( 12 )
   colors[1:x]
}
breaklabel <- function(x){
   labels<- paste0(mybreaks[1:12], "-", mybreaks[2:13])
   labels<- paste0(mybreaks[1:12], "-", mybreaks[2:13])
   labels[1:x]
}
cols=c("purple","darkblue","blue","dodgerblue3","dodgerblue","deepskyblue2","deepskyblue","skyblue","lightskyblue2","lightskyblue3","lightgoldenrod1","khaki1","yellow","goldenrod1","orange","darkorange","firebrick2","red","red2","brown2","indianred","red3")
cols <- rev(cols)
mybreaks <- c(seq(-0.85, 0, length.out = 10), seq(0, 0.6, length=12))
mycolors<- function(x) { cols }
base_size <- 5.5
p <- ggplot(p_stack, aes(Year, Month, z=Anomaly)) + theme_minimal() +
  theme(plot.margin = margin(base_size / 2, 10*base_size, base_size / 2, base_size / 2, unit = "pt")) +
  geom_contour_filled(breaks=mybreaks, show.legend = FALSE) +
  scale_fill_manual(name="Anomaly", values=mycolors(1), drop=FALSE) #, breaks=seq(-0.85,0.6,length=length(cols)+2)) +
p <- p + scale_y_continuous(breaks=1:12, label=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
p <- p + scale_x_continuous(breaks=seq(1900,2020, by=10), label=c("1900","","1920","","1940","","1960","","1980","","2000","","2020"))
print(p + facet_wrap(~factor(Sector, levels=c("Total","King Haakon VII","Ross-Amundsen","East Antarctica","Weddell","Bellingshausen Amundsen")), ncol=2))
dev.off()

# plot individual draws from posterior distribution
pdf(file = "../Figures/SuppFig1.pdf", width = 8, height = 8)
par(mfrow=c(3,2))
# c(bottom, left, top, right)
par(mar=c(0,4.2,0,0))

plot(nsidc_data$tdate[1:957], t(total_fits[c(100),1:957]), type = "l", lty = 1, col = c("steelblue"),
     xlim = c(1899, 2021), ylim = c(-2.2, 2.2),
     axes=F,
     xlab = "Time (year)")
axis(2)
lines(nsidc_data$tdate, nsidc_data$total, col = "steelblue")
abline(v = 1978.792, lty = 2)

par(mar=c(0,0,0,2))
plot(nsidc_data$tdate[1:957], t(total_fits[c(500),1:957]), type = "l", lty = 1, col = c("steelblue"),
     xlim = c(1899, 2021), ylim = c(-2.2, 2.2),
     axes=F,
     xlab = "Time (year)")
axis(4)
lines(nsidc_data$tdate, nsidc_data$total, col = "steelblue")
abline(v = 1978.792, lty = 2)

par(mar=c(0,4.2,0,0))
plot(nsidc_data$tdate[1:957], t(total_fits[c(900),1:957]), type = "l", lty = 1, col = c("steelblue"),
     xlim = c(1899, 2021), ylim = c(-2.2, 2.2),
     axes=F,
     main = "", ylab = "Sea ice anomaly in millions of sq. km", xlab = "Time (year)")
lines(nsidc_data$tdate, nsidc_data$total, col = "steelblue")
abline(v = 1978.792, lty = 2)
axis(2)

par(mar=c(0,0,0,2))
plot(nsidc_data$tdate[1:957], t(total_fits[c(1200),1:957]), type = "l", lty = 1, col = c("steelblue"),
     xlim = c(1899, 2021), ylim = c(-2.2, 2.2),
     axes=F,
     xlab = "Time (year)")
lines(nsidc_data$tdate, nsidc_data$total, col = "steelblue")
abline(v = 1978.792, lty = 2)
axis(4)

par(mar=c(4.2,4.2,0,0))
plot(nsidc_data$tdate[1:957], t(total_fits[c(1600),1:957]), type = "l", lty = 1, col = c("steelblue"),
     xlim = c(1899, 2021), ylim = c(-2.2, 2.2),
     axes=F,
     xlab = "Time (year)")
lines(nsidc_data$tdate, nsidc_data$total, col = "steelblue")
abline(v = 1978.792, lty = 2)
axis(1)
axis(2)

par(mar=c(4.2,0,0,2))
plot(nsidc_data$tdate[1:957], t(total_fits[c(2000),1:957]), type = "l", lty = 1, col = c("steelblue"),
     xlim = c(1899, 2021), ylim = c(-2.2, 2.2),
     axes=F,
     xlab = "Time (year)")
lines(nsidc_data$tdate, nsidc_data$total, col = "steelblue")
abline(v = 1978.792, lty = 2)
axis(1)
axis(4)
dev.off()

pdf(file = "../Figures/Nature-s43247-025-02107-5-PSIS.pdf", width = 8, height = 8)

# These are the components of Supplementary Figure 2
# King Haakon
ggplot(decadal_avg) +
  geom_raster(aes(x = Decade, y = Month, fill = King_Haakon), show.legend = FALSE) +
  # geom_text(aes(x = Decade, y = Month, label = round(King_Haakon, 2))) +
  geom_text(aes(x = Decade, y = Month, label = King_Haakon_text)) +
  scale_fill_gradient2("Anomaly",
                       limits = c(-0.2, 0.2), oob = scales::squish,
                       low = "firebrick",
                       mid = "white",
                       high = "steelblue") +
  scale_y_continuous(breaks = 1:12) +
  ggtitle("King Haakon VII") +
  theme_minimal()
ggsave("../Figures/decadal_avg_king_haakon_R2.png", device = "png", height = 4, width = 7)

# Ross
ggplot(decadal_avg) +
  geom_raster(aes(x = Decade, y = Month, fill = Ross), show.legend = FALSE) +
  geom_text(aes(x = Decade, y = Month, label = Ross_text)) +
  scale_fill_gradient2(limits = c(-0.2, 0.2), oob = scales::squish,
                       low = "firebrick",
                       mid = "white",
                       high = "steelblue") +
  scale_y_continuous(breaks = 1:12) +
  ggtitle("Ross-Amundsen Sea") +
  theme_minimal()
ggsave("../Figures/decadal_avg_ross_R2.png", device = "png", height = 4, width = 7)

# East Antarctica
ggplot(decadal_avg) +
  geom_raster(aes(x = Decade, y = Month, fill = East_Antarctica), show.legend = FALSE) +
  geom_text(aes(x = Decade, y = Month, label = East_Antarctica_text)) +
  scale_fill_gradient2("Anomaly",
                       limits = c(-0.2, 0.2), oob = scales::squish,
                       low = "firebrick",
                       mid = "white",
                       high = "steelblue") +
  ggtitle("East Antarctica") +
  scale_y_continuous(breaks = 1:12) +
  theme_minimal()
ggsave("../Figures/decadal_avg_east_antarctica_R2.png", device = "png", height = 4, width = 7)

# Weddell Sea
ggplot(decadal_avg) +
  geom_raster(aes(x = Decade, y = Month, fill = Weddell), show.legend = FALSE) +
  geom_text(aes(x = Decade, y = Month, label = Weddell_text)) +
  scale_fill_gradient2("Anomaly",
                       limits = c(-0.4, 0.4), oob = scales::squish,
                       low = "firebrick",
                       mid = "white",
                       high = "steelblue") +
  ggtitle("Weddell Sea") +
  scale_y_continuous(breaks = 1:12) +
  theme_minimal()
ggsave("../Figures/decadal_avg_weddell_R2.png", device = "png", height = 4, width = 7)

# Bellingshausen Amundsen
ggplot(decadal_avg) +
  geom_raster(aes(x = Decade, y = Month, fill = Bellingshausen_Amundsen), show.legend = FALSE) +
  geom_text(aes(x = Decade, y = Month, label = Bellingshausen_Amundsen_text)) +
  scale_fill_gradient2("Anomaly",
                       limits = c(-0.2, 0.2), oob = scales::squish,
                       low = "firebrick",
                       mid = "white",
                       high = "steelblue") +
  scale_y_continuous(breaks = 1:12) +
  ggtitle("Bellingshausen Amundsen Sea") +
  theme_minimal()
ggsave("../Figures/decadal_avg_bellingshausen_amundsen_R2.png", device = "png", height = 4, width = 7)

# Total
ggplot(decadal_avg) +
  geom_raster(aes(x = Decade, y = Month, fill = Total), show.legend = FALSE) +
  geom_text(aes(x = Decade, y = Month, label = Total_text)) +
  scale_fill_gradient2("Anomaly",
                       limits = c(-0.4, 0.4), oob = scales::squish,
                       low = "firebrick",
                       mid = "white",
                       high = "steelblue") +
  scale_y_continuous(breaks = 1:12) +
  ggtitle("Total") +
  theme_minimal()
ggsave("../Figures/decadal_avg_total_R2.png", device = "png", height = 4, width = 7)

#load(file="../Data/pos_pred_fits.RData")

pdf(file = "../Figures/Fig1_Supplemental_1.pdf", width = 8, height = 8)

start <- 959
end <- 1500

par(mfrow=c(3,2))
# c(bottom, left, top, right)
par(mar=c(0,4.2,0,0))
# Total
matplot(all_dat$tdate[start:end], t(total_PP[c(T,rep(T,5)),start:end]), type = "l", lty = 1, col = alpha(1, 0.01),
        xlab = "Time (years)", ylab = "",
        ylim = c(-3, 3), axes=F,
        xlim = c(1978, 2024))
lines(all_dat$tdate[start:end], colMeans(total_fits)[start:end], col = 2)
lines(all_dat$tdate[start:end], all_dat$total[start:end], col = "steelblue", lwd=2)
axis(2)
text(x=2020, y=3.0, labels="Total", col=2)
abline(h=0, col="gray")
#
# King_Hakon
par(mar=c(0,0,0,2))
matplot(all_dat$tdate[start:end], t(sector_PP[c(T,rep(T,5)),start:end,1]), type = "l", lty = 1, col = alpha(1, 0.01),
        xlab = "Time (years)", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'),
        ylim = c(-1.8, 1.8), axes=FALSE,
        xlim = c(1978, 2024))
lines(all_dat$tdate[start:end], colMeans(sector_PP[,start:end,1]), col = 2)
lines(all_dat$tdate[start:end], all_dat$King_Hakon[start:end], col = "steelblue", lwd=2)
axis(4)
text(x=2007, y=1.8, labels="King Haakon VII", col=2)
abline(h=0, col="gray")

# Ross
# c(bottom, left, top, right)
par(mar=c(0,4.2,0,0))
matplot(all_dat$tdate[start:end], t(sector_PP[c(T,rep(F,5)),start:end,2]), type = "l", lty = 1, col = alpha(1, 0.01),
        xlab = "Time (years)", ylab = "Sea ice anomaly in millions of sq. km", 
        ylim = c(-1.2, 1.2), axes=FALSE,
        xlim = c(1978, 2024))
lines(all_dat$tdate[start:end], colMeans(sector_PP[,start:end,2]), col = 2)
lines(all_dat$tdate[start:end], all_dat$Ross[start:end], col = "steelblue", lwd=2)
axis(2)
text(x=1987, y=1.0, labels="Ross-Amundsen Sea", col=2)
abline(h=0, col="gray")

# East_Antarctica
par(mar=c(0,0,0,2))
matplot(all_dat$tdate[start:end], t(sector_PP[c(T,rep(F,5)),start:end,3]), type = "l", lty = 1, col = alpha(1, 0.01),
        xlab = "Time (years)", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'),
        ylim = c(-0.8, 0.8), axes=F,
        xlim = c(1978, 2024))
lines(all_dat$tdate[start:end], colMeans(sector_PP[,start:end,3]), col = 2)
lines(all_dat$tdate[start:end], all_dat$East_Antarctica[start:end], col = "steelblue", lwd=2)
axis(4)
text(x=2005, y=0.7, labels="East Antarctica", col=2)
abline(h=0, col="gray")

# Weddell
par(mar=c(4.2,4.2,0,0))
matplot(all_dat$tdate[start:end], t(sector_PP[c(T,rep(F,5)),start:end,4]), type = "l", lty = 1, col = alpha(1, 0.01),
        xlab = "Time (years)", ylab = "",
        ylim = c(-1.0, 1.0), axes=F,
        xlim = c(1978, 2024))
lines(all_dat$tdate[start:end], colMeans(sector_PP[,start:end,4]), col = 2)
lines(all_dat$tdate[start:end], all_dat$Weddell[start:end], col = "steelblue", lwd=2)
text(x=2007, y=-1.0, labels="Weddell", col=2)
axis(1)
axis(2, at=c(-1.5, -1,-0.5,0,0.5,1,1.5), labels=c("-1.5", "-1","-0.5","0","0.5","1",""))
abline(h=0, col="gray")

# Bellingshausen_Amundsen_Sea
par(mar=c(4.2,0,0,2))
matplot(all_dat$tdate[start:end], t(sector_PP[c(T,rep(F,5)),start:end,5]), type = "l", lty = 1, col = alpha(1, 0.01),
        xlab = "Time (years)", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'),
        ylim = c(-0.6, 0.6), axes=F,
        xlim = c(1978, 2024))
lines(all_dat$tdate[start:end], colMeans(sector_PP[,start:end,5]), col = 2)
lines(all_dat$tdate[start:end], all_dat$Bellingshausen_Amundsen_Sea[start:end], col = "steelblue", lwd=2)
axis(1)
axis(4)
text(x=1989.0, y=-0.6, labels="Bellinghausen Amundsen Sea", col=2, cex=0.81)
abline(h=0, col="gray")
dev.off()

#load(file="pos_pred_total.RData")
#total_fits <- r_total
#load(file="pos_pred.RData")
#sector_fits <- r
load(file="../Data/pos_pred_fits_brief.RData")

# before 2015 and 2015 onward
end0 <- end - 12*9
start1 <- end0 + 1
end01 <- end0 - start + 1
start11 <- start1 - start + 1
end1 <- end - start + 1

pdf(file = "../Figures/Fig1_Supplemental_2.pdf", width = 8, height = 8)
par(mfrow=c(3,2))
# c(bottom, left, top, right)
par(mar=c(0,4.2,0,0))
# Total
matplot(all_dat$tdate[start:end], 100*total_PP_brief[start:end], type = "l", lty = 1, col = alpha(1, 0.01),
        xlab = "Time (years)", ylab = "",
        ylim = c(0, 100), axes=F,
        xlim = c(1978, 2024))
lines(all_dat$tdate[start:end0], 100*total_PP_brief[1:end01], col = "steelblue", lwd=2)
lines(all_dat$tdate[start1:end], 100*total_PP_brief[start11:end1], col = "darkgreen", lwd=2)
axis(2)
text(x=2020, y=100*1.0, labels="Total", col=2)
abline(h=0, col="gray")
abline(h=100*c(0.05, 0.95, 0), lty=2, col=2)
#
# King_Hakon
par(mar=c(0,0,0,2))
matplot(all_dat$tdate[start:end], 100*sector_PP_brief[,1], type = "l", lty = 1, col = alpha(1, 0.01),
        xlab = "Time (years)", ylab = parse(text='"Percentile Rank of Sea ice anomaly"'),
        ylim = 100*c(0,1), axes=FALSE,
        xlim = c(1978, 2024))
lines(all_dat$tdate[start:end0], 100*sector_PP_brief[1:end01,1], col = "steelblue", lwd=2)
lines(all_dat$tdate[start1:end], 100*sector_PP_brief[start11:end1,1], col = "darkgreen", lwd=2)
axis(4)
text(x=2003, y=100*1.0, labels="King Haakon VII", col=2)
abline(h=0, col="gray")
abline(h=100*c(0.05, 0.95, 0), lty=2, col=2)

# Ross
# c(bottom, left, top, right)
par(mar=c(0,4.2,0,0))
matplot(all_dat$tdate[start:end], 100*sector_PP_brief[,2], type = "l", lty = 1, col = alpha(1, 0.01),
        xlab = "Time (years)", ylab = "Percentile Rank of Sea ice anomaly",
        ylim = 100*c(0,1), axes=FALSE,
        xlim = c(1978, 2024))
lines(all_dat$tdate[start:end0], 100*sector_PP_brief[1:end01,2], col = "steelblue", lwd=2)
lines(all_dat$tdate[start1:end], 100*sector_PP_brief[start11:end1,2], col = "darkgreen", lwd=2)
axis(2)
text(x=2003, y=100*1.0, labels="Ross-Amundsen Sea", col=2)
abline(h=0, col="gray")
abline(h=100*c(0.05, 0.95, 0), lty=2, col=2)

# East_Antarctica
par(mar=c(0,0,0,2))
matplot(all_dat$tdate[start:end], 100*sector_PP_brief[,3], type = "l", lty = 1, col = alpha(1, 0.01),
        xlab = "Time (years)", ylab = parse(text='"Percentile Rank of Sea ice anomaly"'),
        ylim = 100*c(0,1), axes=FALSE,
        xlim = c(1978, 2024))
lines(all_dat$tdate[start:end0], 100*sector_PP_brief[1:end01,3], col = "steelblue", lwd=2)
lines(all_dat$tdate[start1:end], 100*sector_PP_brief[start11:end1,3], col = "darkgreen", lwd=2)
axis(4)
text(x=2001, y=100*1.0, labels="East Antarctica", col=2)
abline(h=0, col="gray")
abline(h=100*c(0.05, 0.95, 0), lty=2, col=2)

# Weddell
par(mar=c(4.2,4.2,0,0))
matplot(all_dat$tdate[start:end], 100*sector_PP_brief[,4], type = "l", lty = 1, col = alpha(1, 0.01),
        xlab = "Time (years)", ylab = "",
        ylim = 100*c(0,1), axes=FALSE,
        xlim = c(1978, 2024))
lines(all_dat$tdate[start:end0], 100*sector_PP_brief[1:end01,4], col = "steelblue", lwd=2)
lines(all_dat$tdate[start1:end], 100*sector_PP_brief[start11:end1,4], col = "darkgreen", lwd=2)
text(x=1985, y=100*1.0, labels="Weddell", col=2)
axis(2)
axis(1)
abline(h=0, col="gray")
abline(h=100*c(0.05, 0.95, 0), lty=2, col=2)

# Bellingshausen_Amundsen_Sea
par(mar=c(4.2,0,0,2))
matplot(all_dat$tdate[start:end], 100*sector_PP_brief[,5], type = "l", lty = 1, col = alpha(1, 0.01),
        xlab = "Time (years)", ylab = parse(text='"Percentile Rank of Sea ice anomaly"'),
        ylim = 100*c(0,1), axes=FALSE,
        xlim = c(1978, 2024))
lines(all_dat$tdate[start:end0], 100*sector_PP_brief[1:end01,5], col = "steelblue", lwd=2)
lines(all_dat$tdate[start1:end], 100*sector_PP_brief[start11:end1,5], col = "darkgreen", lwd=2)
axis(1)
axis(4)
text(x=2004.5, y=100*1.0, labels="Bellinghausen Amundsen Sea", col=2, cex=0.81)
abline(h=0, col="gray")
abline(h=100*c(0.05, 0.95, 0), lty=2, col=2)
dev.off()

avg_prediction = data.frame("King_Haakon" = NA, "Ross" = NA,
                            "East_Antarctica" = NA, "Weddell" = NA,
                            "Bellingshausen_Amundsen" = NA)
for (i in 1:dim(fits)[2]) {
  avg_prediction[i, 1:5] = colMeans(fits[,i,])
}
avg_prediction$Total = rowSums(avg_prediction)

# think about these sds
# compute sds
for (i in 1:dim(fits)[2]) {
  avg_prediction[i, 7:11] = apply(fits[,i,], 2, sd)
}
avg_prediction[, 12] = apply(total_fits[], 2, sd)
names(avg_prediction)[7:12] = c("King_Haakon_sd", "Ross_sd", "East_Antarctica_sd", "Weddell_sd",
                                "Bellingshausen_Amundsen_sd", "Total_sd")
avg_prediction_e <- data.frame(
            "Total"= c(avg_prediction$Total, all_dat$total[1465:1500]),
            "King_Haakon" = c(avg_prediction$King_Haakon, all_dat$King_Hakon[1465:1500]),
            "Ross" = c(avg_prediction$Ross, all_dat$Ross[1465:1500]),
            "East_Antarctica" = c(avg_prediction$East_Antarctica, all_dat$East_Antarctica[1465:1500]),
            "Weddell" = c(avg_prediction$Weddell, all_dat$Weddell[1465:1500]),
            "Bellingshausen_Amundsen" = c(avg_prediction$Bellingshausen_Amundsen, all_dat$Bellingshausen_Amundsen_Sea[1465:1500]) )
# get a prediction for the total
avg_prediction_e$Month = 1:12
avg_prediction_e$Season = c(1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 1, 1)
avg_prediction_e$Year = rep(1899:2023, each = 12)
avg_prediction_e$Decade = as.character(rep(seq(1900, 2023, by = 10), each = 12*10))[1:1500]

pdf(file = "../Figures/Figure2-no-legend.pdf", width = 8, height = 8)
library(mgcv)
fit_total <- gam(Total ~ s(Month, Year), data=avg_prediction_e, method="ML")
summary(fit_total)

grd <- rep(1899:2023, each = 12/0.2)
grdM <- seq(0,12-0.2,by=0.2)
grd <- data.frame( Month=rep(grdM, length(grd)/length(grdM)), Year=grd )
p_cplot = data.frame("Month"=grd$Month, "Year"=grd$Year,
            "Total"= grd$Month,
            "King_Haakon" = grd$Month,
            "Ross" = grd$Month,
            "East_Antarctica" = grd$Month,
            "Weddell" = grd$Month,
            "Bellingshausen_Amundsen" = grd$Month )

r <- function(x){x}
p_cplot[,"Total"]=r(predict(gam(Total ~ te(Month, Year, bs = c("cc", "cs"), k=c(12,20)), data=avg_prediction_e, method="ML"), newdata=grd))
p_cplot[,"King_Haakon"]=r(predict(gam(King_Haakon ~ te(Month, Year, bs = c("cc", "cs"), k=c(12,20)), data=avg_prediction_e, method="ML"), newdata=grd))
p_cplot[,"Ross"]=r(predict(gam(Ross ~ te(Month, Year, bs = c("cc", "cs"), k=c(12,20)), data=avg_prediction_e, method="ML"),
newdata=grd))
p_cplot[,"East_Antarctica"]=r(predict(gam(East_Antarctica ~ te(Month, Year, bs = c("cc", "cs"), k=c(12,20)), data=avg_prediction_e, method="ML"), newdata=grd))
p_cplot[,"Weddell"]=r(predict(gam(Weddell ~ te(Month, Year, bs = c("cc", "cs"), k=c(12,20)), data=avg_prediction_e, method="ML"), newdata=grd))
p_cplot[,"Bellingshausen_Amundsen"]=r(predict(gam(Bellingshausen_Amundsen ~ te(Month, Year, bs = c("cc", "cs"), k=c(12,20)), data=avg_prediction_e, method="ML"), newdata=grd))

p_stack <- data.frame(Month=rep(p_cplot$Month, 6), Year=rep(p_cplot$Year,6),
  Sector=rep(c("Total","King Haakon VII","Ross-Amundsen","East Antarctica","Weddell","Bellingshausen Amundsen"),each=length(p_cplot$Total)),
  Anomaly=c(p_cplot[,"Total"],p_cplot[,"King_Haakon"],p_cplot[,"Ross"],p_cplot[,"East_Antarctica"],
            p_cplot[,"Weddell"],p_cplot[,"Bellingshausen_Amundsen"]) )

# Heat map colors
range(p_stack$Anomaly)
quantile(p_stack$Anomaly, p=(1:9)/10)
mybreaks <- seq(min(p_stack$Anomaly), max(p_stack$Anomaly), length.out = 13)

mycolors<- function(x) {
   colors<-colorRampPalette(c("darkblue", "yellow"))( 12 )
   colors[1:x]
}
breaklabel <- function(x){
   labels<- paste0(mybreaks[1:12], "-", mybreaks[2:13])
   labels<- paste0(mybreaks[1:12], "-", mybreaks[2:13])
   labels[1:x]
}
cols=c("purple","darkblue","blue","dodgerblue3","dodgerblue","deepskyblue2","deepskyblue","skyblue","lightskyblue2","lightskyblue3","lightgoldenrod1","khaki1","yellow","goldenrod1","orange","darkorange","firebrick2","red","red2","brown2","indianred","red3")
cols <- rev(cols)
mybreaks <- seq(-0.3, 0.6, length.out = length(cols))
mybreaks <- c(seq(-0.3, 0, length.out = 10), seq(0, 0.6, length=12))
mybreaks <- c(seq(-0.5, 0, length.out = 10), seq(0, 0.6, length=12))
mybreaks <- c(seq(-0.85, 0, length.out = 10), seq(0, 0.6, length=12))
mycolors<- function(x) {
cols
}
base_size <- 5.5
p <- ggplot(p_stack, aes(Year, Month, z=Anomaly)) + theme_minimal() +
 #scale_fill_viridis_c(direction=-1) 
# scale_fill_gradientn(colors=c("purple","darkblue","blue","dodgerblue3","dodgerblue","deepskyblue2","deepskyblue","skyblue","lightskyblue2","lightskyblue3","lightgoldenrod1","khaki1","yellow","goldenrod1","orange","darkorange","firebrick2","red","red2","brown2","indianred","red3"), limits=range(p_stack$Anomaly))
 #scale_fill_manual(values=cols, breaks=seq(-0.85,0.6,length=length(cols)+2)) +
 #scale_color_continuous(breaks = c(-0.5, 0, 0.5), labels = c("low", "med", "high")) +
  theme(plot.margin = margin(base_size / 2, 10*base_size, base_size / 2, base_size / 2, unit = "pt")) +
  geom_contour_filled(breaks=mybreaks, show.legend = FALSE) +
  scale_fill_manual(name="Anomaly", values=mycolors(1), drop=FALSE) #, breaks=seq(-0.85,0.6,length=length(cols)+2)) +
p <- p + scale_y_continuous(breaks=1:12, label=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
p <- p + scale_x_continuous(breaks=seq(1900,2020, by=10), label=c("1900","","1920","","1940","","1960","","1980","","2000","","2020"))
print(p + facet_wrap(~factor(Sector, levels=c("Total","King Haakon VII","Ross-Amundsen","East Antarctica","Weddell","Bellingshausen Amundsen")), ncol=2))
dev.off()

tdate.m <- seq(1899+0.5/12,2025,by=1/12)[1:1495]
load(file="../Data/persistence.RData")

pdf(paste0("../Figures/Fig3_persistence.pdf",sep=""),paper="USr",width=11,height=9.5)

#Me=apply(persistence.e,2,mean)
#Be=apply(persistence.e,2,quantile,c(0.05,0.95))
plot(x=tdate.m, y=Me,type='l',xlab="Time (years)",ylab="persistence",ylim=c(min(Be[1,]),max(Be[2,])),
     xlim=c(1900,2025), main="", lwd=3, xaxt="n")
lines(x=tdate.m, y=Be[1,],col=2,lty=2, lwd=3)
lines(x=tdate.m, y=Be[2,],col=2,lty=2, lwd=3)
abline(h=1, lty=1, col="gray")
axis(1, at=c(seq(1900, 2020, by=5), 2024))
abline(h=0, lty=1, col="gray")
abline(v=c(1978.9, 1987.9, 2020, 2024.8), lty=1, col="gray")

#Ms=apply(persistence.s,2,mean)
#Bs=apply(persistence.s,2,quantile,c(0.05,0.95))
plot(x=tdate.m, y=Ms,type='l',xlab="Time (years)",ylab="persistence",
 ylim=c(-1,1),
     xlim=c(1900,2025), main="", lwd=3, xaxt="n")
lines(x=tdate.m, y=Bs[1,],col=2,lty=2, lwd=3)
lines(x=tdate.m, y=Bs[2,],col=2,lty=2, lwd=3)
abline(h=1, lty=1, col="gray")
axis(1, at=c(seq(1900, 2020, by=5), 2024))
abline(h=0, lty=1, col="gray")
abline(v=c(1978.9, 1987.9, 2020, 2024.8), lty=1, col="gray")

dev.off()
#
set.seed(1)
region_names <- c("Total","King Haakon", "Ross", "East Antarctica", "Weddell", "Bellingshausen Amundsen")

pdf("../Figures/Supp_Fig_2_cor_obs.pdf", width=7, height=4)

start <- 1033 # 1/1985
end <- 1404 # 12/2015

# total 
#load(file="../Data/pos_pred_fits.RData")
#load(file="pos_pred_total.RData")
#total_fits <- Y_pp
#all_dat$total <- Y_obs

# sectors 
#load(file="pos_pred.RData")
sie_rec <- abind::abind(total_PP[,start:end], sector_PP[,start:end,])

# the reconstructions during the satelitte era
mth_rec <- all_dat[start:end,"Month"]
year_rec <- all_dat[start:end,"Year"]
all_dat_subset <- all_dat[start:end,]
start <- 1
end <- nrow(all_dat_subset)
# From 1/1985 to 12/2015
sie <- all_dat_subset[start:end,]
mth <- all_dat_subset[start:end,"Month"]
year <- all_dat_subset[start:end,"Year"]
cor_obs = matrix(1.0, ncol=12, nrow=12)
dimnames(cor_obs) <- list(c("J", "F", "M", "A", "M ", " J ", "J ", "A ", "S", "O", "N", "D"),
                          c("J", "F", "M", "A", "M ", " J ", "J ", "A ", "S", "O", "N", "D")
                         )
region <- 5 # Weddell

for (i in 1:12){
  for (j in 1:12){
    if(i!=j){
      v <- NULL
      iyear <- year[mth==i]
      jyear <- year[mth==j]
      imth <- sie[mth==i,2+region]
      jmth <- sie[mth==j,2+region]
      for(k in seq_along(iyear)){
       for(l in seq_along(jyear)){
        if(i < j & jyear[l] == iyear[k]){
         v <- rbind(v, c(imth[k], jmth[l]))
        }
       }
       for(l in seq_along(jyear)){
         if(i > j & jyear[l] == (iyear[k]+1)){
          v <- rbind(v, c(imth[k], jmth[l]))
         }
       }
      }
      cor_obs[i,j] <- cor(v[,1],v[,2])
    }
  }
}
par(mfrow=c(1,1))
melted_cor_obs <- reshape2::melt(round(cor_obs,2))
print( ggplot(data = melted_cor_obs, aes(x=Var2, y=Var1, fill=value)) + 
  xlab("") + ylab("") +
  scale_fill_gradient2("Obs Corr", limits = c(-1, 1), high = "firebrick", low = "steelblue") +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle(region_names[region]) +
geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank() ) )
dev.off()

pdf("../Figures/Supp_Fig_2_cor_pp.pdf", width=7, height=4)
cor_rec = matrix(1.0, ncol=12, nrow=12)
dimnames(cor_rec) <- list(c("J", "F", "M", "A", "M ", " J ", "J ", "A ", "S", "O", "N", "D"),
                          c("J", "F", "M", "A", "M ", " J ", "J ", "A ", "S", "O", "N", "D")
                         )
for (i in 1:12){
  for (j in 1:12){
    if(i!=j){
      v <- NULL
      iyear <- year_rec[mth_rec==i]
      jyear <- year_rec[mth_rec==j]
      imth <- sie_rec[,mth_rec==i,region]
      jmth <- sie_rec[,mth_rec==j,region]
      for(k in seq_along(iyear)){
       for(l in seq_along(jyear)){
        if(i < j & jyear[l] == iyear[k]){
         v <- rbind(v, cbind(imth[,k], jmth[,l]))
        }
       }
       for(l in seq_along(jyear)){
         if(i > j & jyear[l] == (iyear[k]+1)){
          v <- rbind(v, cbind(imth[,k], jmth[,l]))
         }
       }
      }
      cor_rec[i,j] <- cor(v[,1],v[,2])
    # print(nrow(v))
    }
  }
}
# plot correlation matrix
melted_cor_rec <- reshape2::melt(round(cor_rec,2))
print( ggplot(data = melted_cor_rec, aes(x=Var2, y=Var1, fill=value)) + 
  xlab("") + ylab("") +
  scale_fill_gradient2("Rec Corr", limits = c(-1, 1), high = "firebrick", low = "steelblue") +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle(region_names[region]) +
geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank() ) )

dev.off()

pdf("../Figures/Supp_Fig_2_cor_obs_1978-2006.pdf", width=7, height=4)

start <- 959
end <- 1296 # 12/2006

# total 
#load(file="../Data/pos_pred_fits.RData")
#load(file="pos_pred_total_fits.RData")
#total_fits <- Y_pp
#all_dat$total <- Y_obs

# sectors 
#load(file="pos_pred_fits.RData")
sie_rec <- abind::abind(total_PP[,start:end], sector_PP[,start:end,])

# the reconstructions during the satelitte era
mth_rec <- all_dat[start:end,"Month"]
year_rec <- all_dat[start:end,"Year"]
all_dat_subset <- all_dat[start:end,]
start <- 1
end <- nrow(all_dat_subset)
# From 1/1985 to 12/2015
sie <- all_dat_subset[start:end,]
mth <- all_dat_subset[start:end,"Month"]
year <- all_dat_subset[start:end,"Year"]
cor_obs = matrix(1.0, ncol=12, nrow=12)
dimnames(cor_obs) <- list(c("J", "F", "M", "A", "M ", " J ", "J ", "A ", "S", "O", "N", "D"),
                          c("J", "F", "M", "A", "M ", " J ", "J ", "A ", "S", "O", "N", "D")
                         )
region <- 5 # Weddell

for (i in 1:12){
  for (j in 1:12){
    if(i!=j){
      v <- NULL
      iyear <- year[mth==i]
      jyear <- year[mth==j]
      imth <- sie[mth==i,2+region]
      jmth <- sie[mth==j,2+region]
      for(k in seq_along(iyear)){
       for(l in seq_along(jyear)){
        if(i < j & jyear[l] == iyear[k]){
         v <- rbind(v, c(imth[k], jmth[l]))
        }
       }
       for(l in seq_along(jyear)){
         if(i > j & jyear[l] == (iyear[k]+1)){
          v <- rbind(v, c(imth[k], jmth[l]))
         }
       }
      }
      cor_obs[i,j] <- cor(v[,1],v[,2])
    }
  }
}
par(mfrow=c(1,1))
melted_cor_obs <- reshape2::melt(round(cor_obs,2))
print( ggplot(data = melted_cor_obs, aes(x=Var2, y=Var1, fill=value)) + 
  xlab("") + ylab("") +
  scale_fill_gradient2("Obs Corr", limits = c(-1, 1), high = "firebrick", low = "steelblue") +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle(region_names[region]) +
geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank() ) )
dev.off()

pdf("../Figures/Supp_Fig_2_cor_pp_1978-2006.pdf", width=7, height=4)
cor_rec = matrix(1.0, ncol=12, nrow=12)
dimnames(cor_rec) <- list(c("J", "F", "M", "A", "M ", " J ", "J ", "A ", "S", "O", "N", "D"),
                          c("J", "F", "M", "A", "M ", " J ", "J ", "A ", "S", "O", "N", "D")
                         )
for (i in 1:12){
  for (j in 1:12){
    if(i!=j){
      v <- NULL
      iyear <- year_rec[mth_rec==i]
      jyear <- year_rec[mth_rec==j]
      imth <- sie_rec[,mth_rec==i,region]
      jmth <- sie_rec[,mth_rec==j,region]
      for(k in seq_along(iyear)){
       for(l in seq_along(jyear)){
        if(i < j & jyear[l] == iyear[k]){
         v <- rbind(v, cbind(imth[,k], jmth[,l]))
        }
       }
       for(l in seq_along(jyear)){
         if(i > j & jyear[l] == (iyear[k]+1)){
          v <- rbind(v, cbind(imth[,k], jmth[,l]))
         }
       }
      }
      cor_rec[i,j] <- cor(v[,1],v[,2])
    # print(nrow(v))
    }
  }
}
# plot correlation matrix
melted_cor_rec <- reshape2::melt(round(cor_rec,2))
print( ggplot(data = melted_cor_rec, aes(x=Var2, y=Var1, fill=value)) + 
  xlab("") + ylab("") +
  scale_fill_gradient2("Rec Corr", limits = c(-1, 1), high = "firebrick", low = "steelblue") +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle(region_names[region]) +
geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank() ) )

dev.off()

pdf("../Figures/Supp_Fig_2_cor_obs_2007-2023.pdf", width=7, height=4)

start <- 1297 # 1/2007
end <- 1500 # 12/2023

# total 
#load(file="../Data/pos_pred_fits.RData")
#load(file="pos_pred_total_fits.RData")
#total_fits <- Y_pp
#all_dat$total <- Y_obs

# sectors 
#load(file="pos_pred_fits.RData")
sie_rec <- abind::abind(total_PP[,start:end], sector_PP[,start:end,])

# the reconstructions during the satelitte era
mth_rec <- all_dat[start:end,"Month"]
year_rec <- all_dat[start:end,"Year"]
all_dat_subset <- all_dat[start:end,]
start <- 1
end <- nrow(all_dat_subset)
# From 1/1985 to 12/2015
sie <- all_dat_subset[start:end,]
mth <- all_dat_subset[start:end,"Month"]
year <- all_dat_subset[start:end,"Year"]
cor_obs = matrix(1.0, ncol=12, nrow=12)
dimnames(cor_obs) <- list(c("J", "F", "M", "A", "M ", " J ", "J ", "A ", "S", "O", "N", "D"),
                          c("J", "F", "M", "A", "M ", " J ", "J ", "A ", "S", "O", "N", "D")
                         )
region <- 5 # Weddell

for (i in 1:12){
  for (j in 1:12){
    if(i!=j){
      v <- NULL
      iyear <- year[mth==i]
      jyear <- year[mth==j]
      imth <- sie[mth==i,2+region]
      jmth <- sie[mth==j,2+region]
      for(k in seq_along(iyear)){
       for(l in seq_along(jyear)){
        if(i < j & jyear[l] == iyear[k]){
         v <- rbind(v, c(imth[k], jmth[l]))
        }
       }
       for(l in seq_along(jyear)){
         if(i > j & jyear[l] == (iyear[k]+1)){
          v <- rbind(v, c(imth[k], jmth[l]))
         }
       }
      }
      cor_obs[i,j] <- cor(v[,1],v[,2])
    }
  }
}
par(mfrow=c(1,1))
melted_cor_obs <- reshape2::melt(round(cor_obs,2))
print( ggplot(data = melted_cor_obs, aes(x=Var2, y=Var1, fill=value)) + 
  xlab("") + ylab("") +
  scale_fill_gradient2("Obs Corr", limits = c(-1, 1), high = "firebrick", low = "steelblue") +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle(region_names[region]) +
geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank() ) )
dev.off()

pdf("../Figures/Supp_Fig_2_cor_pp_2007-2023.pdf", width=7, height=4)
cor_rec = matrix(1.0, ncol=12, nrow=12)
dimnames(cor_rec) <- list(c("J", "F", "M", "A", "M ", " J ", "J ", "A ", "S", "O", "N", "D"),
                          c("J", "F", "M", "A", "M ", " J ", "J ", "A ", "S", "O", "N", "D")
                         )
for (i in 1:12){
  for (j in 1:12){
    if(i!=j){
      v <- NULL
      iyear <- year_rec[mth_rec==i]
      jyear <- year_rec[mth_rec==j]
      imth <- sie_rec[,mth_rec==i,region]
      jmth <- sie_rec[,mth_rec==j,region]
      for(k in seq_along(iyear)){
       for(l in seq_along(jyear)){
        if(i < j & jyear[l] == iyear[k]){
         v <- rbind(v, cbind(imth[,k], jmth[,l]))
        }
       }
       for(l in seq_along(jyear)){
         if(i > j & jyear[l] == (iyear[k]+1)){
          v <- rbind(v, cbind(imth[,k], jmth[,l]))
         }
       }
      }
      cor_rec[i,j] <- cor(v[,1],v[,2])
    # print(nrow(v))
    }
  }
}

# plot correlation matrix
melted_cor_rec <- reshape2::melt(round(cor_rec,2))
print( ggplot(data = melted_cor_rec, aes(x=Var2, y=Var1, fill=value)) + 
  xlab("") + ylab("") +
  scale_fill_gradient2("Rec Corr", limits = c(-1, 1), high = "firebrick", low = "steelblue") +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle(region_names[region]) +
geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank() ) )

dev.off()

# persistence plots
pdf(file = "../Figures/Supp_Fig_8.pdf", width = 8, height = 8)
tdate.m <- seq(1899+0.5/12,2025,by=1/12)[1:1495]
par(mfrow=c(3,2))
# c(bottom, left, top, right)
par(mar=c(0,4.2,0,0))
#
# King_Hakon
load("../Data/persistence_1.RData")
par(mar=c(0,2,0,0))
B=apply(persistence.e,2,quantile,c(0.05,0.95))
plot(x=tdate.m, y=apply(persistence.e,2,mean),type='l',xlab="Time (years)",ylab="persistence",ylim=c(0.75,1),
     xlim=c(1900,2025), main="", lwd=3, axes=FALSE)
lines(x=tdate.m, y=B[1,],col=2,lty=2, lwd=3)
lines(x=tdate.m, y=B[2,],col=2,lty=2, lwd=3)
abline(h=1, lty=1, col="gray")
axis(1, at=c(seq(1900, 2020, by=5), 2024))
abline(h=0, lty=1, col="gray")
abline(v=c(1978.9, 1987.9, 2020, 2024.8), lty=1, col="gray")
axis(2)
text(x=1920, y=0.98, labels="King Haakon VII", col=2)
abline(h=0, col="gray")

# Ross
par(mar=c(0,4.2,0,0))
load("../Data/persistence_2.RData")
par(mar=c(0,0,0,2))
B=apply(persistence.e,2,quantile,c(0.05,0.95))
plot(x=tdate.m, y=apply(persistence.e,2,mean),type='l',xlab="Time (years)",ylab="persistence",ylim=c(0.75,1),
     xlim=c(1900,2025), main="", lwd=3, axes=FALSE)
lines(x=tdate.m, y=B[1,],col=2,lty=2, lwd=3)
lines(x=tdate.m, y=B[2,],col=2,lty=2, lwd=3)
abline(h=1, lty=1, col="gray")
axis(1, at=c(seq(1900, 2020, by=5), 2024))
abline(h=0, lty=1, col="gray")
abline(v=c(1978.9, 1987.9, 2020, 2024.8), lty=1, col="gray")
axis(4)
text(x=2007, y=0.98, labels="Ross Sea", col=2)
abline(h=0, col="gray")

# East_Antarctica
load("../Data/persistence_3.RData")
par(mar=c(0,2,0,0))
B=apply(persistence.e,2,quantile,c(0.05,0.95))
plot(x=tdate.m, y=apply(persistence.e,2,mean),type='l',xlab="Time (years)",ylab="persistence",ylim=c(0.75,1),
     xlim=c(1900,2025), main="", lwd=3, axes=FALSE)
lines(x=tdate.m, y=B[1,],col=2,lty=2, lwd=3)
lines(x=tdate.m, y=B[2,],col=2,lty=2, lwd=3)
abline(h=1, lty=1, col="gray")
axis(1, at=c(seq(1900, 2020, by=5), 2024))
abline(h=0, lty=1, col="gray")
abline(v=c(1978.9, 1987.9, 2020, 2024.8), lty=1, col="gray")
axis(2)
text(x=1920, y=0.75, labels="East Antarctica", col=2)
abline(h=0, col="gray")

# Weddell
load("../Data/persistence_4.RData")
par(mar=c(0,0,0,2))
B=apply(persistence.e,2,quantile,c(0.05,0.95))
plot(x=tdate.m, y=apply(persistence.e,2,mean),type='l',xlab="Time (years)",ylab="persistence",ylim=c(0.75,1),
     xlim=c(1900,2025), main="", lwd=3, axes=FALSE)
lines(x=tdate.m, y=B[1,],col=2,lty=2, lwd=3)
lines(x=tdate.m, y=B[2,],col=2,lty=2, lwd=3)
abline(h=1, lty=1, col="gray")
axis(1, at=c(seq(1900, 2020, by=5), 2024))
abline(h=0, lty=1, col="gray")
abline(v=c(1978.9, 1987.9, 2020, 2024.8), lty=1, col="gray")
axis(4)
text(x=1910, y=0.75, labels="Weddell", col=2)
abline(h=0, col="gray")

# Bellingshausen_Amundsen_Sea
load("../Data/persistence_5.RData")
par(mar=c(0,2,0,0))
B=apply(persistence.e,2,quantile,c(0.05,0.95))
plot(x=tdate.m, y=apply(persistence.e,2,mean),type='l',xlab="Time (years)",ylab="persistence",ylim=c(0.75,1),
     xlim=c(1900,2025), main="", lwd=3, axes=FALSE)
lines(x=tdate.m, y=B[1,],col=2,lty=2, lwd=3)
lines(x=tdate.m, y=B[2,],col=2,lty=2, lwd=3)
abline(h=1, lty=1, col="gray")
axis(1, at=c(seq(1900, 2020, by=5), 2024))
abline(h=0, lty=1, col="gray")
abline(v=c(1978.9, 1987.9, 2020, 2024.8), lty=1, col="gray")
axis(2)
text(x=1930, y=0.76, labels="Bellinghausen Amundsen Sea", col=2, cex=0.95)
abline(h=0, col="gray")
dev.off()

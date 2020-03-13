library(dplyr)
library(pracma)
library(ggplot2)
library(reshape2)
library(tidyr)
library(gtools)
library(nlme)
library(ggpubr)
library(ncdf4)
library(chron)
library(maps)
library(fields)
library(maptools)
library(mapdata)
library(zoo)
library(oce)

# # I have catch data in spearate files for each management area
# # clean up and combine data sets
# 
# ch <- read.csv("Chignik catch.csv")
# 
# # drop Chinook and change to tall
# ch <- ch %>%
#   select(-Chinook) %>%
#   gather(species, catch, -Year)
# 
# ch$area <- "Chignik"
# 
# kd <- read.csv("Kodiak catch.csv")
# 
# # drop Chinook and change to tall
# kd <- kd %>%
#   select(-Chinook) %>%
#   gather(species, catch, -Year)
# 
# kd$area <- "Kodiak"
# 
# sp <- read.csv("South Alaska Peninsula catch.csv")
# 
# colnames(sp) <- c("Year", "Chinook", "Sockeye", "Coho", "Pink", "Chum")
# 
# # drop Chinook and change to tall
# sp <- sp %>%
#   select(-Chinook) %>%
#   gather(species, catch, -Year)
# 
# sp$area <- "S. Peninsula"
# 
# ci <- read.csv("Cook Inlet catch.csv")
# 
# colnames(ci) <- c("Year", "Sockeye", "Coho", "Pink", "Chum")
# 
# # change to tall
# ci <- ci %>%
#   gather(species, catch, -Year)
# 
# ci$area <- "Cook Inlet"
# 
# pws <- read.csv("Prince William Sound catch.csv")
# 
# # drop Chinook and change to tall
# pws <- pws %>%
#   select(-Chinook) %>%
#   gather(species, catch, -Year)
# 
# pws$area <- "Prince William Sound"
# 
# se <- read.csv("Southeast Alaska catch.csv")
# 
# # drop Chinook and change to tall
# se <- se %>%
#   select(-Chinook) %>%
#   gather(species, catch, -Year)
# 
# se$area <- "Southeast"
# 
# # combine and limit to 1965-present
# raw.dat <- rbind(sp, ch, kd, ci, pws, se)
# raw.dat <- filter(raw.dat, Year >= 1965)
# 
# # change 0s (2019) to NA!
# change <- raw.dat$catch==0
# raw.dat$catch[change] <- NA
# 
# # but 2018 Chignik catch is really 0!
# raw.dat$catch[raw.dat$Year==2018 & raw.dat$area=="Chignik"] <- 0
# 
# raw.dat <- raw.dat %>%
#   group_by(species, Year) %>%
#   summarise(log(sum(catch), 10)) 
# colnames(raw.dat)[3] <- "log.catch"
# 
# ggplot(raw.dat, aes(Year, log.catch, color=species)) +
#   theme_bw() +
#   geom_line()
# 
# # remove Chum
# raw.dat <- raw.dat %>%
#   filter(species != "Chum")
# 
# raw.dat$even.odd <- ifelse(odd(raw.dat$Year)==T,"odd", "even")
# 
# raw.dat$species.plot <- paste(raw.dat$species, raw.dat$even.odd, sep="-")
# 
# raw.dat$species.plot <- ifelse(raw.dat$species=="Pink", raw.dat$species.plot, raw.dat$species)
# 
# ggplot(raw.dat, aes(Year, log.catch, color=species.plot)) +
#   theme_bw() +
#   geom_line()
# 
# # don't know why I'm having touble with tidyverse!
# raw.dat$species <- as.factor(raw.dat$species)
# raw.dat$species.plot <- as.factor(raw.dat$species.plot)
# 
# spp <- levels(raw.dat$species.plot)
# raw.dat <- raw.dat %>%
#   select(Year, species.plot, log.catch)
# 
# raw.dat <- raw.dat[,2:4]
# 
# raw.dat <- raw.dat %>%
#   spread(species.plot, log.catch)
# 
# raw.dat[,2:5] <- scale(raw.dat[,2:5])
# 
# raw.dat <- raw.dat %>%
#   gather(key="species", value="catch", -Year)
# 
# raw.dat <- na.omit(raw.dat)
# 
# # now lag to mean ocean entry for each species
# # catch-1 for pink/coho, catch-2 for sockeye
# raw.dat$entry.year <- ifelse(raw.dat$species %in% c("Pink-even", "Pink-odd", "Coho"), raw.dat$Year-1, raw.dat$Year-2)

raw.dat <- read.csv("salmon.and.covariate.data.csv")

# note that years in this file are already lagged to ocean entry!

# load ERSST
# uncomment these lines to download data
# # identify latest year and month needed
# year <- 2019
# month <- "07"
# 
# URL <- paste("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1854-01-01):1:(", year, "-", month, "-01T00:00:00Z)][(0.0):1:(0.0)][(20):1:(70)][(120):1:(250)]", sep="")
# 
# download.file(URL, "data/North.Pacific.ersst")

# open netcdf file of SST data
nc <- nc_open("/Users/MikeLitzow 1/Documents/R/climate-data/data/North.Pacific.ersst")

# extract dates
ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
sst.d <- dates(h, origin = c(1,1,1970))

sst.x <- ncvar_get(nc, "longitude")
sst.y <- ncvar_get(nc, "latitude")

SST <- ncvar_get(nc,  "sst")
# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SST <- aperm(SST, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
sst.lat <- rep(sst.y, length(sst.x))   
sst.lon <- rep(sst.x, each = length(sst.y))   
dimnames(SST) <- list(as.character(sst.d), paste("N", sst.lat, "E", sst.lon, sep=""))

# plot to check
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=tim.colors(64))
contour(sst.x, sst.y, z, add=T) 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# a-ok!

# now get NDJFM means for each cell
m <- months(sst.d)
yr <- as.numeric(as.character(years(sst.d)))

win.months <- m[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
win.yrs <- ifelse(m %in% c("Nov", "Dec"), yr+1, yr)
win.yrs <- win.yrs[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]

# get seperate matrix for winter
win.SST <- SST[m %in% win.months,]

f <- function(x) tapply(x, win.yrs, mean)
win.SST <- apply(win.SST, 2, f)

# and smooth each cell as 3-yr running mean
# I'm doing this b/c analysis shows that catch most strongly 
# responds to SST over 3-yr time scales
f <- function(x) rollmean(x, 3, fill=NA)
win.SST <- apply(win.SST, 2, f)

# plot to check
SST.mean <- colMeans(win.SST, na.rm=T)
z <- t(matrix(SST.mean,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=tim.colors(64))
contour(sst.x, sst.y, z, add=T) 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# all good!
# now make correlation plots for each spp. in each era...

# first GOA sockeye
winter.sockeye.65.88 <- winter.sockeye.89.13 <- winter.sockeye.14.19 <- NA

for(j in 1:ncol(win.SST)){
  # note that we are using catch year as the nominal year matching the proposed eras
   # j <- 1
  winter.sockeye.65.88[j] <- 
    cor(win.SST[rownames(win.SST) %in% 1963:1986, j], 
        raw.dat$catch[raw.dat$species=="Sockeye" & 
                                      raw.dat$Year %in% 1963:1986])
  
  winter.sockeye.89.13[j] <- 
    cor(win.SST[rownames(win.SST) %in% 1987:2011, j], 
        raw.dat$catch[raw.dat$species=="Sockeye" & 
                        raw.dat$Year %in% 1987:2011])
  
  winter.sockeye.14.19[j] <- 
    cor(win.SST[rownames(win.SST) %in% 2012:2017, j], 
        raw.dat$catch[raw.dat$species=="Sockeye" & 
                        raw.dat$Year %in% 2012:2017])
}

# now coho!
winter.coho.65.88 <- winter.coho.89.13 <- winter.coho.14.19 <- NA

for(j in 1:ncol(win.SST)){
  # note that we are using catch year as the nominal year matching the proposed eras
  # j <- 1
  winter.coho.65.88[j] <- 
    cor(win.SST[rownames(win.SST) %in% 1964:1987, j], 
        raw.dat$catch[raw.dat$species=="Coho" & 
                        raw.dat$Year %in% 1964:1987])
  
  winter.coho.89.13[j] <- 
    cor(win.SST[rownames(win.SST) %in% 1988:2012, j], 
        raw.dat$catch[raw.dat$species=="Coho" & 
                        raw.dat$Year %in% 1988:2012])
  
  winter.coho.14.19[j] <- 
    cor(win.SST[rownames(win.SST) %in% 2013:2018, j], 
        raw.dat$catch[raw.dat$species=="Coho" & 
                        raw.dat$Year %in% 2013:2018])
}


# now pink!
# combining odd and even for our purposes...
winter.pink.65.88 <- winter.pink.89.13 <- winter.pink.14.19 <- NA

pink.combined <- raw.dat %>%
  filter(species %in% c("Pink-odd", "Pink-even")) %>%
  select(Year, species, catch) %>%
  arrange(Year)

pink.combined <- na.omit(pink.combined)
ggplot(pink.combined, aes(Year, catch)) +
  geom_line()

for(j in 1:ncol(win.SST)){
  # note that we are using catch year as the nominal year matching the proposed eras
  # j <- 1
  winter.pink.65.88[j] <- 
    cor(win.SST[rownames(win.SST) %in% 1964:1987, j], 
        pink.combined$catch[pink.combined$Year %in% 1964:1987])
  
  winter.pink.89.13[j] <- 
    cor(win.SST[rownames(win.SST) %in% 1988:2012, j], 
        pink.combined$catch[pink.combined$Year %in% 1988:2012])
  
  winter.pink.14.19[j] <- 
    cor(win.SST[rownames(win.SST) %in% 2013:2018, j], 
        pink.combined$catch[pink.combined$Year %in% 2013:2018])
}


png("correlations with 3-yr smoothed winter sst sockeye pink coho.png", 6, 6, units="in", res=300)

par(mfrow=c(3,3), mar=c(0,0.5,1.5,0.5), oma=c(2,2,2,0))
# new.col <- tim.colors(64)
# grays <- c("gray90", "gray89", "gray88",
#            "gray87","gray86")
# 
# new.col[27:36] <- c(grays[5:1], grays[1:5])

new.col <- oceColorsPalette(64)

lim <- c(-1,1)
# first, sockeye
z <- t(matrix(winter.sockeye.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")
# map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("Sockeye", outer=T, cex=1.2, side=2, adj=0.88)
mtext("1965-1988", outer=T, cex=1.2, side=3, adj=0.1, line=-1)

z <- t(matrix(winter.sockeye.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")
mtext("1989-2013", outer=T, cex=1.2, side=3, adj=0.5, line=-1)

z <- t(matrix(winter.sockeye.14.19,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")
mtext("2014-2019", outer=T, cex=1.2, side=3, adj=0.9, line=-1)

# now pink
par(mar=c(0.5,0.5,0.5,0.5))

z <- t(matrix(winter.pink.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")
mtext("Pink", outer=T, cex=1.2, side=2, adj=0.5)

z <- t(matrix(winter.pink.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

z <- t(matrix(winter.pink.14.19,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

# coho
par(mar=c(1.5,0.5,0,0.5))
z <- t(matrix(winter.coho.65.88,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")
mtext("Coho", outer=T, cex=1.2, side=2, adj=0.15)

z <- t(matrix(winter.coho.89.13,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

# add legend strip
mt.cex <- 1.1
l.mar <- 3
l.cex <- 1
l.l <- 1.2
tc.l <- -0.2

z[1,1] <- -1; z[66,26] <- 1
image.plot(z, legend.only=TRUE, horizontal =TRUE,  legend.lab = "r", 
           smallplot = c(0,1,0.03,0.08), 
           legend.cex=1, col=new.col,
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0))) 

z <- t(matrix(winter.coho.14.19,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=new.col, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")
contour(sst.x,sst.y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")


dev.off()

library(oce)

col <- oceColorsPalette(64)

z <- t(matrix(winter.coho.14.19,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
imagep(sst.x,sst.y,z, zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")

image(sst.x,sst.y,z, col=col, 
      zlim=lim, ylim=c(20,66), yaxt="n", xaxt="n", xlab="", ylab="")


# code for processing CMIP5 preindustrial SST estimates
# and comparing these to observations

library(tidyverse)
library(ncdf4)
library(maps)
library(mapdata)
library(fields)
library(chron)

# load CMIP5 summaries from John Walsh

dat <- read.csv("CMIP5 GOA SST.csv")

dat <- dat %>%
  gather(model, anomaly, -Year, -Era)

ggplot(dat, aes(Year, anomaly, color=model)) +
  theme_bw() +
  geom_line() + 
  facet_wrap(~Era, scales="free_x")

# now load ERSSTv5 for the same box, and calculate annual anomalies

# 50º-60ºN, 210º-230ºE, 1854-present

# identify latest year and month needed
year <- 2018
month <- "12"

URL <- paste("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1854-01-01):1:(", year, "-", month, "-01T00:00:00Z)][(0.0):1:(0.0)][(50):1:(60)][(210):1:(230)]", sep="")

download.file(URL, "GOA.box.ersst")

# process
nc <- nc_open("GOA.box.ersst")

# extract dates
ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

# extract study area
x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")

SST <- ncvar_get(nc, "sst", verbose = F)

# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SST <- aperm(SST, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot to check
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=tim.colors(64), xlim=c(190,240), ylim=c(40,66))
contour(x, y, z, add=T) 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# get anomaly for 1981:2010
yr <- as.numeric(as.character(years(d)))

weights <-  sqrt(cos(lat*pi/180))
ff <- function(x) weighted.mean(x, w=weights, na.rm=T)
weighted.mean <- apply(SST, 1, ff)

annual.sst <- tapply(weighted.mean, yr, mean)

mu <- mean(annual.sst[names(annual.sst) %in% 1981:2010])
sd <- sd(annual.sst[names(annual.sst) %in% 1981:2010])

annual.anomaly <- (annual.sst - mu)/sd

annual.anomaly <- annual.anomaly[names(annual.anomaly) >=1900]

plot(names(annual.anomaly), annual.anomaly, type="o", pch=19)

# now compare thes observations with each model!
compare.dat <- dat %>%
  filter(Era=="present", Year <=2005) %>%
  mutate(observed=rep(annual.anomaly[names(annual.anomaly) %in% 1987:2005],5))

ggplot(compare.dat, aes(observed, anomaly)) +
  theme_bw() +
  geom_point() +
  facet_wrap(~model)

model.means <- tapply(compare.dat$anomaly, compare.dat$model, mean)
observed.mean <- mean(annual.anomaly[names(annual.anomaly) %in% 1987:2005])

bias <- model.means - observed.mean 

# so...we subtract this bias from each model's preindustrial simulations
# to calculate unbiased estimates of preindustrial distributions for comparison with observations

envelope <- dat %>%
  filter(Era=="preindustrial") %>%
  mutate(debiased=anomaly-bias[match(model, names(bias))])

tapply(envelope$debiased, envelope$model, mean)

envelope$debiased.3.yr <- NA

mods <- unique(envelope$model)

for(i in 1:length(mods)){
 # i <- 1
  temp <- envelope %>%
    filter(model==mods[i])
  
  envelope$debiased.3.yr[envelope$model==mods[i]] <- zoo::rollmean(temp$debiased, 3, fill=NA)
  
  envelope$anomaly.3.yr[envelope$model==mods[i]] <- zoo::rollmean(temp$anomaly, 3, fill=NA)
}

fmin <- function(x) min(x, na.rm=T)
fmax <- function(x) max(x, na.rm=T)

ranges <- data.frame(min.anomaly=fmin(envelope$anomaly),
                     max.anomaly=fmax(envelope$anomaly),
                     min.anomaly3=fmin(envelope$anomaly.3.yr),
                     max.anomaly3=fmax(envelope$anomaly.3.yr),
                     min.debiased=fmin(envelope$debiased),
                     max.debiased=fmax(envelope$debiased),
                     min.debiased3=fmin(envelope$debiased.3.yr),
                     max.debiased3=fmax(envelope$debiased.3.yr))

# now plot...
annual <- data.frame(year=names(annual.anomaly), anomaly=annual.anomaly, 
                     anomaly3=zoo::rollmean(annual.anomaly, 3, fill=NA))
annual$year <- as.numeric(as.character(annual$year))
ggplot(filter(annual, year >= 1965), aes(year, anomaly)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  geom_line(aes(year, anomaly3), color="red") +
  geom_hline(yintercept = c(ranges$min.debiased, ranges$max.debiased), lty=2) +
  geom_hline(yintercept = c(ranges$min.debiased3, ranges$max.debiased3), lty=2, color="red")

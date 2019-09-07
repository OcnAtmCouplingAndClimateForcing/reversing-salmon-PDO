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

ggplot(filter(dat, Era=="preindustrial"), aes(Year, anomaly, color=model)) +
  theme_bw() +
  geom_line()

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
m <- months(d[yr %in% 1981:2010])
# m <- months(d)
f <- function(x) tapply(x, m, mean)
mu <- apply(SST[yr %in% 1981:2010,], 2, f)	# Compute monthly means for each time series (location)
# mu <- apply(SST, 2, f)	

mu <- mu[rep(1:12, length(d)/12),] 

anom <- rowMeans(SST - mu, na.rm=T)   # Compute matrix of anomalies!

annual <- tapply(anom, yr, mean)
annual <- annual[names(annual) >=1900]
plot(names(annual), annual, type="o", pch=19)

# now compare thes observations with each model!
compare.dat <- dat %>%
  filter(Era=="present", Year <=2018) %>%
  mutate(observed=rep(annual[names(annual) %in% 1987:2018],5))

ggplot(compare.dat, aes(observed, anomaly)) +
  theme_bw() +
  geom_point() +
  facet_wrap(~model)

model.means <- tapply(compare.dat$anomaly, compare.dat$model, mean)
observed.mean <- mean(annual[names(annual) %in% 1987:2018])

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
  
}


# 
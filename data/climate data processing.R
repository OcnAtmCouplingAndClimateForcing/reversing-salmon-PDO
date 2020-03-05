library(ncdf4)
library(zoo)
library(gplots)
library(dplyr)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(tidyr)
library(nlme)
library(ggplot2)
library(oce)

# using monthly NCEP/NCAR!
nc.slp <- nc_open("/Users/MikeLitzow/Documents/R/climate-data/data/North.Pacific.NCEP.NCAR.slp")

# load ERSST
nc <- nc_open("/Users/MikeLitzow/Documents/R/climate-data/data/North.Pacific.ersst")

# process sst first
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

# now limit to GOA
# extract study area
# 54-62 deg. N, 200-226 deg. E
keep1 <- sst.lon %in% 198:226 
keep2 <- sst.lat %in% 54:62

sst.lon <- sst.lon[keep1]
sst.lat <- sst.lat[keep2]

sst.x <- sst.x[sst.x %in% 198:226]
sst.y <- sst.y[sst.y %in% 54:62]

SST <- SST[,keep1 & keep2]

# need to drop Bristol Bay cells
BB <- c("N58E200", "N58E202", "N56E200", "N56E198", "N58E198", "N60E198")
SST[,BB] <- NA

# plot to check
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=tim.colors(64))
contour(sst.x, sst.y, z, add=T) 
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)


# now get monthly anomalies wrt 1951-1980
sst.yr <- as.numeric(as.character(years(sst.d)))

sst.m <- months(sst.d)
m <- months(sst.d[sst.yr %in% 1951:1980])

f <- function(x) tapply(x, m, mean)
mu <- apply(SST[sst.yr %in% 1981:2010,], 2, f)

mu <- mu[rep(1:12, floor(length(sst.d)/12)),] 

xtra <- 12*((length(sst.d)/12)-floor(length(sst.d)/12))

mu <- rbind(mu, mu[1:xtra,])

sst.anom <- rowMeans((SST - mu), na.rm=T)   # Compute matrix of anomalies!

# plot to check
dec.yr <- sst.yr + (as.numeric(sst.m)-0.5)/12
plot(dec.yr, sst.anom, type="l", xlim=c(1951, max(dec.yr)), ylim=c(-2,2.5)) # looks good from ~1950 onwards!

# restrict SST to winter months
sst.win.yr <- sst.yr
sst.win.yr[sst.m %in% c("Nov", "Dec")] <- sst.win.yr[sst.m %in% c("Nov", "Dec")] +1

win <- c("Nov", "Dec", "Jan", "Feb", "Mar")
sst.anom <- sst.anom[sst.m %in% win]
sst.win.yr <- sst.win.yr[sst.m %in% win]
sst.anom <- tapply(sst.anom, sst.win.yr, mean) # get winter averages!

# restrict to 1960:2019
sst.anom <- sst.anom[names(sst.anom) %in% 1960:2019]

# plot to check
plot(1960:2019, sst.anom, type="l") # looks right...
lines(1960:2019, rollmean(sst.anom, 5, fill=NA), col="red", lwd=1.5)
abline(h=0, col="grey")

# now process SLP data - first, extract dates
raw <- ncvar_get(nc.slp, "time")  # seconds since 1-1-1970
h <- raw/(24*60*60)
slp.d <- dates(h, origin = c(1,1,1970))

slp.x <- ncvar_get(nc.slp, "longitude")
slp.y <- ncvar_get(nc.slp, "latitude")

SLP <- ncvar_get(nc.slp, "slp", verbose = F)
# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SLP <- aperm(SLP, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SLP <- matrix(SLP, nrow=dim(SLP)[1], ncol=prod(dim(SLP)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
slp.lat <- rep(slp.y, length(slp.x))   
slp.lon <- rep(slp.x, each = length(slp.y))   
dimnames(SLP) <- list(as.character(slp.d), paste("N", slp.lat, "E", slp.lon, sep=""))

# plot to check
SLP.mean <- colMeans(SLP)
z <- t(matrix(SLP.mean,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=tim.colors(64))
contour(slp.x, slp.y, z, add=T) 
map('world2Hires',fill=F,add=T, lwd=2)

slp.m <- months(slp.d)
slp.yr <- as.numeric(as.character(years(slp.d)))
slp.win.yr <- ifelse(slp.m %in% c("Nov", "Dec"), slp.yr+1, slp.yr)

# calculate NDJFM mean SLP for the whole area

slp.NDJFM <- SLP[slp.m %in% c("Nov", "Dec", "Jan", "Feb", "Mar"),]
slp.yr.NDJFM <- slp.win.yr[slp.m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
f <- function(x) tapply(x, slp.yr.NDJFM, mean)
slp.NDJFM <- apply(slp.NDJFM, 2, f)

# calculate NDJFM NPI - 
# 30º-65ºN, 160ºE-140ºW 

NPI <- slp.NDJFM 

drop <- slp.lat < 30
NPI[,drop] <- NA

drop <- slp.lat > 65
NPI[,drop] <- NA

drop <- slp.lon < 160
NPI[,drop] <- NA

drop <- slp.lon > 220
NPI[,drop] <- NA

# plot to check
NPI.mean <- colMeans(NPI)
z <- t(matrix(NPI.mean,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=tim.colors(64))
contour(slp.x, slp.y, z, add=T) 
map('world2Hires',fill=F,add=T, lwd=2)

# now weight by latitude
weight <- sqrt(cos(slp.lat*pi/180))
ff <- function(x) weighted.mean(x, w=weight, na.rm=T)

NPI <- apply(NPI, 1, ff)
# NPI <- tapply(NPI,slp.yr.NDJFM,mean)

# save and leave out 1948, which is incomplete
salmon.covariates <- data.frame(year=names(NPI[2:length(NPI)]),
                                NPI.NDJFM=NPI[2:length(NPI)])

salmon.covariates$SSTa.NDJFM <- NA
salmon.covariates$SSTa.NDJFM <- 
  sst.anom[match(salmon.covariates$year, names(sst.anom))]

# now FMA 20 m salinity from GAK1

# now GAK1 salinity....
data <- read.csv("/Users/MikeLitzow/Documents/R/climate-data/data/GAK1-6.27.19.csv")

data$year <- floor(data$dec.year)

require(lubridate)
data$day <- yday(date_decimal(data$dec.year)) # get Julian date

sal.set <- as.data.frame(tapply(data$Sal, list(data$dec.year, data$Depth), mean))
colnames(sal.set) <- c("d0", "d10", "d20", "d30", "d50", "d75", "d100", "d150", "d200", "d250")

sal.set$year <- floor(as.numeric(rownames(sal.set))) 
sal.set$day <- yday(date_decimal(as.numeric(rownames(sal.set))))

unloadNamespace("lubridate")

salFMA <- filter(sal.set, day>31 & day <=120)

sal20mu <- tapply(salFMA$d20, salFMA$year, mean, na.rm=T)
plot(names(sal20mu), sal20mu, type="b")

salmon.covariates$FMA.20m.GAK1.sal <- NA
salmon.covariates$FMA.20m.GAK1.sal <- 
  sal20mu[match(salmon.covariates$year, names(sal20mu))]


# papa advection 
papa <- read.csv("/Users/MikeLitzow/Documents/R/NSF-GOA/paper 1 final/submitted/papa.csv")
papa.vector <- papa$latitude
names(papa.vector) <- papa$year

salmon.covariates$Papa <- 
  papa.vector[match(salmon.covariates$year, names(papa.vector))]

# coastal ssh
# need to create synthetic time series combining GODAS and SODA time series!

# first SODA
nc <- nc_open("/Users/MikeLitzow/Documents/R/NSF-GOA/paper 1 final/submitted/hawaii_3e19_7ccd_16ff_e908_64db_63e9.nc")
# view dates (middle of month):
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

ssh <- ncvar_get(nc, "ssh") # get all the data!
x <- ncvar_get(nc, "longitude")     # view longitudes (degrees East)
y <- ncvar_get(nc, "latitude")     # view latitudes
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# process!
ssh <- aperm(ssh, 3:1)  # First, reverse order of dimensions ("transpose" array)

ssh <- matrix(ssh, nrow=dim(ssh)[1], ncol=prod(dim(ssh)[2:3]))  # Change to matrix

dimnames(ssh) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# block out evertyhing but the area we want
ssh.mean <- colMeans(ssh)

drop <- ssh.mean < 0.07
ssh[,drop] <- NA

ssh[,lon < 198] <- NA
ssh[,lat < 55] <- NA

# check
SSH.mean <- colMeans(ssh)
z <- t(matrix(SSH.mean,length(y)))  
image(x,y,z, col=tim.colors(64), xlim=c(180,250), ylim=c(50,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

m <- months(d)
yr <- years(d)

ssh.soda.fma <- ssh.mean.soda[m %in% c("Feb", "Mar", "Apr")]
yr.fma <- yr[m %in% c("Feb", "Mar", "Apr")]
ssh.soda.fma <- tapply(ssh.soda.fma, yr.fma, mean)

# now GODAS
nc <- nc_open("/Users/MikeLitzow/Documents/R/climate-data/data/North.Pacific.godas.sshgsfc")
# view dates (middle of month):
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

ssh <- ncvar_get(nc, "sshgsfc") # get all the data!
x <- ncvar_get(nc, "longitude")     # view longitudes (degrees East)
y <- ncvar_get(nc, "latitude")     # view latitudes
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# process!
ssh <- aperm(ssh, 3:1)  # First, reverse order of dimensions ("transpose" array)

ssh <- matrix(ssh, nrow=dim(ssh)[1], ncol=prod(dim(ssh)[2:3]))  # Change to matrix

dimnames(ssh) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# block out evertyhing but the area we want
ssh.mean <- colMeans(ssh, na.rm=T)
drop <- ssh.mean < 0.07
ssh[,drop] <- NA

ssh[,lon < 198] <- NA
ssh[,lat < 55] <- NA

# check
SSH.mean <- colMeans(ssh, na.rm=T)
z <- t(matrix(SSH.mean,length(y)))  
image.plot(x,y,z, col=tim.colors(64), xlim=c(180,250), ylim=c(50,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# save coordinate and means for study site Fig.
# save x and y for study site fig
x.ssh <- x
y.ssh <- y

plot.ssh <- colMeans(ssh)

ssh.mean.soda <- rowMeans(ssh, na.rm=T)

ssh.mean.godas <- rowMeans(ssh, na.rm=T)

m <- months(d)
yr <- years(d)

ssh.godas.fma <- ssh.mean.godas[m %in% c("Feb", "Mar", "Apr")]
yr.fma <- yr[m %in% c("Feb", "Mar", "Apr")]
ssh.godas.fma <- tapply(ssh.godas.fma, yr.fma, mean)

plot(ssh.godas.fma[names(ssh.godas.fma <= 2010)], ssh.soda.fma[names(ssh.godas.fma >=1980)])
cor(ssh.godas.fma[names(ssh.godas.fma <= 2010)], ssh.soda.fma[names(ssh.godas.fma >=1980)], use="p") # 0.778

# estimate pre-1980 GODAS values from SODA
mod <- lm(ssh.godas.fma[names(ssh.godas.fma <= 2010)] ~ ssh.soda.fma[names(ssh.godas.fma >=1980)], 
          na.action="na.exclude")
godas.estimate.fma <- coef(mod)[2]*ssh.soda.fma + coef(mod)[1]

plot(1949:2010, godas.estimate.fma, type="l", col="red", xlim=c(1949,2019), ylim=c(0.02,0.22))
lines(1980:2019, ssh.godas.fma, col="blue")

godas.synthetic.fma <- c(godas.estimate.fma[names(godas.estimate.fma) <1980], ssh.godas.fma)

plot(1949:2019, godas.synthetic.fma, type="l")

salmon.covariates$FMA.ssh.GODAS <-
  ssh.godas.fma[match(salmon.covariates$year, names(ssh.godas.fma))]

salmon.covariates$FMA.ssh.SODA <- 
  ssh.soda.fma[match(salmon.covariates$year, names(ssh.soda.fma))] 

salmon.covariates$FMA.ssh.COMBINED <- 
  godas.synthetic.fma[match(salmon.covariates$year, names(godas.synthetic.fma))] 


# total GOA FMA wind stress
# as with SSH, we will make both SODA and GODAS time series
# first SODA
nc <- nc_open("/Users/MikeLitzow/Documents/R/NSF-GOA/paper 1 final/submitted/hawaii_3e19_7ccd_16ff_e908_64db_63e9.nc")
# view dates (middle of month):
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

taux <- ncvar_get(nc, "taux") # get all the data!
tauy <- ncvar_get(nc, "tauy")

x <- ncvar_get(nc, "longitude")     # view longitudes (degrees East)
y <- ncvar_get(nc, "latitude")     # view latitudes
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# process!
taux <- aperm(taux, 3:1)  # First, reverse order of dimensions ("transpose" array)
taux <- matrix(taux, nrow=dim(taux)[1], ncol=prod(dim(taux)[2:3]))  # Change to matrix

# same for tauy
tauy <- aperm(tauy, 3:1)  # First, reverse order of dimensions ("transpose" array)
tauy <- matrix(tauy, nrow=dim(tauy)[1], ncol=prod(dim(tauy)[2:3]))

# and convert to total stress
stress <- sqrt(taux^2+tauy^2)
dimnames(stress) <-  list(as.character(d), paste("N", lat, "E", lon, sep=""))

# block out evertyhing but the area we want
# drop a series of points to mark out the AK Peninsula and Bristol Bay

stress[,lat<55] <- NA

px <- c(197, 199, 200.5, 201.5, 203, 203.5, 204)
py <- c(55, 55.5, 56, 56.5, 57, 58, 59)

for(i in 1:length(px)){
  
  stress[,lat > py[i] & lon < px[i]] <- NA
}

# check
z <- t(matrix(colMeans(stress),length(y)))  
image(x,y,z, col=tim.colors(64), xlim=c(180,250), ylim=c(50,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

stress.mean.soda <- rowMeans(stress, na.rm=T)

m <- months(d)
yr <- years(d)

stress.soda.fma <- stress.mean.soda[m %in% c("Feb", "Mar", "Apr")]
yr.fma <- yr[m %in% c("Feb", "Mar", "Apr")]
stress.soda.fma <- tapply(stress.soda.fma, yr.fma, mean)

plot(1949:2010, stress.soda.fma, type="l")

salmon.covariates$FMA.stress.SODA <-
  stress.soda.fma[match(salmon.covariates$year, names(stress.soda.fma))]


# now GODAS
# first v stress
nc <- nc_open("/Users/MikeLitzow/Documents/R/climate-data/data/North.Pacific.godas.vflxsfc")
# view dates (middle of month):
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

v.stress <- ncvar_get(nc, "vflxsfc") # get all the data!
x <- ncvar_get(nc, "longitude")     # view longitudes (degrees East)
y <- ncvar_get(nc, "latitude")     # view latitudes
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# process
v.stress <- aperm(v.stress, 3:1)  
v.stress <- matrix(v.stress, nrow=dim(v.stress)[1], ncol=prod(dim(v.stress)[2:3])) 

nc <- nc_open("/Users/MikeLitzow/Documents/R/climate-data/data/North.Pacific.godas.uflxsfc")
u.stress <- ncvar_get(nc, "uflxsfc")
u.stress <- aperm(u.stress, 3:1)  
u.stress <- matrix(u.stress, nrow=dim(u.stress)[1], ncol=prod(dim(u.stress)[2:3])) 

stress <- sqrt(v.stress^2 + u.stress^2)

dimnames(stress) <-  list(as.character(d), paste("N", lat, "E", lon, sep=""))

# check
z <- t(matrix(colMeans(stress, na.rm=T),length(y)))  
image(x,y,z, col=tim.colors(64), xlim=c(180,250), ylim=c(50,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# block out evertyhing but the area we want
# drop a series of points to mark out the AK Peninsula and Bristol Bay

stress[,lat<55] <- NA
stress[,lon<200] <- NA

# check
z <- t(matrix(colMeans(stress, na.rm=T),length(y)))  
image(x,y,z, col=tim.colors(64), xlim=c(180,250), ylim=c(50,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# save coordinate and means for study site Fig.
# save x and y for study site fig
x.stress <- x
y.stress <- y

plot.stress <- colMeans(stress, na.rm = T)
sum(!is.na(plot.stress)) # just to check - 213 positive cells, which looks right

stress.mean.godas <- rowMeans(stress, na.rm=T)

m <- months(d)
yr <- years(d)

stress.godas.fma <- stress.mean.godas[m %in% c("Feb", "Mar", "Apr")]
yr.fma <- yr[m %in% c("Feb", "Mar", "Apr")]
stress.godas.fma <- tapply(stress.godas.fma, yr.fma, mean)

plot(1980:2019, stress.godas.fma, type="l", col="blue", xlim = c(1949,2019), ylim=c(0.05, 0.14))
lines(1949:2010, stress.soda.fma, col="red")

mod <- lm(stress.godas.fma[names(stress.godas.fma) %in% 1980:2010] ~
            stress.soda.fma[names(stress.soda.fma) %in% 1980:2010])

fma.stress.GODAS.predicted <- coef(mod)[2]*stress.soda.fma + coef(mod)[1]
fma.stress.SYNTHETIC <- c(fma.stress.GODAS.predicted[names(fma.stress.GODAS.predicted) < 1980],
                          stress.godas.fma)

plot(1949:2019, fma.stress.SYNTHETIC, type="l")

salmon.covariates$FMA.stress.GODAS <-
  stress.godas.fma[match(salmon.covariates$year, names(stress.godas.fma))]

salmon.covariates$FMA.stress.COMBINED <-
  fma.stress.SYNTHETIC[match(salmon.covariates$year, names(fma.stress.SYNTHETIC))]

# make a study site figure

png("study site figure.png", 8,5.5, units="in", res=300)
par(mfrow=c(2,2), tcl=-0.2, cex.lab=0.8, cex.axis=0.8, mar=c(1,2,1,0.5),mgp=c(1.5,0.3,0), oma=c(0,0,0,2))
mt.cex <- 1.1
l.mar <- 3
l.cex <- 1.3
l.l <- 0.1
land.col <- "lightyellow3"
xlim <- c(190, 233)
ylim <- c(49,62)
new.col <- oceColorsPalette(64)


z <- t(matrix(plot.stress,length(y.stress)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
zlim <- c(0, max(z, na.rm=T))
image.plot(x.stress,y.stress,z, col=new.col,  yaxt="n", xaxt="n", zlim=zlim, xlim=xlim, ylim=ylim, legend.cex=l.cex, xlab="", ylab="", legend.line = l.l)
map('world2Hires', c('Canada', 'usa'), fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col=land.col)
gak.y <- 59 + (50.7/60)
gak.x <- 360-149 + (28/60)
points(gak.x, gak.y, pch=21, bg="#CC79A7", col="black", cex=1.7)
text(gak.x+4, gak.y-0.2, "GAK1", col="#CC79A7", cex=1.2)
mtext("a) Wind stress (Pa) and GAK1 site", adj=0, cex=1)


z <- t(matrix(plot.ssh,length(y.ssh)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image.plot(x.ssh,y.ssh,z, col=new.col, zlim=c(-0.181, 0.181), yaxt="n", xaxt="n", xlim=xlim, ylim=ylim, xlab="", ylab="")
map('world2Hires', c('Canada', 'usa'), fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col=land.col)
points(360-145,50, pch=21, bg="#56B4E9", col="black", cex=1.5)
text(360-145, 50, "Ocean Station Papa ", pos=2, cex=1.2, col="#56B4E9")
mtext("b) Sea surface height (m) and Ocean Station Papa", adj=0, cex=1)


z <- t(matrix(colMeans(SST), length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image.plot(sst.x,sst.y,z, col=new.col, yaxt="n", xaxt="n", xlim=xlim, ylim=ylim, zlim=c(6,11), xlab="", ylab="")
box()
map('world2Hires', c('Canada', 'usa'), fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col=land.col)

mtext("c) Sea surface temperature (°C)", adj=0, cex=1)

# location...
# study site box
x1 <- c(190,190,233,233,190)
y1 <- c(49,62,62,49,49)

# NPI box
# calculate NDJFM NPI - 
# 30º-65ºN, 160ºE-140ºW 
npi.x <- c(160-1.25,160-1.25,221.25,221.25,160-1.25)
npi.y <- c(30-1.25,66.25,66.25,30-1.25,30-1.25)
par(mar=c(1,2,1,4))
plot(10,10, type="n", xlim=c(130,250), ylim=c(20,66), xlab="", ylab="", xaxt="n", yaxt="n")
map('world2Hires', c('Canada', 'usa', 'Mexico', 'hawaii', 'ussr', 'China', 'Japan', 'South Korea', 'North Korea'), fill=T, add=T, lwd=0.5, col=land.col)

# map('world2Hires',fill=F,add=T, lwd=1)
lines(x1,y1, lwd=2, col="#D55E00")
lines(npi.x,npi.y, lwd=2, col="#0072B2")

text(130, 24, "Study site", cex=1.2, col="#D55E00", pos=4)
text(130, 20.5, "North Pacific Index", cex=1.2, col="#0072B2", pos=4)
mtext("d) Study site and North Pacific Index", adj=0, cex=1)
dev.off()
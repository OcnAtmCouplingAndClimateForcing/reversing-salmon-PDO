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

# restrict to 1964:2019
sst.anom <- sst.anom[names(sst.anom) %in% 1964:2019]

# plot to check
plot(1964:2019, sst.anom, type="l") # looks right... but WTF!
lines(1964:2019, rollmean(sst.anom, 5, fill=NA), col="red", lwd=1.5)
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

# now separate into three months means for each cell - 
# NDJ, FMA, and MJJ
slp.m <- months(slp.d)
slp.yr <- as.numeric(as.character(years(slp.d)))
slp.win.yr <- ifelse(slp.m %in% c("Nov", "Dec"), slp.yr+1, slp.yr)

slp.NDJ <- SLP[slp.m %in% c("Nov", "Dec", "Jan"),]
slp.yr.NDJ <- slp.win.yr[slp.m %in% c("Nov", "Dec", "Jan")]
f <- function(x) tapply(x, slp.yr.NDJ, mean)
slp.NDJ <- apply(slp.NDJ, 2, f)

slp.FMA <- SLP[slp.m %in% c("Feb", "Mar", "Apr"),]
slp.yr.FMA <- slp.yr[slp.m %in% c("Feb", "Mar", "Apr")]
f <- function(x) tapply(x, slp.yr.FMA, mean)
slp.FMA <- apply(slp.FMA, 2, f)

slp.MJJ <- SLP[slp.m %in% c("May", "Jun", "Jul"),]
slp.yr.MJJ <- slp.yr[slp.m %in% c("May", "Jun", "Jul")]
f <- function(x) tapply(x, slp.yr.MJJ, mean)
slp.MJJ <- apply(slp.MJJ, 2, f)

# plot to check
SLP.mean <- colMeans(slp.MJJ)
z <- t(matrix(SLP.mean,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=tim.colors(64))
contour(slp.x, slp.y, z, add=T) 
map('world2Hires',fill=F,add=T, lwd=2)

# plot to check
xx <- c(191.25, 218.75, 218.75, 191.25, 191.25)
yy <- c(58.75, 58.75, 51.25, 51.25, 58.75)

SLP.cv <- apply(slp.MJJ, 2, sd)/SLP.mean
z <- t(matrix(SLP.cv,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=tim.colors(64), main="MJJ SLP Coefficient of variation")
contour(slp.x, slp.y, z, add=T) 
map('world2Hires',fill=F,add=T, lwd=2)
lines(xx, yy, col="white")

# that looks like a reasonable area to average summer slp over summertime

# compare with FMA values
# change box
xx <- c(183.75, 211.25, 211.25, 183.75, 183.75)
yy <- c(58.75, 58.75, 43.75, 43.75, 58.75)

SLP.cv <- apply(slp.FMA, 2, sd)/SLP.mean
z <- t(matrix(SLP.cv,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=tim.colors(64), main="FMA SLP Coefficient of variation")
contour(slp.x, slp.y, z, add=T) 
map('world2Hires',fill=F,add=T, lwd=2)
lines(xx, yy, col="white")


# summarize for regression
spring.slp <- slp.FMA
drop <- slp.lat < min(yy) # using values defined for the box in the plot 
spring.slp[,drop] <- NA

drop <- slp.lat > max(yy)
spring.slp[,drop] <- NA

drop <- slp.lon > max(xx) 
spring.slp[,drop] <- NA  

drop <- slp.lon < min(xx)
spring.slp[,drop] <- NA

# check
SLP.mean <- colMeans(spring.slp)
z <- t(matrix(SLP.mean,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=tim.colors(64))
contour(slp.x, slp.y, z, add=T) 
map('world2Hires',fill=F,add=T, lwd=2)

# and calculate area-weighted mean to regress on PDO
# now weight by latitude
weight <- sqrt(cos(slp.lat*pi/180))
ff <- function(x) weighted.mean(x, w=weight, na.rm=T)

spring.slp <- apply(spring.slp, 1, ff)
plot(names(spring.slp), spring.slp, type="l")

#######################
# summarize MJJ data for regression
summer.slp <- slp.MJJ
drop <- slp.lat < 51.25 
summer.slp[,drop] <- NA
  
drop <- slp.lat > 58.75
summer.slp[,drop] <- NA
  
drop <- slp.lon > 218.75 
summer.slp[,drop] <- NA  

drop <- slp.lon < 191.25
summer.slp[,drop] <- NA

# check
SLP.mean <- colMeans(summer.slp)
z <- t(matrix(SLP.mean,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=tim.colors(64))
contour(slp.x, slp.y, z, add=T) 
map('world2Hires',fill=F,add=T, lwd=2)

# and calculate area-weighted mean to regress on PDO
# now weight by latitude
weight <- sqrt(cos(slp.lat*pi/180))
ff <- function(x) weighted.mean(x, w=weight, na.rm=T)

summer.slp <- apply(summer.slp, 1, ff)
plot(names(summer.slp), summer.slp, type="l")

# plot to check
SLP.mean <- colMeans(slp.MJJ[rownames(slp.MJJ) %in% 2014:2019,])
z <- t(matrix(SLP.mean,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=tim.colors(64))
contour(slp.x, slp.y, z, add=T) 
map('world2Hires',fill=F,add=T, lwd=2)

# examine winter patterns
# plot to check
SLP.mean <- colMeans(slp.NDJ)
z <- t(matrix(SLP.mean,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=tim.colors(64))
contour(slp.x, slp.y, z, add=T) 
map('world2Hires',fill=F,add=T, lwd=2)

SLP.mean <- colMeans(slp.NDJ[rownames(slp.NDJ) %in% 2014:2019,])
z <- t(matrix(SLP.mean,length(slp.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(slp.x,slp.y,z, col=tim.colors(64))
contour(slp.x, slp.y, z, add=T) 
map('world2Hires',fill=F,add=T, lwd=2)
r1 <- r2 <- r3 <- NA # vectors to catch regression coefficients

# window start years for all 6-yr windows that can be fit into 1950-2019
# win.start <- 1950:2014
win.start <- seq(1954,2014,by=6)

# now we'll go through and regress the three SLP periods
# on winter SST for each of the 6-yr windows and plot!

# do the regressions first and get the range before plotting!

NDJ.regr <- FMA.regr <- MJJ.regr <- matrix(nrow = ncol(SLP), ncol = length(win.start))

for(i in 1:length(win.start)){
  # i <- 1
  window <- win.start[i]:(win.start[i]+5)
  x <- sst.anom[names(sst.anom) %in% window]
  
  # loop through each cell
  for(j in 1:ncol(SLP)){
   # j <- 1
    mod <- lm(slp.NDJ[rownames(slp.NDJ) %in% window,j]~x)
    NDJ.regr[j,i] <- coef(mod)[2]
    
    mod <- lm(slp.FMA[rownames(slp.FMA) %in% window,j]~x)
    FMA.regr[j,i] <- coef(mod)[2]
    
    mod <- lm(slp.MJJ[rownames(slp.MJJ) %in% window,j]~x)
    MJJ.regr[j,i] <- coef(mod)[2]
  }}
  
lim <- range(NDJ.regr, FMA.regr, MJJ.regr) # range limits for plotting



par(mar=c(0.5,0.5,1,0.5), mfrow=c(6,3))

new.col <- tim.colors(64)
grays <- c("gray90", "gray89", "gray88",
           "gray87","gray86")

new.col[27:36] <- c(grays[5:1], grays[1:5])

for(j in 1:ncol(NDJ.regr)){
 # j <- 1
  z <- t(matrix(NDJ.regr[,j],length(slp.y)))
  image(slp.x,slp.y,z, col=new.col, zlim=c(lim[1], -lim[1]), xlim=c(150,240), ylim=c(40,70),
        yaxt="n", xaxt="n", ylab="", xlab="")
  contour(slp.x, slp.y, z, add=T, col="grey") 
  map('world2Hires',fill=F,add=T, lwd=2)
  mtext(paste("NDJ SLP ", win.start[j], "-", win.start[j]+5, sep=""), cex=0.8)
  
  z <- t(matrix(FMA.regr[,j],length(slp.y)))
  image(slp.x,slp.y,z, col=new.col, zlim=c(lim[1], -lim[1]), xlim=c(150,240), ylim=c(40,70),
        yaxt="n", xaxt="n", ylab="", xlab="")
  contour(slp.x, slp.y, z, add=T, col="grey") 
  map('world2Hires',fill=F,add=T, lwd=2)
  mtext(paste("FMA SLP ", win.start[j], "-", win.start[j]+5, sep=""), cex=0.8) 
  
  z <- t(matrix(MJJ.regr[,j],length(slp.y)))
  image(slp.x,slp.y,z, col=new.col, zlim=c(lim[1], -lim[1]), xlim=c(150,240), ylim=c(40,70),
        yaxt="n", xaxt="n", ylab="", xlab="")
  contour(slp.x, slp.y, z, add=T, col="grey") 
  map('world2Hires',fill=F,add=T, lwd=2)
  mtext(paste("MJJ SLP ", win.start[j], "-", win.start[j]+5, sep=""), cex=0.8)
}

#####
# the 2nd and 3rd periods are tough to interpret
# look at just NDJFM SLP to simplify and see if this idea holds merit

slp.NDJFM <- SLP[slp.m %in% c("Nov", "Dec", "Jan", "Feb", "Mar"),]
slp.yr.NDJFM <- slp.win.yr[slp.m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
f <- function(x) tapply(x, slp.yr.NDJFM, mean)
slp.NDJFM <- apply(slp.NDJFM, 2, f)

# recalculate / replot
win.start <- seq(1954,2014,by=6)

NDJFM.regr <- matrix(nrow = ncol(SLP), ncol = length(win.start))

for(i in 1:length(win.start)){
  # i <- 1
  window <- win.start[i]:(win.start[i]+5)
  x <- sst.anom[names(sst.anom) %in% window]
  
  # loop through each cell
  for(j in 1:ncol(SLP)){
    # j <- 1
    mod <- lm(slp.NDJFM[rownames(slp.NDJFM) %in% window,j]~x)
    NDJFM.regr[j,i] <- coef(mod)[2]

  }}

lim <- range(NDJFM.regr) # range limits for plotting



par(mar=c(0.5,0.5,1,0.5), mfrow=c(4,3))

new.col <- tim.colors(64)
grays <- c("gray90", "gray89", "gray88",
           "gray87","gray86")

new.col[27:36] <- c(grays[5:1], grays[1:5])

for(j in 1:ncol(NDJFM.regr)){
  # j <- 1
  z <- t(matrix(NDJFM.regr[,j],length(slp.y)))
  image(slp.x,slp.y,z, col=new.col, zlim=c(lim[1], -lim[1]), xlim=c(150,240), ylim=c(40,70),
        yaxt="n", xaxt="n", ylab="", xlab="")
  contour(slp.x, slp.y, z, add=T, col="grey") 
  map('world2Hires',fill=F,add=T, lwd=2)
  mtext(paste("NDJFM SLP ", win.start[j], "-", win.start[j]+5, sep=""), cex=0.8)
  
}
dev.off()

##
# last thing for now - compare 1950-1988, 1989-2013, and 2014-2019

# recalculate / replot
win1 <- 1950:1988
win2 <- 1989:2013
win3 <- 2014:2019

era.regr <- matrix(nrow = ncol(SLP), ncol = 3)


# loop through each cell
for(j in 1:ncol(SLP)){

  x <- sst.anom[names(sst.anom) %in% win1]
  mod <- lm(slp.NDJFM[rownames(slp.NDJFM) %in% win1,j]~x)
  era.regr[j,1] <- coef(mod)[2]
  
  x <- sst.anom[names(sst.anom) %in% win2]
  mod <- lm(slp.NDJFM[rownames(slp.NDJFM) %in% win2,j]~x)
  era.regr[j,2] <- coef(mod)[2]
  
  x <- sst.anom[names(sst.anom) %in% win3]
  mod <- lm(slp.NDJFM[rownames(slp.NDJFM) %in% win3,j]~x)
  era.regr[j,3] <- coef(mod)[2]
    
  }

lim <- range(era.regr) # range limits for plotting


par(mar=c(0.5,0.5,1,0.5), mfrow=c(1,3))

new.col <- tim.colors(64)
grays <- c("gray90", "gray89", "gray88",
           "gray87","gray86")

new.col[27:36] <- c(grays[5:1], grays[1:5])

  z <- t(matrix(era.regr[,1],length(slp.y)))
  image(slp.x,slp.y,z, col=new.col, zlim=c(lim[1], -lim[1]), xlim=c(150,240), ylim=c(40,70),
        yaxt="n", xaxt="n", ylab="", xlab="")
  contour(slp.x, slp.y, z, add=T, col="grey") 
  map('world2Hires',fill=F,add=T, lwd=2)
  mtext("NDJFM SLP 1950-1988", cex=0.8)
  
  z <- t(matrix(era.regr[,2],length(slp.y)))
  image(slp.x,slp.y,z, col=new.col, zlim=c(lim[1], -lim[1]), xlim=c(150,240), ylim=c(40,70),
        yaxt="n", xaxt="n", ylab="", xlab="")
  contour(slp.x, slp.y, z, add=T, col="grey") 
  map('world2Hires',fill=F,add=T, lwd=2)
  mtext("NDJFM SLP 1989-2013", cex=0.8)
  
  z <- t(matrix(era.regr[,3],length(slp.y)))
  image(slp.x,slp.y,z, col=new.col, zlim=c(lim[1], -lim[1]), xlim=c(150,240), ylim=c(40,70),
        yaxt="n", xaxt="n", ylab="", xlab="")
  contour(slp.x, slp.y, z, add=T, col="grey") 
  map('world2Hires',fill=F,add=T, lwd=2)
  mtext("NDJFM SLP 2014-2019", cex=0.8)



#######
# I think the better thing to do, rather than trying to compare
# a slp-sst regression based on only 6 points, 
# would be to plot the time series of SLPa values from the area
# (47–57 N, 153–169 W)


slp.m <- months(slp.d)

f <- function(x) tapply(x, slp.m, mean)
mu <- apply(SLP, 2, f)

mu <- mu[rep(1:12, floor(length(slp.d)/12)),] 

xtra <- 12*((length(slp.d)/12)-floor(length(slp.d)/12))

mu <- rbind(mu, mu[1:xtra,])

slp.anom <- SLP - mu

# now restrict to our "Aleutian Low" area
keep.lat <- slp.lat >= 47 & slp.lat <= 56
keep.lon <- slp.lon >= 192 & slp.lon <= 206
slp.check <- SLP[,keep.lat & keep.lon]
plot.y <- slp.y[slp.y >= 47 & slp.y <= 56]
plot.x <- slp.x[slp.x >= 192 & slp.x <= 206]

slp.anom <- slp.anom[,keep.lat & keep.lon]

# check to plot!
# plot to check
z <- t(matrix(colMeans(slp.check),length(plot.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(plot.x,plot.y,z, col=tim.colors(64), xlim=c(130,250), ylim=c(20,70))
contour(plot.x, plot.y, z, add=T) 
map('world2Hires',fill=F,add=T, lwd=2) # looks good!

# anom.NDJFM <- slp.anom[slp.m %in% c("Nov", "Dec", "Jan", "Feb", "Mar"),]
# slp.yr.NDJFM <- slp.win.yr[slp.m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
# f <- function(x) tapply(x, slp.yr.NDJFM, mean)
# anom.NDJFM <- apply(anom.NDJFM, 2, f)
# anom.NDJFM <- rowMeans(anom.NDJFM)

# naming as AL (=="Aleutian Low" area)

AL.slp.NDJFM <- slp.check[slp.m %in% c("Nov", "Dec", "Jan", "Feb", "Mar"),]
slp.yr.NDJFM <- slp.win.yr[slp.m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
f <- function(x) tapply(x, slp.yr.NDJFM, mean)
AL.slp.NDJFM <- apply(AL.slp.NDJFM, 2, f)
AL.slp.NDJFM <- rowMeans(AL.slp.NDJFM)

# drop pre-64
#  anom.NDJFM <- anom.NDJFM[names(anom.NDJFM) %in% 1965:2019]
AL.slp.NDJFM <- AL.slp.NDJFM[names(AL.slp.NDJFM) %in% 1964:2019]

# plot(1965:2019, scale(anom.NDJFM), type="l", col="blue", ylim=c(-2,3))
plot(1965:2019, scale(slp.NDJFM), type="l", col="blue", ylim=c(-2,3))
lines(1965:2019, scale(sst.anom), col="red")

plot(1965:2019, slp.NDJFM, type="l", col="blue")

cor(anom.NDJFM[names(anom.NDJFM) %in% 1965:1988], sst.anom[names(sst.anom) %in% 1965:1988]) # -0.6184721
cor(anom.NDJFM[names(anom.NDJFM) %in% 1989:2013], sst.anom[names(sst.anom) %in% 1989:2013]) #  -0.6187042!
cor(anom.NDJFM[names(anom.NDJFM) %in% 2014:2019], sst.anom[names(sst.anom) %in% 2014:2019]) # -0.5708142

cor(slp.NDJFM[names(slp.NDJFM) %in% 1965:1988], sst.anom[names(sst.anom) %in% 1965:1988]) # -0.6184721
mod <- lm(slp.NDJFM[names(slp.NDJFM) %in% 1965:1988] ~ sst.anom[names(sst.anom) %in% 1965:1988])
coef(mod) # -5.879752

cor(slp.NDJFM[names(slp.NDJFM) %in% 1989:2013], sst.anom[names(sst.anom) %in% 1989:2013]) #  -0.6187042!
mod <- lm(slp.NDJFM[names(slp.NDJFM) %in% 1989:2013] ~ sst.anom[names(sst.anom) %in% 1989:2013])
coef(mod) # -5.0201

cor(slp.NDJFM[names(slp.NDJFM) %in% 2014:2019], sst.anom[names(sst.anom) %in% 2014:2019]) # -0.5708142
mod <- lm(slp.NDJFM[names(slp.NDJFM) %in% 2014:2019] ~ sst.anom[names(sst.anom) %in% 2014:2019])
coef(mod) # -7.395536

sd(slp.NDJFM[names(slp.NDJFM) %in% 1965:1988]) # 5.12285
sd(slp.NDJFM[names(slp.NDJFM) %in% 1989:2013]) # 3.847468
sd(slp.NDJFM[names(slp.NDJFM) %in% 2014:2019]) # 5.92545

sd(sst.anom[names(sst.anom) %in% 1965:1988]) # 0.538856
sd(sst.anom[names(sst.anom) %in% 1989:2013]) # 0.4741827
sd(sst.anom[names(sst.anom) %in% 2014:2019]) # 0.4573477

# starting to feel my way into this...
# no difference in the correlations between sst and slp among eras,
# some indication of chaning slope?

# use NDJFM SLP in paper, 
# and also latitude and longitude of peak minimum SLP...

# recalculate SLPa
slp.anom <- SLP - mu

# restict to NDJFM averages
anom.NDJFM <- slp.anom[slp.m %in% c("Nov", "Dec", "Jan", "Feb", "Mar"),]
anom.yr.NDJFM <- slp.win.yr[slp.m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
f <- function(x) tapply(x, anom.yr.NDJFM, mean)
anom.NDJFM <- apply(anom.NDJFM, 2, f)

# plot all 72 years to check!

pdf("annual NP winter SLP anomalies.pdf", 11,8)

par(mar=c(0.5,0.5,1,0.5), mfrow=c(4,3))
lim <- range(anom.NDJFM)
new.col <- tim.colors(64)
grays <- c("gray90", "gray89", "gray88",
           "gray87","gray86")

new.col[27:36] <- c(grays[5:1], grays[1:5])

for(i in 1:nrow(anom.NDJFM)){
  # j <- 1
  z <- t(matrix(anom.NDJFM[i,],length(slp.y)))
  image(slp.x,slp.y,z, col=new.col, zlim=c(lim[1], -lim[1]), xlim=c(150,240), ylim=c(40,70),
        yaxt="n", xaxt="n", ylab="", xlab="")
  contour(slp.x, slp.y, z, add=T, col="grey") 
  map('world2Hires',fill=F,add=T, lwd=2)
  mtext(paste("NDJFM SLPa ", rownames(anom.NDJFM)[i], sep=""), cex=0.8)
  
}
dev.off()

# so SLPa isn't what we want...

# try SLP!

# calculate NDJFM mean SLP for the whole area

slp.NDJFM <- SLP[slp.m %in% c("Nov", "Dec", "Jan", "Feb", "Mar"),]
slp.yr.NDJFM <- slp.win.yr[slp.m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
f <- function(x) tapply(x, slp.yr.NDJFM, mean)
slp.NDJFM <- apply(slp.NDJFM, 2, f)

# pdf("annual NP winter SLP values.pdf", 11,8)
# 
# par(mar=c(0.5,0.5,1,0.5), mfrow=c(4,3))
# lim <- range(slp.NDJFM)
# 
# new.col[27:36] <- c(grays[5:1], grays[1:5])
# 
# for(i in 1:nrow(slp.NDJFM)){
#   # j <- 1
#   z <- t(matrix(slp.NDJFM[i,],length(slp.y)))
#   image(slp.x,slp.y,z, col=tim.colors(64), zlim=lim, xlim=c(150,240), ylim=c(40,70),
#         yaxt="n", xaxt="n", ylab="", xlab="")
#   contour(slp.x, slp.y, z, add=T, col="grey") 
#   map('world2Hires',fill=F,add=T, lwd=2)
#   mtext(paste("NDJFM SLP ", rownames(anom.NDJFM)[i], sep=""), cex=0.8)
#   
# }
# dev.off()

# capture time series of the minimum SLP value and its lat/long for each year!
# only need them for 1964:2019

# salmon.covariates <- data.frame(year=1964:2019, 
#                                 NDJFM.sst.anom=sst.anom,
#                                 NDJFM.AL.slp=AL.slp.NDJFM,
#                                 NDJFM.slp.min=NA,
#                                 NDJFM.slp.min.lat=NA,
#                                 NDJFM.slp.min.lon=NA)
# 
# salmon.covariates$era <- as.factor(ifelse(
#   salmon.covariates$year<=1988, 1, 
#   ifelse(salmon.covariates$year %in% 1989:2013, 2, 3)
# ))
# 
# # loop through each year of interest
# 
# for(i in 1964:2019){
#   # i <- 1964
#   temp <- slp.NDJFM[rownames(slp.NDJFM) == i,]
#   
#   salmon.covariates$NDJFM.slp.min[salmon.covariates$year==i] <- min(temp)
#   salmon.covariates$NDJFM.slp.min.lat[salmon.covariates$year==i] <- slp.lat[which.min(temp)]
#   salmon.covariates$NDJFM.slp.min.lon[salmon.covariates$year==i] <- slp.lon[which.min(temp)]
# }

# and!
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


# saving with NPI - this got lost somehow previously!
# salmon.covariates <- read.csv("salmon.covariates.csv", row.names = 1)
# 
# salmon.covariates$NDJFM.NPI <-
#   NPI[match(salmon.covariates$year, names(NPI))]
# 
# write.csv(salmon.covariates, "salmon.covariates.csv", row.names = F)

# 
# ggplot(salmon.covariates, aes(NDJFM.slp.min.lon, NDJFM.slp.min.lat, color=era)) +
#   theme_linedraw() +
#   geom_point() + 
#   geom_smooth(method="lm", se=F)
# 
# ggplot(salmon.covariates, aes(NDJFM.slp.min, NDJFM.slp.min.lat, color=era)) +
#   theme_linedraw() +
#   geom_point() + 
#   geom_smooth(method="lm", se=F)
# 
# ggplot(salmon.covariates, aes(NDJFM.slp.min, NDJFM.slp.min.lon, color=era)) +
#   theme_linedraw() +
#   geom_point() + 
#   geom_smooth(method="lm", se=F)
# 
# ggplot(salmon.covariates, aes(NDJFM.slp.min, NDJFM.sst.anom, color=era)) +
#   theme_linedraw() +
#   geom_point() + 
#   geom_smooth(method="lm", se=F) # that's interesting - the slope is very similar among
# # eras, but the intercept (sst) is very different
# 
# ggplot(salmon.covariates, aes(NDJFM.slp.min.lat, NDJFM.sst.anom, color=era)) +
#   theme_linedraw() +
#   geom_point() + 
#   geom_smooth(method="lm", se=F) 
# 
# ggplot(salmon.covariates, aes(NDJFM.slp.min.lon, NDJFM.sst.anom, color=era)) +
#   theme_linedraw() +
#   geom_point() + 
#   geom_smooth(method="lm", se=F) 
# 
# ggplot(salmon.covariates, aes(NDJFM.slp.min, NDJFM.AL.slp, color=era)) +
#   theme_linedraw() +
#   geom_point() + 
#   geom_smooth(method="lm", se=F) # so that relationship is constant
# 
# ggplot(salmon.covariates, aes(NDJFM.AL.slp, NDJFM.sst.anom, color=era)) +
#   theme_linedraw() +
#   geom_point() + 
#   geom_smooth(method="lm") 

# continue populating with the GOA physical vars we've used before

# first, Ketchikan-Seward SLP difference
# Extract GOA SLP, 47.5-65 deg. N, 190-235 deg. E:

ks.diff <- SLP[, "N55E230"] - SLP[, "N60E210"] # these coordinates are correct!

# now limit to NDJ
ks.diff <- ks.diff[slp.m %in% c("Nov", "Dec", "Jan")]
yr <- years(d)

# advance yr to match Jan
temp.slp.yr <- slp.yr
temp.slp.yr[m %in% c("Nov", "Dec")] <- 1 + temp.slp.yr[m %in% c("Nov", "Dec")]
temp.slp.yr <- temp.slp.yr[m %in% c("Nov", "Dec", "Jan")]

ks.diff <- tapply(ks.diff, temp.slp.yr, mean)
plot(names(ks.diff), ks.diff, type="l")
salmon.covariates$NDJ.slp.grad=ks.diff[names(ks.diff) %in% 1964:2019]

ggplot(salmon.covariates, aes(NDJFM.slp.min, NDJ.slp.grad, color=era)) +
  theme_linedraw() +
  geom_point() + 
  geom_smooth(method="lm") 

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

ggplot(salmon.covariates, aes(NDJFM.AL.slp, FMA.20m.GAK1.sal, color=era)) +
  theme_linedraw() +
  geom_point() + 
  geom_smooth(method="lm", se=F) 

ggplot(salmon.covariates, aes(NDJFM.sst.anom, FMA.20m.GAK1.sal, color=era)) +
  theme_linedraw() +
  geom_point() + 
  geom_smooth(method="lm", se=F) 

# papa advection 
papa <- read.csv("/Users/MikeLitzow/Documents/R/NSF-GOA/paper 1 final/submitted/papa.csv")
papa.vector <- papa$latitude
names(papa.vector) <- papa$year

salmon.covariates$Papa <- 
  papa.vector[match(salmon.covariates$year, names(papa.vector))]

# coastal ssh
# need to create synthetic time series combining GODAS and SODA time series!

# first SODA
nc <- nc_open("/Users/MikeLitzow 1/Documents/R/NSF-GOA/paper 1 final/submitted/hawaii_3e19_7ccd_16ff_e908_64db_63e9.nc")
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

ssh.mean.soda <- rowMeans(ssh, na.rm=T)

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

# add freshwater and river output time series...
hill <- read.csv("/Users/MikeLitzow 1/Documents/R/NSF-GOA/paper 1 final/submitted/hill.original.corrected.dates.csv") 

# compute Hill's anomalies
mu <- tapply(hill$input.km.3, hill$month, mean)
sd <- tapply(hill$input.km.3, hill$month, sd)

mu <- rep(mu, length.out=nrow(hill))
sd <- rep(sd, length.out=nrow(hill))

hill$total.anomaly <- (hill$input.km.3-mu)/sd

# limit to FMA and FMA
fma.hill <- hill[hill$month %in% 2:4,]
fma.fw <- tapply(fma.hill$total.anomaly, fma.hill$year, mean)

salmon.covariates$FMA.FW.HILL <-
  fma.fw[match(salmon.covariates$year, names(fma.fw))]

# Susitna river
su <- read.csv("/Users/MikeLitzow 1/Documents/R/KSWG/susitna r at gold creek.csv") 
su$log.flow <- log(su$mean_va)

su <- su %>%
  filter(month_nu %in% 2:4)


ff <- function(x) sum(!is.na(x))

su_n <- tapply(su$log.flow, su$year_nu, ff)

# using monthly anomalies for now b/c we have incomplete # of obs for 2017
# can quickly change back to raw observations if needed when we have complete data

mu <- tapply(su$log.flow, su$month_nu, mean)
su$mu <- c(rep(mu, (length(unique(su$year_nu))-1)), mu[1])
su$anom <- su$log.flow-su$mu

su_anom <- tapply(su$anom, su$year_nu, mean)
plot(names(su_anom), su_anom, type='b')

salmon.covariates$FMA.Susitna.Flow <- 
  su_anom[match(salmon.covariates$year, names(su_anom))]


# total GOA FMA wind stress
# as with SSH, we will make both SODA and GODAS time series
# first SODA
nc <- nc_open("/Users/MikeLitzow 1/Documents/R/NSF-GOA/paper 1 final/submitted/hawaii_3e19_7ccd_16ff_e908_64db_63e9.nc")
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
nc <- nc_open("/Users/MikeLitzow 1/Documents/R/climate-data/data/North.Pacific.godas.vflxsfc")
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

nc <- nc_open("/Users/MikeLitzow 1/Documents/R/climate-data/data/North.Pacific.godas.uflxsfc")
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

# finally, summer upwelling
ui <- read.csv("/Users/MikeLitzow 1/Documents/R/climate-data/data/upwelling.csv")

# change months to numbers for convenience
colnames(ui)[3:14] <- 1:12

# use the 3 northernmost positions
ui1 <- subset(ui, POSITION=="60N.149W")
ui2 <- subset(ui, POSITION=="60N.146W")
ui3 <- subset(ui, POSITION=="57N.137W")

# stack so that we can turn into a chronology
upwell <- stack(ui1, select=colnames(ui1)[3:14])
upwell$year <- rep(1946:2019, times=12)
colnames(upwell)[1:2] <- c("N60.W149", "month")

ui2 <- stack(ui2, select=colnames(ui2)[3:14])
ui3 <- stack(ui3, select=colnames(ui3)[3:14])

# add the other locations to "upwell"
upwell$N60.W146 <- ui2$values
upwell$N57.W137 <- ui3$values
# and order
upwell <- upwell[order(upwell$year),]

upwell$mean <- rowMeans(upwell[,c(1,4,5)])

MJJ.up <- upwell[upwell$month %in% 5:7,]

MJJ.uw <- tapply(MJJ.up$mean, MJJ.up$year, mean)

plot(1946:2019, MJJ.uw, type="l")

salmon.covariates$MJJ.upwelling <- 
  MJJ.uw[match(salmon.covariates$year, names(MJJ.uw))]

# and....add PDO!
names <- read.table("/Users/MikeLitzow 1/Documents/R/FATE2/non-som/~pdo", skip=30, nrows=1, as.is = T)
pdo <- read.table("/Users/MikeLitzow 1/Documents/R/FATE2/non-som/~pdo", skip=31, nrows=119, fill=T, col.names = names)
pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
pdo <- pdo %>%
  gather(month, value, -YEAR) %>%
  arrange(YEAR)

# load estimated values!
estimated.pdo <- read.csv("estimated PDO through May 2019.csv", row.names = 1)

# and plug these into the pdo we're using (for now!)
more.rows <- data.frame(YEAR=2019, month=c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL"), value=NA)
pdo <-rbind(pdo, more.rows)

pdo$value <- estimated.pdo$real.pdo

win.pdo <- pdo[pdo$month %in% c("NOV", "DEC", "JAN", "FEB", "MAR"),]
win.pdo$win.yr <- ifelse(win.pdo$month %in% c("NOV", "DEC"), win.pdo$YEAR+1, win.pdo$YEAR)
win.pdo <- tapply(win.pdo$value, win.pdo$win.yr, mean)


# use this if salmon covariates not loaded!
salmon.covariates <- read.csv("salmon.covariates.csv", row.names = 1)

salmon.covariates$NDJFM.PDO <-
  win.pdo[match(salmon.covariates$year, names(win.pdo))]

# and add the summer slp!
salmon.covariates$MJJ.SLP.EBS.GOA <- 
  summer.slp[match(salmon.covariates$year, names(summer.slp))]

# and the spring slp
salmon.covariates$FMA.SLP.EBS.GOA <- 
  spring.slp[match(salmon.covariates$year, names(spring.slp))]

# now save for use in the DFA model!
write.csv(salmon.covariates, "salmon.covariates.csv")

library(ncdf4)
library(zoo)
library(gplots)
library(dplyr)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(ggplot2)
library(oce)


# script for calculating GOA sst anomalies 
# download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1950-01-01):1:(2020-12-01T00:00:00Z)][(0.0):1:(0.0)][(54):1:(62)][(200):1:(226)]", "~temp")

# paste into browser for windows!


# load and process SST data
# nc <- nc_open("~temp")

nc <- nc_open("./data/nceiErsstv5_b864_8959_89dd.nc")

# extract dates

ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

# extract study area
# 54-62 deg. N, 200-226 deg. E
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

# drop Bristol Bay cells
BB <- c("N58E200", "N58E202", "N56E200")
SST[,BB] <- NA

# drop eastern and southern areas
SST[,lat == 54] <- NA
SST[,lon > 210] <- NA

# and check
temp.mean <- colMeans(SST, na.rm=T)
z <- t(matrix(temp.mean,length(y)))  
image.plot(x,y,z, col=oceColorsPalette(64), xlim=c(195,215), ylim=c(54,62))
map('world2Hires',c('Canada', 'usa'), fill=T,xlim=c(130,250), ylim=c(20,66),add=T, lwd=1, col="lightyellow3")

# calculate monthly anomaly

SST <- rowMeans(SST, na.rm = T)

# and Jan-Jun annual mean
yr <- as.numeric(as.character(years(d)))
m <- months(d)
jan.jun.m <- m[m %in% c("Jan", "Feb", "Mar", "Apr", "May", "Jun")]
jan.jun.yr <- yr[m %in% c("Jan", "Feb", "Mar", "Apr", "May", "Jun")]

jan.jun.SST <- SST[m %in% c("Jan", "Feb", "Mar", "Apr", "May", "Jun")]
jan.jun.SST <- tapply(jan.jun.SST, jan.jun.yr, mean)

anom <- (jan.jun.SST - mean(jan.jun.SST[names(jan.jun.SST) %in% 1981:2010])) / 
  sd(jan.jun.SST[names(jan.jun.SST) %in% 1981:2010])

xprt <- data.frame(year = 1950:2020,
                   jan.jun.SST = jan.jun.SST,
                   anom = anom)

ggplot(xprt, aes(year, jan.jun.SST)) +
  geom_line() +
  geom_point()

ggplot(xprt, aes(year, anom)) +
  geom_line() +
  geom_point()

write.csv(xprt, "./data/goa.jan.jun.sst.csv", row.names = F)

## winter (NDJFM) means-----------
# and Jan-Jun annual mean
yr <- as.numeric(as.character(years(d)))
m <- months(d)

win.yr <- yr
win.yr <- if_else(m %in% c("Nov", "Dec"), yr+1, yr)  

ndjfm.m <- m[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
ndjfm.yr <- win.yr[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]

ndjfm.SST <- SST[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
ndjfm.SST <- tapply(ndjfm.SST, ndjfm.yr, mean)

# export
xprt <- data.frame(year = 1951:2020,
                   ndjfm.SST = ndjfm.SST[names(ndjfm.SST) %in% 1951:2020])

write.csv(xprt, "./data/goa.ndjfm.sst.csv", row.names = F)

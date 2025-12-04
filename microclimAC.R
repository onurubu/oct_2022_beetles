##### Coskun Kucukkaragoz
#### 22nd of August 2025
### R script exploring a mechanistic model to fill in the Cederberg temp data

# set-up ####

# packages
rm(list = ls())
library(tidyverse)
library(microclima)
library(micropoint)
library(lubridate)

library(micropoint)


# testing
mout <- runpointmodel(climdata, reqhgt = 1, forestparams, paii = NA,  groundparams, lat =49.96807, long= -5.215668)
tme <- as.POSIXct(climdata$obs_time)
plot(mout$tair ~ tme, type="l", xlab = "MOnth", ylab = "Air emperature",
     ylim = c(-5, 35), col = "red")
par(new = TRUE)
plot(climdata$temp ~ tme, type="l", xlab = "", ylab = "",
     ylim = c(-5, 35), col = rgb(0, 0, 0, 0.5))

# Inspect radiation and radiation at time of greatest temperature anomaly
anom <- mout$tair - climdata$temp
hr <- which.max(anom)
par(mfrow = c(2, 1))
hist(climdata$swdown, main = "", xlab = "Downward shortwave radiation")
abline(v = climdata$swdown[hr], col = "red")
hist(climdata$windspeed, main = "", xlab = "Wind speed")
abline(v = climdata$windspeed[hr], col = "red")
# plot vertical temperature profile for that hour
par(mfrow = c(1, 1))
out <- plotprofile(climdata, hr = hr, "tair", forestparams, paii = NA, groundparams, lat =49.96807, long= -5.215668)






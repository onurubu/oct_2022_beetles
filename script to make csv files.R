
rm(list = ls())


xx <- rep(1:17, each = 40)
yy <- rep(rep(1:4, each = 10), 17)
zz <- rep(1:10, 68)

dat <- data.frame(site = xx, replicate = yy, trap = zz)
dat <- rbind(dat, dat)

years <- rep(c(2002,2022), each = 680)
dat <- cbind(dat, years)
dat <- dat %>% rename("year" = "years")

write.csv(dat, "./veg_cover.csv", row.names = FALSE)

# detrend (and deseasonalize) NDVI

library(data.table)
library(mgcv)
library(magrittr)
library(ggplot2)

dat_ndvi <- readRDS("data/gimms-ndvi-01.rds")


f_dets <- function(dat1){
  
  dat_out <- copy(dat1)
  gam1 <- gam(ndvi ~ s(xx_half_month, bs = "cc") + year, data = dat1)
  
  newdat1 <- dat1[, .(year = mean(year), xx_half_month)]
  # fit_detrend <- predict(gam1, newdat1)
  # fit_detrend_deseason <- resid(gam1)
  
  newdat2 <- dat1[, .(year, xx_half_month = mean(xx_half_month))]
  fit_trend <- predict(gam1, newdat2)
  
  dat_out[, ndvi_avg_season := fit_detrend]
  dat_out[!is.na(ndvi), ndvi_resid := resid(gam1)]
  dat_out[, ndvi_detrend := ndvi_avg_season + ndvi_resid]
  dat_out[, ndvi_trend := fit_trend]
  dat_out
}

# interpolate NA's
setkey(dat_ndvi, site, year, xx_half_month)
dat_ndvi[, xx_obs := 1:.N, site]
dat_ndvi[is.na(ndvi)] # -> at most 2 consecutive missing values
dat_ndvi[, ndvi_f := approx(xx_obs, ndvi, xout = xx_obs)$y, site]
dat_ndvi[is.na(ndvi), ndvi := ndvi_f]

dat_ndvi2 <- dat_ndvi[, f_dets(.SD), .(site)]
dat_ndvi2

saveRDS(dat_ndvi2, file = "data/gimms-ndvi-03-detrend.rds")

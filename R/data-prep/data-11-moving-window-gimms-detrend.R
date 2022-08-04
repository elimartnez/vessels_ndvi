# make moving window data



library(magrittr)
library(data.table)
library(foreach)
library(fs)
library(purrr)
library(lubridate)


# gimms - ndvi ------------------------------------------------------------

dat_ndvi <- readRDS("data/gimms-ndvi-03-detrend.rds")

dat_ndvi[, xx_half_month] %>% table
dat_ndvi[, year] %>% table
dat_ndvi[, site] %>% table

dat_ndvi_mov <- copy(dat_ndvi)
setkey(dat_ndvi_mov, site, year, xx_half_month)

# already in detrend
# dat_ndvi_mov[, xx_obs := 1:.N, site]
# # interpolate NA's
# dat_ndvi_mov[is.na(ndvi)] # -> at most 2 consecutive missing values
# dat_ndvi_mov[, ndvi_f := approx(xx_obs, ndvi, xout = xx_obs)$y, site]
# dat_ndvi_mov[is.na(ndvi), ndvi := ndvi_f]

all_xx_half_month <- sort(unique(dat_ndvi_mov$xx_half_month))
all_years <- sort(unique(dat_ndvi_mov$year))
all_periods <- 1:48 - 1 # 2 years

dat_mov <- foreach(
  i_xx_half_month = all_xx_half_month,
  .final = rbindlist
) %do% {

  foreach(
    i_year = all_years,
    .final = rbindlist
  ) %do% {

    foreach(
      i_period = all_periods,
      .final = rbindlist
    ) %do% {

      xx_obs_end <- dat_ndvi_mov[year == i_year & xx_half_month == i_xx_half_month,
                                   unique(xx_obs)]
      xx_obs_start <- xx_obs_end - i_period

      dat_i_mean <- dat_ndvi_mov[xx_obs >= xx_obs_start &
                                   xx_obs <= xx_obs_end,
                                 .(mean_ndvi = mean(ndvi, na.rm = F),
                                   mean_ndvi_detrend = mean(ndvi_detrend, na.rm = F),
                                   mean_ndvi_resid = mean(ndvi_resid, na.rm = F),
                                   nn = .N),
                                 site]

      dat_i_mean[, ":="(xx_half_month = i_xx_half_month, year = i_year, period = i_period)]
      dat_i_mean
    }
  }

}

saveRDS(dat_mov, file = "data/gimms-ndvi-04-detrend-moving.rds")

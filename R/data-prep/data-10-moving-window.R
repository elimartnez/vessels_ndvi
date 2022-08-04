# make moving window data



library(magrittr)
library(data.table)
library(foreach)
library(fs)
library(purrr)
library(lubridate)


# gimms - ndvi ------------------------------------------------------------

dat_ndvi <- readRDS("data/gimms-ndvi-01.rds")

dat_ndvi[, xx_half_month] %>% table
dat_ndvi[, year] %>% table
dat_ndvi[, site] %>% table

dat_ndvi_mov <- copy(dat_ndvi)
setkey(dat_ndvi_mov, site, year, xx_half_month)
dat_ndvi_mov[, xx_obs := 1:.N, site]

# interpolate NA's
dat_ndvi_mov[is.na(ndvi)] # -> at most 2 consecutive missing values
dat_ndvi_mov[, ndvi_f := approx(xx_obs, ndvi, xout = xx_obs)$y, site]
dat_ndvi_mov[is.na(ndvi), ndvi := ndvi_f]

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
                                   nn = .N),
                                 site]

      dat_i_mean[, ":="(xx_half_month = i_xx_half_month, year = i_year, period = i_period)]
      dat_i_mean
    }
  }

}

saveRDS(dat_mov, file = "data/gimms-ndvi-02-moving.rds")



# modis ndvi -------------------------------------------------------------

dat <- readRDS("data/modis-ndvi-01.rds")

dat[, yday(date)] %>% table
dat[, year(date)] %>% table
dat[, site] %>% table

dat_mov <- copy(dat)
setkey(dat_mov, site, date)
dat_mov[, year := year(date)]
dat_mov[, doy := yday(date)]

all_doy <- sort(unique(dat_mov$doy))
all_years <- sort(unique(dat_mov$year))
all_periods <- 1:45 - 1 # 2 years (16day freq)

dat_out <- foreach(
  i_doy = all_doy,
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
      
      end_date <- dat_mov[year == i_year & doy == i_doy, unique(date)]
      
      dat_i_mean <- dat_mov[date >= end_date - i_period * 16 &
                                  date <= end_date,
                                .(mean_ndvi = mean(ndvi, na.rm = T),
                                  mean_ndvi_f = mean(ndvi_f, na.rm = T),
                                  nn = .N),
                                site]
      
      dat_i_mean[, ":="(doy = i_doy, year = i_year, period = i_period)]
      dat_i_mean
    }
  }
  
}

saveRDS(dat_out, file = "data/modis-ndvi-02-moving.rds")



# gpp lai fpar et -------------------------------------------------------------

dat <- readRDS("data/gpp-lai-fpar-et-01.rds")

dat[, yday(date)] %>% table
dat[, year(date)] %>% table
dat[, site] %>% table

dat_mov <- copy(dat)
setkey(dat_mov, site, date)
dat_mov[, year := year(date)]
dat_mov[, doy := yday(date)]

all_doy <- sort(unique(dat_mov$doy))
all_years <- sort(unique(dat_mov$year))
all_periods <- 1:90 - 1 # 2 years

dat_out <- foreach(
  i_doy = all_doy,
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

      end_date <- dat_mov[year == i_year & doy == i_doy, unique(date)]

      dat_i_mean <- dat_mov[date >= end_date - i_period * 8 &
                              date <= end_date,
                            .(mean_lai = mean(lai, na.rm = T),
                              mean_lai_f = mean(lai_f, na.rm = T),
                              mean_fpar = mean(fpar, na.rm = T),
                              mean_fpar_f = mean(fpar_f, na.rm = T),
                              mean_et = mean(et, na.rm = T),
                              mean_et_f = mean(et_f, na.rm = T),
                              mean_gpp = mean(gpp, na.rm = T),
                              mean_gpp_f = mean(gpp_f, na.rm = T),
                              nn = .N),
                            site]

      dat_i_mean[, ":="(doy = i_doy, year = i_year, period = i_period)]
      dat_i_mean
    }
  }

}

saveRDS(dat_out, file = "data/gpp-lai-fpar-et-02-moving.rds")


# EOF ---------------------------------------------------------------------




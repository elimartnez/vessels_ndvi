# read in MODIS data


library(magrittr)
library(data.table)
library(foreach)
library(fs)
library(purrr)
library(lubridate)
source("R/functions/number2binary.R")



# gpp/lai/fpar/et ---------------------------------------------------------------------

# gpp (has same QC as lai and fpar)
dir_ls("data-modis/gpp/") %>% 
  map(fread, colClasses = "character") %>% 
  rbindlist -> dat

dat_gpp <- dat[, .(site, date = ymd(calendar_date), gpp = as.numeric(value), scale = as.numeric(scale))]
dat_gpp[gpp > 30000, gpp := NA]
dat_gpp[, gpp := gpp*scale]
summary(dat_gpp)

# lai
dir_ls("data-modis/lai/") %>% 
  map(fread, colClasses = "character") %>% 
  rbindlist -> dat

dat_lai <- dat[, .(site, date = ymd(calendar_date), lai = as.numeric(value), scale = as.numeric(scale))]
dat_lai[lai > 100, lai := NA]
dat_lai[, lai := lai*scale]
summary(dat_lai)

# fpar
dir_ls("data-modis/fpar/") %>% 
  map(fread, colClasses = "character") %>% 
  rbindlist -> dat

dat_fpar <- dat[, .(site, date = ymd(calendar_date), fpar = as.numeric(value), scale = as.numeric(scale))]
dat_fpar[fpar > 100, fpar := NA]
dat_fpar[, fpar := fpar*scale]
summary(dat_fpar)

# et
dir_ls("data-modis/et/") %>% 
  map(fread, colClasses = "character") %>% 
  rbindlist -> dat

dat_et <- dat[, .(site, date = ymd(calendar_date), et = as.numeric(value), scale = as.numeric(scale))]
dat_et[et > 32700 | et < -32767, table(site)]
dat_et[et > 32700 | et < -32767, et := NA]
dat_et[, et := et*scale]
summary(dat_et)


# merge all three
dat_lai[, .(site, date, lai)] %>% 
  merge(dat_fpar[, .(site, date, fpar)], all.x = T) %>% 
  merge(dat_et[, .(site, date, et)], all.x = T)  %>% 
  merge(dat_gpp[, .(site, date, gpp)], all.x = T) -> dat_vals


# remove ZUR, since NA anyway
dat_vals[site == "ZUR"] %>% summary
dat_vals <- dat_vals[site != "ZUR"]

# add qc layer
# see mod15 user guide table 5 for bit meaning
# dir_ls("data-modis/lai_fpar_et_qc/") %>% 
#   map(fread, colClasses = "character") %>% 
#   rbindlist -> dat_qc
dir_ls("data-modis/gpp_qc/") %>%
  map(fread, colClasses = "character") %>%
  rbindlist -> dat_qc


dat_qc_lai <- dat_qc[, .(site, date = ymd(calendar_date), lai_qc = as.numeric(value))]

bit_cols <- paste0("b_", 7:0)
dat_qc_lai[, (bit_cols) := as.data.table(t(sapply(lai_qc, number2binary, 8)))]

# aggregate bits for cloud and 5 level score
dat_qc_lai[, ":="(b_43 = paste0(b_4, b_3),
                  b_765 = paste0(b_7, b_6, b_5))]

# merge qc to values
dat_out <- merge(
  dat_qc_lai[, .(site, date, b_0, b_2, b_43, b_765)],
  dat_vals,
  by = c("site", "date")
)

# add qc subset
# what to mask out?
#       b_0 == 1 (non-good quality) -> not really, because b_0 == 1 looks actually ok
#       b_43 == 01 (significant clouds)
#       b_43 == 10 (mixed clouds)

dat_out[, ":="(lai_qc = lai,
               fpar_qc = fpar,
               et_qc = et,
               gpp_qc = gpp)]
dat_out[b_43 == "10" | b_43 == "01", 
        ":="(lai_qc = NA,
             fpar_qc = NA,
             et_qc = NA,
             gpp_qc = NA)]

# interpolate missing value
setkey(dat_out, site, date)
dat_out[, xx_obs := 1:.N, site]

# interpolate NA's
dat_out[, lai_f := approx(xx_obs, lai_qc, xout = xx_obs)$y, site]
dat_out[, fpar_f := approx(xx_obs, fpar_qc, xout = xx_obs)$y, site]
dat_out[, et_f := approx(xx_obs, et_qc, xout = xx_obs)$y, site]
dat_out[, gpp_f := approx(xx_obs, gpp_qc, xout = xx_obs)$y, site]



saveRDS(dat_out, file = "data/gpp-lai-fpar-et-01.rds")



# ndvi --------------------------------------------------------------------

dir_ls("data-modis/ndvi/") %>% 
  map(fread, colClasses = "character") %>% 
  rbindlist -> dat

dat_ndvi <- dat[, .(site, date = ymd(calendar_date), ndvi = as.numeric(value), scale = as.numeric(scale))]
dat_ndvi[ndvi < -2000 | ndvi > 10000, ndvi := NA]
dat_ndvi[, ndvi := ndvi*scale]
summary(dat_ndvi)

# add qc layer
dir_ls("data-modis/ndvi_qc/", regexp = "pixel") %>% 
  map(fread, colClasses = "character") %>% 
  rbindlist -> dat_qc

dat_qc_ndvi <- dat_qc[, .(site, date = ymd(calendar_date), ndvi_pixel_rel = as.numeric(value))]

# merge
dat_out <- merge(dat_ndvi[, .(site, date, ndvi)],
                 dat_qc_ndvi)

# add qc subset
# mask out snow and clouds, 2 and 3
# keep 0 good and 1 marginal

dat_out[, ndvi_qc := ndvi]
dat_out[ndvi_pixel_rel > 1, ndvi_qc := NA]

# interpolate missing value
setkey(dat_out, site, date)
dat_out[, xx_obs := 1:.N, site]

# interpolate NA's
dat_out[, ndvi_f := approx(xx_obs, ndvi_qc, xout = xx_obs)$y, site]



saveRDS(dat_out, file = "data/modis-ndvi-01.rds")




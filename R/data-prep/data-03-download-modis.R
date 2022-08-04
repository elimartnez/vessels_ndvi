library(MODISTools)
library(magrittr)
library(data.table)
library(foreach)

# read chrono locations
dat_coord_all <- readRDS("data/chrono-00-meta.rds")
setnames(dat_coord_all, "site", "site_name") 


# some MODIS info ---------------------------------------------------------


dat_prod <- mt_products()
# gpp
# MOD17A2H

# net ET
# MOD16A2

# LAI/FPAR
# MOD15A2H

# NDVI
# MOD13Q1

mt_bands("MOD17A2H")
mt_bands("MOD16A2")
mt_bands("MOD15A2H")
mt_bands("MOD13Q1")



# download MODIS data -----------------------------------------------------

# lai fpar gpp all the same
# ndvi pixel reliability is high level summary

tibble::tribble(
  ~variable, ~product, ~band, 
  # ----------------------- #
  "gpp", "MOD17A2H", "Gpp_500m",
  "et", "MOD16A2", "ET_500m",
  "lai",  "MOD15A2H", "Lai_500m",
  "fpar", "MOD15A2H", "Fpar_500m",
  "ndvi", "MOD13Q1", "250m_16_days_NDVI",
  "ndvi_qc", "MOD13Q1", "250m_16_days_pixel_reliability",
  "gpp_qc", "MOD17A2H", "Psn_QC_500m",
  # "lai_fpar_et_qc",  "MOD15A2H", "FparLai_QC",
  # "ndvi_qc", "MOD13Q1", "250m_16_days_VI_Quality"
) %>% data.table -> dat_modis


dat_modis[, start := "2000-01-01"]
dat_modis[, end := "2019-12-31"]
dat_modis[, internal := F]
dat_modis[, out_dir := paste0("data-modis/", variable, "/")]
dat_modis

# create dirs
fs::dir_create(dat_modis$out_dir)


# download missing
dat_coord_all <- dat_coord_all[c(21:24)]
dat_modis <- dat_modis[c(5)]

# parallel?
# mitmatmisc::init_parallel_ubuntu()

# download all
zz <- foreach(i = 1:nrow(dat_modis)) %do% {
  mt_batch_subset(df = dat_coord_all,
                  product = dat_modis[i, product],
                  band = dat_modis[i, band],
                  start = dat_modis[i, start],
                  end = dat_modis[i, end],
                  internal = dat_modis[i, internal],
                  out_dir = dat_modis[i, out_dir])

}









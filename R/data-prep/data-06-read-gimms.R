# extract the time series from GIMMS

library(gimms)
library(fs)
library(magrittr)
library(sf)
library(foreach)
library(data.table)



# check extent of chronos -------------------------------------------------
dat_coord <- readRDS("data/chrono-00-meta.rds")
sf_coord <- st_as_sf(dat_coord,
                     coords = c("lon", "lat"),
                     crs = "+proj=longlat +datum=WGS84")
mapview::mapview(sf_coord)
ext_custom <- extent(-12, 50, 27, 55) #xmin xmax ymin ymax





# gimms 2 raster ----------------------------------------------------------


gimms_files <- fs::dir_ls("/mnt/CEPH_PROJECTS/CLIRSNOW/gimms/nc/")
tif_files <- path("/mnt/CEPH_PROJECTS/CLIRSNOW/gimms/tif-treestudy/",
                  path_ext_remove(path_file(gimms_files)),
                  ext = "tif")

# do it in manual batches (need at least 2 files, otherwise error in gimms function)
# i_batch <- 65:69
# 
# rr <- rasterizeGimms(gimms_files[i_batch], 
#                      filename = tif_files[i_batch],
#                      ext = ext_custom,
#                      keep = c(0,1),
#                      overwrite = F)

# plot(rr[[1:9]])

# rm(rr)


# extract time series of chronos ------------------------------------------

dt_gimms <- foreach(
  i_file = tif_files,
  .final = rbindlist
) %do% {
  
  rr <- stack(i_file)
  
  xx <- raster::extract(rr, sf_coord)
  
  data.table(site = sf_coord$site,
             xx) %>% 
    melt(id.vars = "site")
  
}

dt_gimms[, c("year", "half_year", "half_year_i") := tstrsplit(variable, "_|[.]", keep = 4:6)]
dt_gimms[, ":="(year = as.numeric(year), 
                half_year_i = as.numeric(half_year_i))]
dt_gimms[half_year == "0712", half_year_i := half_year_i + 12]

# library(ggplot2)
# dt_gimms[site == "AST"] %>%
#   ggplot(aes(half_year_i, value))+
#   geom_point()+
#   facet_wrap(~year)


dt_gimms_out <- dt_gimms[, .(site,
                             year,
                             xx_half_month = half_year_i,
                             ndvi = value)]

saveRDS(dt_gimms_out, file = "data/gimms-ndvi-01.rds")




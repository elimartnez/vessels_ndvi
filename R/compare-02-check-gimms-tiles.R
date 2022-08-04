# extract the time series from GIMMS

library(gimms)
library(fs)
library(magrittr)
library(sf)
library(foreach)
library(data.table)
library(mapview)


# check extent of chronos -------------------------------------------------
dat_coord <- readRDS("data/chrono-00-meta.rds")
sf_coord <- st_as_sf(dat_coord,
                     coords = c("lon", "lat"),
                     crs = "+proj=longlat +datum=WGS84")
mapview(sf_coord)
ext_custom <- extent(-12, 50, 27, 55) #xmin xmax ymin ymax





# gimms 2 raster ----------------------------------------------------------


gimms_files <- fs::dir_ls("/mnt/CEPH_PROJECTS/CLIRSNOW/gimms/nc/")
tif_files <- path("/mnt/CEPH_PROJECTS/CLIRSNOW/gimms/tif-treestudy/",
                  path_ext_remove(path_file(gimms_files)),
                  ext = "tif")


# mapview --------------------------------------------------------------------

rr <- raster(tif_files[1])
plot(rr)

cells_sites <- cellFromXY(rr, as(sf_coord, "Spatial"))
rr2 <- rr
rr2[] <- NA
rr2[cells_sites] <- rr[cells_sites]

bmap <- mapviewGetOption("basemaps")[c(4, 1:3, 5)]

mapview(list(sf_coord, rr2), map.types = bmap)

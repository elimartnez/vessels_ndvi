

library(fs)
library(magrittr)
library(sf)
library(foreach)
library(data.table)
library(mapview)
library(MODISTools)

# check extent of chronos -------------------------------------------------

dat_coord <- readRDS("data/chrono-00-meta.rds")
sf_coord <- st_as_sf(dat_coord,
                     coords = c("lon", "lat"),
                     crs = "+proj=longlat +datum=WGS84")
mapview(sf_coord)
ext_custom <- extent(-12, 50, 27, 55) #xmin xmax ymin ymax





# modis 2 raster ----------------------------------------------------------

# dir_ls("data-modis/ndvi/") %>% 
#   lapply(fread, colClasses = "character") %>% 
#   rbindlist -> dat
# dat[, nrows := as.numeric(nrows)]
# dat[, ncols := as.numeric(ncols)]

dir_ls("data-modis/ndvi/") %>% 
  lapply(fread) %>% 
  rbindlist -> dat


dat1 <- dat[, .SD[1], site]

bb <- apply(dat1, 1, function(x){
  mt_bbox(xllcorner = x['xllcorner'],
          yllcorner = x['yllcorner'],
          cellsize = x['cellsize'],
          nrows = x['nrows'],
          ncols = x['ncols'])
})

lapply(bb, st_sf) %>% 
  do.call(rbind, .) -> sf_modis

# sf_modis[["id"]] <- as.character(1:24)
sf_modis[["id"]] <- "modis_pixels"
# st_geometry(sf_modis)

sf_modis %>% st_combine() -> sf_modis2


# mapview --------------------------------------------------------------------


bmap <- mapviewGetOption("basemaps")[c(4, 1:3, 5)]

mapview(list(sf_modis, sf_coord), label = "site", map.types = bmap)

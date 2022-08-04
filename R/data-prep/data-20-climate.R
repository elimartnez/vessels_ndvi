# read in climate data (eli - chelsa)

library(data.table)
library(magrittr)

dat1 <- fread("data-eli/v3/vessels_clima2.csv")
dat1 %>% summary
setnames(dat1, c("variable", "tmea", "spei"), c("site", "tmean", "SPEI"))

dat2 <- fread("data-eli/v5/vessels_clima2_second2.csv")
dat2 %>% summary
dat2 %>% str
dat2[, V1 := NULL]
setnames(dat2, c("variable", "tmea", "spei"), c("site", "tmean", "SPEI"))

dat <- rbindlist(list(dat1, dat2), use.names = T)
dat[, .N, site]
dat[, .N, year]

dat_coord <- readRDS("data/chrono-00-meta.rds")
dat[, .N, site] %>% merge(dat_coord) # -> all good

saveRDS(dat, "data/climate-01.rds")

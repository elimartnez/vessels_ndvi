library(magrittr)
library(data.table)
library(readxl)
library(foreach)

# fn <- "data-eli/v7/MODIFIED_Kh_cronos_Kh1.xlsx"
# fn_out <- "data/chrono-02-kh1.rds"

fn <- "data-eli/v7/MODIFIED_Kh_cronos_Kh2.xlsx"
fn_out <- "data/chrono-02-kh2.rds"

# meta --------------------------------------------------------


# 
tbl_meta <- read_excel(fn, "Metadata")
tbl_meta %>%
  data.table %>%
  .[, .(site = ID, lon = long_num, lat = lat_num, species = Species)] -> dat_coord

saveRDS(dat_coord, file = "data/chrono-00-meta.rds")
# 
# 
# # plot map
# sf_coord <- sf::st_as_sf(dat_coord,
#                          coords = c("lon", "lat"),
#                          crs = "+proj=longlat +datum=WGS84")
# mapview::mapview(sf_coord)


# chrono kh1 ------------------------------------------------------------------


# eli
# dat_eli <- fread("data-eli/v1/eli_data.csv", skip = 1)
# fread("data-eli/v1/eli_data.csv", nrows = 1) %>% 
#   as.character() %>% tolower() -> names(dat_eli)
# 
# melt(dat_eli, id.vars = c("site", "variable"),
#      variable.name = "year", variable.factor = F) -> dat_eli
# dat_eli[, year := as.numeric(year)]


# eli & others

excel_sheets(fn)


l_data <- foreach(i_sheet = excel_sheets(fn)[-1]) %do% {
  read_excel(fn, i_sheet) %>% 
    .[, 1:4] %>% 
    data.table -> dat
  # setnames(dat, 1, "year")
  # dat_long <- melt(dat, id.vars = "year", variable.factor = F)
  # dat_long[, site := i_sheet]
  # dat_long
  dat[, site := i_sheet]
  dat
}
# l_data %>% str

# dat_others <- rbindlist(l_data)
dat_others <- rbindlist(l_data, use.names = T, fill = T)
dat_others[, .(site, year, TRW, MVA, Kh)] %>% 
  melt(id.vars = c("site", "year")) -> dat_chrono
# dat_chrono <- rbindlist(list(dat_others_long, dat_eli), use.names = T)
# dat_chrono <- dat_chrono[!is.na(value)]

# check names and variables
dat_chrono[, .N, variable]
dat_chrono[, .N, site]
dat_chrono[, .N, site] %>% merge(dat_coord, all = T)
# --> all good

saveRDS(dat_chrono, file = fn_out)




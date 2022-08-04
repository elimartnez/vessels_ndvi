


# backup download ---------------------------------------------------------


# gpp
mt_batch_subset(df = dat_coord_sub,
                product = "MOD17A2H",
                band = "Gpp_500m",
                start = "2000-01-01",
                end = "2012-12-31",
                internal = F,
                out_dir = "data-modis/gpp/")

# et
mt_batch_subset(df = dat_coord_sub,
                product = "MOD16A2",
                band = "ET_500m",
                start = "2000-01-01",
                end = "2012-12-31",
                internal = F,
                out_dir = "data-modis/et/")

# lai
mt_batch_subset(df = dat_coord_sub,
                product = "MOD15A2H",
                band = "Lai_500m",
                start = "2000-01-01",
                end = "2012-12-31",
                internal = F,
                out_dir = "data-modis/lai/")

# fpar
mt_batch_subset(df = dat_coord_sub,
                product = "MOD15A2H",
                band = "Fpar_500m",
                start = "2000-01-01",
                end = "2012-12-31",
                internal = F,
                out_dir = "data-modis/fpar/")

# ndvi
mt_batch_subset(df = dat_coord_sub,
                product = "MOD13Q1",
                band = "250m_16_days_NDVI",
                start = "2000-01-01",
                end = "2012-12-31",
                internal = F,
                out_dir = "data-modis/ndvi/")



# test gimms extract ------------------------------------------------------


library(gimms)
library(fs)
library(magrittr)
library(sf)

dat_coord <- readRDS("data/chrono-00-meta.rds")
sf_coord <- st_as_sf(dat_coord,
                     coords = c("lon", "lat"),
                     crs = "+proj=longlat +datum=WGS84")
mapview::mapview(sf_coord)
ext_custom <- extent(-12, 50, 27, 55)

gimms_files <- fs::dir_ls("/mnt/CEPH_PROJECTS/CLIRSNOW/gimms/nc/")
tif_files <- path("/mnt/CEPH_PROJECTS/CLIRSNOW/gimms/tif-treestudy/",
                  path_ext_remove(path_file(gimms_files)),
                  ext = "tif")

rr <- rasterizeGimms(gimms_files[1:2], 
                     filename = tif_files[1:2],
                     ext = ext_custom,
                     keep = c(0,1),
                     overwrite = T)


# plot(rr[[1:9]])

xx <- raster::extract(rr, sf_coord)
xx



# temp plots --------------------------------------------------------------

dat[site == "BOI"] %>% 
  ggplot(aes(date, gpp))+
  geom_point()+
  facet_wrap(~ year(date), scales = "free_x")+
  scale_x_date(date_breaks = "3 months", date_labels = "%b")+
  theme_bw()


dat1 <- dat[site == "AST"]
dat1[, table(b_0, b_43)]
dat[, table(b_0, b_43)] 

dat1 %>% 
  ggplot(aes(date, lai, colour = b_43, shape = as.factor(b_0)))+
  geom_point()+
  facet_wrap(~year(date), scales = "free_x")+
  scale_x_date(date_breaks = "3 months", date_labels = "%b")+
  scale_shape_manual(values = c(19, 1))+
  theme_bw()




# lai - check qc bits ---------------------------------------------------------------------


dat <- readRDS("data/lai-01.rds")

all_sites <- sort(unique(dat$site))

pdf(file = "fig/raw-ts/lai-qc.pdf",
    width = 14, height = 8)
walk(all_sites, function(i_site){
  
  gg <- dat[site == i_site] %>% 
    ggplot(aes(date, lai, colour = b_43, shape = as.factor(b_0)))+
    geom_point()+
    facet_wrap(~year(date), scales = "free_x")+
    scale_x_date(date_breaks = "3 months", date_labels = "%b")+
    scale_shape_manual(values = c(19, 1))+
    theme_bw()+
    ggtitle(i_site)
  
  print(gg)
  
  
})
dev.off()

# -> need to mask out :
#       b_0 == 1 (non-good quality)
#       b_43 == 01 (significant clouds)




# check negative corr ndvi vs trw for SAL, ROB, AL2 -----------------------



dat_chrono[variable == "TRW"]

eli_var <- "TRW"
eli_var <- "MVA"
dat_merge <- merge(dat_gimms[nn == (period + 1)], 
                   dat_chrono[variable == eli_var],
                   by = c("site", "year"))

dat_merge[site %in% c("SAL", "ROB", "AL2"), 
          .(site, year, value)] %>% 
  unique %>% 
  ggplot(aes(year, value))+
  geom_point()+
  facet_wrap(~site)


dat_merge[, 
          .(site, year, value)] %>% 
  unique %>% 
  ggplot(aes(year, value))+
  geom_point()+
  facet_wrap(~site, scales = "free_y")



dat_merge[xx_half_month == 16 & period == 0, 
          .(site, year, mean_ndvi, value)] %>% 
  unique %>% 
  ggplot(aes(year, mean_ndvi))+
  geom_point()+
  facet_wrap(~site, scales = "free_y")



dat_merge[xx_half_month == 13 & period == 2, 
          .(site, year, mean_ndvi, value)] %>% 
  unique %>% 
  ggplot(aes(year, mean_ndvi))+
  geom_point()+
  facet_wrap(~site, scales = "free_y")



dat_merge[xx_half_month == 13 & period == 2, 
          .(site, year, mean_ndvi, value)] %>% 
  unique %>% 
  ggplot(aes(value, mean_ndvi))+
  geom_point()+
  facet_wrap(~site, scales = "free")



# compare the two pca's ---------------------------------------------------


pca1 %>% str
pca_site %>% str

plot(pca1$sdev, pca_site$sdev)

pca1$x %>% str

plot(pca_site$rotation[, "PC1"], pca1$x[, "PC1"])
plot(pca_site$rotation[, "PC2"], pca1$x[, "PC2"])
plot(pca_site$rotation[, "PC3"], pca1$x[, "PC3"])


pca_corr_mva <- pca1

eli_var <- "TRW"
pca_corr_trw <- pca1


pca_corr_mva %>% str
pca_corr_trw %>% str

plot(pca_corr_mva$sdev, pca_corr_trw$sdev)

plot(pca_corr_mva$rotation[, "PC1"], pca_corr_trw$rotation[, "PC1"])
plot(pca_corr_mva$rotation[, "PC2"], pca_corr_trw$rotation[, "PC2"])
plot(pca_corr_mva$rotation[, "PC3"], pca_corr_trw$rotation[, "PC3"])

plot(pca_corr_mva$x[, "PC1"], pca_corr_trw$x[, "PC1"])
plot(pca_corr_mva$x[, "PC2"], pca_corr_trw$x[, "PC2"])
plot(pca_corr_mva$x[, "PC3"], pca_corr_trw$x[, "PC3"])



# detrend deseason --------------------------------------------------------

# detrend (and deseasonalize) NDVI

library(data.table)
library(mgcv)
library(magrittr)
library(ggplot2)

dat_ndvi <- readRDS("data/gimms-ndvi-01.rds")

dat_ndvi
dat1 <- dat_ndvi[site == "CFI"]

ggplot(dat1, aes(xx_half_month, ndvi))+
  geom_point()+
  facet_wrap(~year)

gam1 <- gam(ndvi ~ s(xx_half_month, bs = "cc") + year, data = dat1)
gam1 %>% summary
plot(gam1)
gam.check(gam1)

newdat1 <- dat1[, .(year = mean(year), xx_half_month)]
dat1[, fit1 := predict(gam1, newdat1)]

ggplot(dat1, aes(xx_half_month, ndvi))+
  geom_point()+
  geom_line(aes(y = fit1))+
  facet_wrap(~year)

ggplot(dat1, aes(xx_half_month, ndvi))+
  # geom_point()+
  geom_line(aes(y = fit1))+
  facet_wrap(~year)


dat_ndvi2[xx_half_month == 12] %>% 
  ggplot(aes(year, ndvi_trend))+
  geom_point()+
  facet_wrap(~site)

dat_ndvi2[site == "CFI"] %>% 
  ggplot(aes(xx_half_month, ndvi))+
  geom_point()+
  geom_line(aes(y = ndvi_detrend))+
  facet_wrap(~year)




# compare two Kh for PAO and ORI ------------------------------------------

dat <- data.table(read_excel(fn, "PAO"))
dat <- data.table(read_excel(fn, "ORI"))

dat %>% pairs



# compare chrono corr: pearson, spearman ----------------------------------

dat1 <- dat_merge[period < agg_period/15, 
                  .(pearson = cor(value, mean_ndvi_detrend, method = "p"),
                    spearman = cor(value, mean_ndvi_detrend, method = "s"),
                    kendall = cor(value, mean_ndvi_detrend, method = "k")),
                  .(site, nyear, xx_half_month, period)]
dat1 %>% 
  ggplot(aes(pearson, spearman))+
  geom_abline()+
  geom_point()+
  facet_wrap(~site)

dat1 %>% 
  ggplot(aes(pearson, kendall))+
  geom_abline()+
  geom_point()+
  facet_wrap(~site)

dat1 %>% 
  ggplot(aes(site, spearman - pearson))+
  geom_hline(yintercept = 0)+
  geom_boxplot()+
  theme_bw()

dat1 %>% 
  ggplot(aes(pearson, spearman - pearson))+
  geom_hline(yintercept = 0)+
  geom_point()+
  facet_wrap(~site)+
  theme_bw()

dat1 %>% 
  ggplot(aes(pearson, kendall - pearson))+
  geom_hline(yintercept = 0)+
  geom_point()+
  facet_wrap(~site)+
  theme_bw()


dat_merge[period < agg_period/15 & site == "DOL" & xx_half_month >=5 & xx_half_month <= 20] %>% 
  ggplot(aes(value, mean_ndvi_detrend))+
  geom_point()+
  facet_grid(period ~ xx_half_month)


# understand PCA decomposition --------------------------------------------


pca1_rot %>% dim
pca1_x %>% dim

mat_pca %>% dim

mat_pca2 <- pca1_rot %*% t(pca1_x)
plot(mat_pca, mat_pca2)

mat_pca_full <- pca1$rotation %*% t(pca1$x)
plot(mat_pca, mat_pca_10)


# test plot clim pca ------------------------------------------------------

# grouping: PC (4x1) - climvar x month
dat_pca_clim[tree_var == "Kh"] %>% 
  split(by = "variable") %>% 
  lapply(function(x){
    x %>% 
      split(by = "climvar_fct") %>% 
      lapply(function(y){
        y %>% 
          ggplot(aes(climvar_value, value))+
          geom_point(size = 1, colour = "grey50")+ 
          # geom_smooth(se = F, colour = "black", linetype = "dashed")+
          # geom_smooth(method = lm, formula = y ~ x, se = F, colour = "black", linetype = "dashed")+
          geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = F, colour = "black", linetype = "dashed")+
          scale_y_continuous(n.breaks = 4)+
          facet_grid(climvar_fct ~ month_f)+
          theme_bw()+
          theme(panel.grid = element_blank())+
          ylab("PC site coefficient")+
          xlab(NULL)
      }) %>% 
      cowplot::plot_grid(plotlist = ., ncol = 1)
  }) %>% 
  cowplot::plot_grid(plotlist = ., ncol = 1,
                     labels = paste0("(", letters[1:4], ")", " ", "PC", 1:4))



# test gam expl variance --------------------------------------------------

load("data/paper-fig-01.rda")

dat1 <- dat_pca_clim[tree_var == "Kh" & variable == "PC1" & climvar == "tmean" & month == 3]

gm1 <- gam(value ~ s(climvar_value, k = 3), data = dat1)
gm1 %>% summary
broom::tidy(gm1)
broom::glance(gm1)

dat_gam_summ <- dat_pca_clim[tree_var == "Kh", 
                             broom::tidy(gam(value ~ s(climvar_value, k = 3), select = T)),
                             .(variable, climvar, month)]

dat_gam_summ

# EOF ---------------------------------------------------------------------



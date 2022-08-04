# compare MODIS and GIMMS NDVI to check for homogenous forests


library(data.table)
library(magrittr)
library(ggplot2)
library(lubridate)
library(patchwork)



# prep data ---------------------------------------------------------------


dat_modis <- readRDS("data/modis-ndvi-01.rds")
dat_gimms <- readRDS("data/gimms-ndvi-01.rds")

dat_modis[, .N, site]
dat_gimms[, .N, site]

dat_modis[, yday(date)] %>% unique %>% sort %>% length
dat_gimms$xx_half_month %>% table

dat_modis[, year := year(date)]
dat_modis[, yday := yday(date)]
dat_modis[, yday_int := as.numeric(as.factor(yday))]


dat_merge <- merge(
  dat_gimms[, .(site, year, xx = xx_half_month, ndvi_gimms = ndvi)],
  dat_modis[, .(site, year, xx = yday_int, ndvi_modis = ndvi_f)]
)

# correlations MODIS vs GIMMS ---------------------------------------------

# dat_merge %>%
# dat_merge[xx %in% c(7:19)] %>% 
dat_merge[xx %in% c(10:16)] %>%
  ggplot(aes(ndvi_modis, ndvi_gimms))+
  geom_point()+
  facet_wrap(~site, scales = "free")

dat_merge[, cor(ndvi_modis, ndvi_gimms, use = "p"), site]
dat_merge[xx %in% c(7:19), cor(ndvi_modis, ndvi_gimms, use = "p"), site]
dat_merge[xx %in% c(10:16), cor(ndvi_modis, ndvi_gimms, use = "p"), site]

# problem: seasonality is causing all the correlation
# if only summer, then maybe. but still doubtfull


# comparison moving window corr plots -------------------------------------

min_month <- 1 #incl
max_month <- 12 #incl
# sites_to_remove <- c("ORI", "PAO")
agg_period <- 120


dat_chrono <- readRDS("data/chrono-01.rds")
# dat_chrono <- dat_chrono[! site %in% sites_to_remove]
dat_gimms_mv <- readRDS("data/gimms-ndvi-02-moving.rds")
dat_modis_mv <- readRDS("data/modis-ndvi-02-moving.rds")


# ** gimms full --------------------------------------------------------------



# 
# dat_merge_gimms <- merge(dat_gimms_mv[nn == (period + 1) & year >= 2000], 
#                    dat_chrono[variable == "MVA"],
#                    by = c("site", "year"))

dat_merge_gimms <- merge(dat_gimms_mv[nn == (period + 1)], 
                         dat_chrono[variable == "MVA"],
                         by = c("site", "year"))

dat_merge_gimms[, nyear := length(unique(year)), site]

# subset to 1 year
dat_cor_gimms <- dat_merge_gimms[period < agg_period/15, 
                                 cor.test(value, mean_ndvi)[c("estimate", "p.value")],
                                 .(site, nyear, xx_half_month, period)]

dat_cor_gimms[, xx_month := xx_half_month/2 + 0.5]
dat_cor_gimms[, yy_period_days := period * 15]
dat_cor_gimms[, yy_plot := yy_period_days + 15/2]

# # subset months
# dat_cor <- dat_cor[xx_month >= min_month & xx_month < (max_month + 1)]
# dat_cor[, .N, .(site)]

# # reorder sites
# dat_cor[, avg_estimate_jja := mean(estimate[xx_month >= 6 & xx_month < 9]), .(site)]
# dat_cor[, avg_estimate := mean(estimate), .(site)]
# dat_cor[, site_f := fct_reorder(site, avg_estimate)]


# ** gimms modis period ---------------------------------------------------



dat_merge_gimms2 <- merge(dat_gimms_mv[nn == (period + 1) & year >= 2000],
                          dat_chrono[variable == "MVA"],
                          by = c("site", "year"))

# dat_merge_gimms <- merge(dat_gimms_mv[nn == (period + 1)], 
#                          dat_chrono[variable == "MVA"],
#                          by = c("site", "year"))

dat_merge_gimms2[, nyear := length(unique(year)), site]

# subset to 1 year
dat_cor_gimms2 <- dat_merge_gimms2[period < agg_period/15, 
                                   cor.test(value, mean_ndvi)[c("estimate", "p.value")],
                                   .(site, nyear, xx_half_month, period)]

dat_cor_gimms2[, xx_month := xx_half_month/2 + 0.5]
dat_cor_gimms2[, yy_period_days := period * 15]
dat_cor_gimms2[, yy_plot := yy_period_days + 15/2]

# # subset months
# dat_cor <- dat_cor[xx_month >= min_month & xx_month < (max_month + 1)]
# dat_cor[, .N, .(site)]

# # reorder sites
# dat_cor[, avg_estimate_jja := mean(estimate[xx_month >= 6 & xx_month < 9]), .(site)]
# dat_cor[, avg_estimate := mean(estimate), .(site)]
# dat_cor[, site_f := fct_reorder(site, avg_estimate)]



# ** modis ----------------------------------------------------------------



dat_merge_modis <- merge(dat_modis_mv[nn == (period + 1)], 
                         dat_chrono[variable == "MVA"],
                         by = c("site", "year"))
dat_merge_modis[, nyear := length(unique(year)), site]



# subset to 1 year
dat_cor_modis <- dat_merge_modis[period < agg_period/16, 
                                 cor.test(value, mean_ndvi_f)[c("estimate", "p.value")],
                                 .(site, nyear, doy, period)]

dat_cor_modis[, xx_date := ymd("1999-01-01")]
yday(dat_cor_modis$xx_date) <- dat_cor_modis$doy

dat_cor_modis[, yy_period_days := period * 16]


# ** loop plot ------------------------------------------------------------



pdf("fig/07-comparison-modis-gimms-MVA/comparison-01.pdf",
    width = 15, height = 4)
for(i_site in sort(unique(dat_chrono$site))){
  
  dat_merge_gimms[site == i_site, range(year)] %>% 
    paste(collapse = "-") -> period_gimms
  dat_merge_gimms2[site == i_site, range(year)] %>% 
    paste(collapse = "-") -> period_gimms2
  dat_merge_modis[site == i_site, range(year)] %>% 
    paste(collapse = "-") -> period_modis
  
  
  gg1 <- 
  dat_cor_gimms[site == i_site] %>% 
    ggplot(aes(xx_month, yy_plot, fill = estimate))+
    geom_raster()+
    scale_fill_gradient2("cor", limits = c(-1,1))+
    coord_cartesian(expand = F)+
    scale_x_continuous(breaks = c(min_month : max_month) - 0.25,
                       labels = month.abb[min_month : max_month])+
    scale_y_continuous(breaks = c(0, 30, 60, 90, 120))+
    facet_wrap(~site)+
    theme_bw()+
    xlab("End of aggregation period")+
    ylab("Number of days used to aggregate NDVI")+
    ggtitle("GIMMS (full period)", period_gimms)
  
  gg2 <- 
    dat_cor_gimms2[site == i_site] %>% 
    ggplot(aes(xx_month, yy_plot, fill = estimate))+
    geom_raster()+
    scale_fill_gradient2("cor", limits = c(-1,1))+
    coord_cartesian(expand = F)+
    scale_x_continuous(breaks = c(min_month : max_month) - 0.25,
                       labels = month.abb[min_month : max_month])+
    scale_y_continuous(breaks = c(0, 30, 60, 90, 120))+
    facet_wrap(~site)+
    theme_bw()+
    xlab("End of aggregation period")+
    ylab("Number of days used to aggregate NDVI")+
    ggtitle("GIMMS (same period as MODIS)", period_gimms2)
  
  gg3 <- 
  dat_cor_modis[site == i_site] %>% 
    ggplot(aes(xx_date, yy_period_days, fill = estimate))+
    geom_raster()+
    scale_fill_gradient2("cor", limits = c(-1,1))+
    coord_cartesian(expand = F)+
    scale_x_date(date_breaks = "month",
                 date_labels = "%b")+
    scale_y_continuous()+
    facet_wrap(~site)+
    theme_bw()+
    xlab("End of aggregation period")+
    ylab("Number of days used to aggregate NDVI")+
    ggtitle("MODIS", period_modis)
  
  
  gg_out <- gg1 + gg2 + gg3 + plot_layout(guides = "collect")
  print(gg_out)
  
}
dev.off()





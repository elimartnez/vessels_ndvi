# PCA stuff


library(data.table)
library(magrittr)
library(ggplot2)
library(purrr)
library(lubridate)
library(patchwork)
library(forcats)
library(effects)

dat_chrono <- readRDS("data/chrono-01.rds")
dat_chrono[, table(variable, site)]


dat_gimms <- readRDS("data/gimms-ndvi-02-moving.rds")
dat_gimms$site %>% table

dat_gimms[nn != (period + 1)]
dat_chrono[, .N, variable]


dat_merge <- merge(dat_gimms[nn == (period + 1)], 
                   dat_chrono[variable == "MVA"],
                   by = c("site", "year"))
dat_merge[, nyear := length(unique(year)), site]

# subset to 1 year
dat_cor <- dat_merge[period < 360/15, 
                     cor.test(value, mean_ndvi)[c("estimate", "p.value")],
                     .(site, nyear, xx_half_month, period)]

dat_cor[, xx_month := xx_half_month/2 + 0.5]
dat_cor[, yy_period_days := period * 15]

dat_cor[, .N, .(site)]

dat_cor[site == "AL1"] %>% 
ggplot(aes(xx_month, yy_period_days, fill = estimate))+
  geom_raster()+
  geom_point(data = dat_cor[site == "AL1" & p.value < 0.1], size = 0.1, colour = "grey20")+
  scale_fill_gradient2(limits = c(-1,1))+
  coord_cartesian(expand = F)+
  scale_x_continuous(breaks = 1:12 - 0.25,
                     labels = month.abb)+
  scale_y_continuous()+
  facet_wrap(~site + nyear, labeller = label_both)+
  theme_bw()+
  theme(legend.position = "none")+
  # ggtitle(paste0("Correlation between different averages of NDVI and annual ", eli_var),
  #         "dots indicate p-value < 0.1; correlation is pearson")+
  # xlab("End of aggregation period")+
  # ylab("Number of days used to aggregate NDVI")
  xlab(NULL)+ylab(NULL)
# gg_cor


dat_cor[site != "ORI" & site != "PAO"] %>% 
  dcast(xx_month + yy_period_days ~ site, value.var = "estimate") -> dat_wide

dat_wide[, -c("xx_month", "yy_period_days"), with = F] %>% 
  as.matrix -> mat_pca


pca1 <- prcomp(mat_pca, scale. = T, rank. = 20)
plot(pca1$sdev^2)
summary(pca1)

dat_plot_pca <- data.table(pca1$rotation[, 1:6])
dat_plot_pca[, site := rownames(pca1$rotation)]
dat_plot_pca_melt <- melt(dat_plot_pca, id.vars = "site")
setkey(dat_plot_pca_melt, variable, value)
dat_plot_pca_melt[, site_f := fct_inorder(site)]


dat_plot_pca_melt %>% 
  ggplot(aes(site_f, value))+
  geom_point()+
  facet_wrap(~variable)


# correlation plot

cbind(dat_wide[, .(xx_month, yy_period_days)], pca1$x[, 1:6]) %>% 
  melt(id.vars = c("xx_month", "yy_period_days")) %>% 
  ggplot(aes(xx_month, yy_period_days, fill = value))+
  geom_raster()+
  scale_fill_gradient2()+
  coord_cartesian(expand = F)+
  scale_x_continuous(breaks = 1:12 - 0.25,
                     labels = month.abb)+
  scale_y_continuous()+
  facet_wrap(~variable)+
  theme_bw()

# maps
dat_meta <- readRDS("data/chrono-00-meta.rds")

dat_plot_pca_melt %>% 
  merge(dat_meta, by = "site") -> dat_map

ggplot(dat_map, aes(lon, lat, colour = value))+
  borders(xlim = range(dat_meta$lon), ylim = range(dat_meta$lat))+
  geom_point(size = 3)+
  facet_wrap(~variable)+
  scale_colour_gradient2()+
  coord_equal(xlim = range(dat_meta$lon), ylim = range(dat_meta$lat))+
  theme_bw()


# correlate PCs to climate
dat_clim <- readRDS("data/climate-01.rds")

dat_gimms[, .(site, year)] %>% 
  unique %>% 
  merge(dat_clim) %>% 
  .[, lapply(.SD, mean), .(site, month)] -> dat_plot_clim

merge(dat_plot_clim,
      dat_plot_pca,
      by = "site") %>% 
  melt(id.vars = c("site", "month", "PPET"), measure.vars = paste0("PC", 1:6)) %>% 
  merge(dat_meta, by = "site") -> dat_pca_clim


dat_pca_clim[variable %in% paste0("PC", 1:4) & month %in% c(1:12)] %>% 
  ggplot(aes(PPET, value))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # geom_point()+ 
  geom_text(aes(label = site, colour = species), size = 2)+
  geom_smooth(method = lm)+
  facet_grid(variable ~ month, scales = "free_x")+
  theme_bw()

dat_pca_clim[variable %in% paste0("PC", 1:4) & month %in% c(1)] %>% 
  ggplot(aes(species, value))+
  # geom_hline(yintercept = 0)+
  # geom_vline(xintercept = 0)+
  # geom_point()+ 
  geom_text(aes(label = site, colour = species), size = 4)+
  geom_smooth(method = lm)+
  theme_bw()+
  facet_wrap(~variable)


biplot(pca1)
biplot(pca1, choices = c(1,3))



lm1 <- lm(value ~ PPET*month*species,
           dat_pca_clim[variable == "PC1"])
lm1 %>% anova
step(lm1)
lm2 <- lm(value ~ PPET*species,
          dat_pca_clim[variable == "PC1"])
allEffects(lm2) %>% plot

# correlate values to PC
cbind(dat_wide, PCx = pca1$x[, 1]) %>% 
  melt(id.vars = c("xx_month", "yy_period_days", "PCx")) %>% 
  ggplot(aes(PCx, value, colour = xx_month))+
  geom_point()+
  facet_wrap(~variable)+
  theme_bw()




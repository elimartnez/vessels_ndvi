# plot trends and average seasonal

library(data.table)
library(magrittr)
library(ggplot2)


dat_ndvi <- readRDS("data/gimms-ndvi-03-detrend.rds")

sites_to_remove <- c("CHL", "EIC", "BOI")
dat_ndvi <- dat_ndvi[!site %in% sites_to_remove]

# trend -------------------------------------------------------------------


dat_ndvi[, .(site, year, ndvi_trend)] %>% unique -> dat_ndvi_trend

gg_trend <- dat_ndvi_trend %>% 
  ggplot(aes(year, ndvi_trend))+
  geom_line()+
  facet_wrap(~site)+
  theme_bw()+
  xlab(NULL)+
  ylab("NDVI trend")

ggsave(gg_trend,
       filename = "fig/sm/trends.png",
       width = 12, height = 8)



# season ------------------------------------------------------------------



dat_ndvi[, .(site, xx_half_month, ndvi_avg_season)] %>% unique -> dat_ndvi_season

gg_season <-
  dat_ndvi_season %>% 
  ggplot(aes(xx_half_month/2 + 0.5, ndvi_avg_season))+
  geom_line()+
  facet_wrap(~site)+
  theme_bw()+
  xlab("month")+
  scale_x_continuous(breaks = seq(1,12,by=2))+
  ylab("NDVI average season")

ggsave(gg_season,
       filename = "fig/sm/season.png",
       width = 12, height = 8)




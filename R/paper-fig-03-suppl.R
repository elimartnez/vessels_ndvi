# paper fig


library(data.table)
library(magrittr)
library(ggplot2)
library(lubridate)
library(forcats)
library(foreach)
library(fs)

library(scico)
library(ggrepel)
library(patchwork)
library(cluster)
library(lemon)


fig_path_suppl <- "C:/Users/MMatiu/Dropbox/NDVI-vessels/Manuscript/fig-table-suppl/"

load("data/paper-fig-01.rda")




# pca - varexp ------------------------------------------------------------

dat_pca_varexp[, PC_num := as.integer(substr(PC, 3, 99))]

gg_varexp <- dat_pca_varexp %>% 
  ggplot(aes(PC_num, prop_variance))+
  geom_point()+
  geom_line(aes(group = variable))+
  geom_vline(xintercept = 4.5, linetype = "dashed")+
  facet_wrap(~tree_var)+
  theme_bw()+
  xlab("PC number")+
  ylab("Proportion of variance explained")

ggsave(gg_varexp,
       filename = path(fig_path_suppl, "pca-variance-explained.png"),
       width = 9, height = 3)


# Kh additional --------------------------------------------


## individual correlation plots --------------------------------------------


gg <-
  dat_cor[tree_var == "Kh"] %>% 
  
  ggplot(aes(xx_month, yy_plot, fill = estimate))+
  geom_raster()+
  scale_fill_gradient2("Pearson\ncorrelation", limits = c(-1,1))+
  coord_cartesian(expand = F)+
  scale_x_continuous(breaks = c(min_month : max_month) - 0.25,
                     labels = month.abb[min_month : max_month])+
  scale_y_continuous(breaks = c(0, 30, 60, 90, 120))+
  facet_wrap(~site_f)+
  theme_bw()+
  xlab("End of aggregation period")+
  ylab("Number of days used to aggregate NDVI")


ggsave(gg,
       file = path(fig_path_suppl, "Kh-individual-correlations.png"), 
       width = 12, height = 8, units = "in")



## climate corr full -------------------------------------------------------


gg_tmean <-
  dat_pca_clim[tree_var == "Kh" & climvar == "tmean"] %>% 
  ggplot(aes(climvar_value, value))+
  geom_smooth(method = lm, formula = y ~ x, se = F, colour = "black", linetype = "dashed")+
  # geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = F, colour = "black", linetype = "dashed")+
  geom_point(aes(colour = species))+ 
  geom_text_repel(aes(label = site, colour = species), size = 2)+
  scale_color_brewer("Species", palette = "Set1")+
  scale_y_continuous(n.breaks = 4)+
  facet_grid(variable ~ month_f, scales = "free")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  ylab("PC site coefficient")+
  xlab("Temperature [deg C]")



ggsave(gg_tmean,
       file = path(fig_path_suppl, "Kh-pca-climate-tmean.png"),
       width = 20, height = 10)



gg_prec <-
  dat_pca_clim[tree_var == "Kh" & climvar == "prec"] %>% 
  ggplot(aes(climvar_value, value))+
  geom_smooth(method = lm, formula = y ~ x, se = F, colour = "black", linetype = "dashed")+
  # geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = F, colour = "black", linetype = "dashed")+
  geom_point(aes(colour = species))+ 
  geom_text_repel(aes(label = site, colour = species), size = 2)+
  scale_color_brewer("Species", palette = "Set1")+
  scale_y_continuous(n.breaks = 4)+
  facet_grid(variable ~ month_f, scales = "free")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  ylab("PC site coefficient")+
  xlab("Precipitation [mm]")



ggsave(gg_prec,
       file = path(fig_path_suppl, "Kh-pca-climate-prec.png"),
       width = 20, height = 10)




gg_ppet <-
  dat_pca_clim[tree_var == "Kh" & climvar == "PPET_ratio"] %>% 
  ggplot(aes(climvar_value, value))+
  geom_smooth(method = lm, formula = y ~ x, se = F, colour = "black", linetype = "dashed")+
  # geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = F, colour = "black", linetype = "dashed")+
  geom_point(aes(colour = species))+ 
  geom_text_repel(aes(label = site, colour = species), size = 2)+
  scale_color_brewer("Species", palette = "Set1")+
  scale_y_continuous(n.breaks = 4)+
  facet_grid(variable ~ month_f, scales = "free")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  ylab("PC site coefficient")+
  xlab("P/PET [ratio]")



ggsave(gg_ppet,
       file = path(fig_path_suppl, "Kh-pca-climate-ppet.png"),
       width = 20, height = 10)





# MVA analysis ------------------------------------------------------------



## pca  --------------------------------------------------------------------



## climate -----------------------------------------------------------------






# TRW analysis ------------------------------------------------------------












# optional ----------------------------------------------------------------



## time series chrono ------------------------------------------------------



# depends




## time series ndvi --------------------------------------------------------



# ~24 pages (only zenodo or so)




## comparison gimms - modis ------------------------------------------------

# ?? did we really check it?







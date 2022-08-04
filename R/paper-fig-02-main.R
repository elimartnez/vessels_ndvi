# paper fig


library(data.table)
library(magrittr)
library(ggplot2)
library(lubridate)
library(forcats)
library(foreach)
library(fs)
library(mgcv)

library(scico)
library(ggrepel)
library(patchwork)
library(cluster)
library(lemon)

library(flextable)
library(officer)

# fig_path_main <- "C:/Users/MMatiu/Dropbox/NDVI-vessels/Manuscript/fig-table-main/"
fig_path_main <- "paper-fig-table/"

load("data/paper-fig-01.rda")



# fig 2 - explanatory corr plot -------------------------------------------

gg_fig2 <-
dat_cor[tree_var == "Kh" & site %in% c("SAL", "HIN", "MUQ", "TU2")] %>% 
  
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

# gg_fig2

ggsave(gg_fig2,
       file = path(fig_path_main, "Figure 2.png"), 
       width = 8, height = 4, units = "in")




# fig 3 - pca with coef ---------------------------------------------------

lim_col <- range(dat_pca_cor$value)
lim_coef <- range(dat_pca_site$value)


gg_fig3_sub1 <-
  dat_pca_cor[tree_var == "Kh"] %>% 
  ggplot(aes(xx_month, yy_plot, fill = value))+
  geom_raster()+
  scale_fill_scico("PC loading",
                   palette = "roma", rescaler = scales::rescale_mid, limits = lim_col)+
  coord_cartesian(expand = F)+
  scale_x_continuous(breaks = c(min_month : max_month) - 0.25,
                     labels = month.abb[min_month : max_month])+
  scale_y_continuous()+
  facet_grid(.~variable)+
  theme_bw()+
  xlab("End of aggregation period")+
  ylab("Number of days \nused to aggregate NDVI")



gg_fig3_sub2 <-
  dat_pca_site[tree_var == "Kh"] %>% 
  ggplot(aes(value, fct_rev(site)))+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_point()+
  facet_grid(.~variable)+
  scale_x_symmetric()+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  ylab("Site")+
  xlab("PC site coefficients")





gg_fig3 <- (gg_fig3_sub1 / gg_fig3_sub2)+
  plot_layout(heights = c(0.25, 0.5))



ggsave(gg_fig3,
       file = path(fig_path_main, "Figure 3.png"), 
       width = 12, height = 6, units = "in")



# fig 4 - climate scatters (gam) ------------------------------------------------

dat_pca_clim[, climvar_fct := fct_recode(
  factor(climvar, levels = c("tmean", "prec", "PPET_ratio")),
  "Temperature [deg C]" = "tmean",
  "Precipitation [mm]" = "prec",
  "P/PET [ratio]" = "PPET_ratio" 
)]


f_gam <- function(dat){
  xx <- seq(min(dat$climvar_value),
            max(dat$climvar_value),
            length.out = 50)
  gm <- gam(value ~ s(climvar_value, k = 3), select = T, data = dat)
  data.table(xx = xx,
             yy = as.vector(predict(gm, data.table(climvar_value = xx))))
}

s_gam <- function(dat){
  gm <- gam(value ~ s(climvar_value, k = 3), select = T, data = dat)
  data.table(deviance.explained = (gm$null.deviance - gm$deviance) / gm$null.deviance,
             broom::tidy(gm))
}

dat_pca_clim_sgam <- dat_pca_clim[tree_var == "Kh", 
                                  s_gam(.SD),
                                  .(variable, climvar_fct, month_f)]
dat_pca_clim_fgam <- dat_pca_clim[tree_var == "Kh", 
                                  f_gam(.SD),
                                  .(variable, climvar_fct, month_f)]
dat_pca_clim_gamplot <- merge(dat_pca_clim_sgam,
                              dat_pca_clim_fgam)
dat_pca_clim_gamplot[, sig_yn := ifelse(p.value < 0.1, "p<0.1", "n.s.")]

gg_fig4 <- lapply(sort(unique(dat_pca_clim$variable)), function(x_pc){
  
  ggplot(dat_pca_clim[tree_var == "Kh" & variable == x_pc],
         aes(climvar_value, value))+
    geom_point(size = 1, colour = "grey40")+ 
    geom_line(data = dat_pca_clim_gamplot[variable == x_pc], 
              inherit.aes = F,
              aes(xx, yy, colour = sig_yn, linetype = sig_yn))+
    scale_y_continuous(n.breaks = 4)+
    scale_color_manual("", values = c("p<0.1" = "black", "n.s." = "grey50"))+
    scale_linetype_manual("", values = c("p<0.1" = "solid", "n.s." = "dashed"))+
    facet_grid(month_f ~ climvar_fct, scales = "free", switch = "x")+
    theme_bw()+
    theme(panel.grid = element_blank(),
          strip.background.x = element_blank(),
          strip.placement = "outside",
          legend.position = "none")+
    ylab("PC site coefficient")+
    xlab(NULL)+
    ggtitle("")
  
}) %>% 
  cowplot::plot_grid(plotlist = ., 
                     labels = paste0("(", letters[1:4], ")", " ", "PC", 1:4))


ggsave(gg_fig4, 
       filename = path(fig_path_main, "Figure 4_gam.png"),
       # filename = "test.png",
       height = 15, width = 10)


# fig 4 - climate scatters (lm) ------------------------------------------------

dat_pca_clim[, climvar_fct := fct_recode(
  factor(climvar, levels = c("tmean", "prec", "PPET_ratio")),
  "Temperature [deg C]" = "tmean",
  "Precipitation [mm]" = "prec",
  "P/PET [ratio]" = "PPET_ratio" 
)]


f_lm <- function(dat){
  xx <- seq(min(dat$climvar_value),
            max(dat$climvar_value),
            length.out = 50)
  fit <- lm(value ~ climvar_value, data = dat)
  data.table(xx = xx,
             yy = as.vector(predict(fit, data.table(climvar_value = xx))))
}

s_lm <- function(dat){
  fit <- lm(value ~ climvar_value, data = dat)
  data.table(broom::tidy(fit))
}

dat_pca_clim_slm <- dat_pca_clim[tree_var == "Kh", 
                                  s_lm(.SD),
                                  .(variable, climvar_fct, month_f)]
dat_pca_clim_flm <- dat_pca_clim[tree_var == "Kh", 
                                  f_lm(.SD),
                                  .(variable, climvar_fct, month_f)]
dat_pca_clim_lmplot <- merge(dat_pca_clim_slm[term == "climvar_value"],
                              dat_pca_clim_flm)
dat_pca_clim_lmplot[, sig_yn := ifelse(p.value < 0.1, "p<0.1", "n.s.")]

gg_fig4 <- lapply(sort(unique(dat_pca_clim$variable)), function(x_pc){
  
  ggplot(dat_pca_clim[tree_var == "Kh" & variable == x_pc],
         aes(climvar_value, value))+
    geom_point(size = 1, colour = "grey40")+ 
    geom_line(data = dat_pca_clim_lmplot[variable == x_pc], 
              inherit.aes = F,
              aes(xx, yy, colour = sig_yn, linetype = sig_yn))+
    scale_y_continuous(n.breaks = 4)+
    scale_color_manual("", values = c("p<0.1" = "black", "n.s." = "grey50"))+
    scale_linetype_manual("", values = c("p<0.1" = "solid", "n.s." = "dashed"))+
    facet_grid(month_f ~ climvar_fct, scales = "free", switch = "x")+
    theme_bw()+
    theme(panel.grid = element_blank(),
          strip.background.x = element_blank(),
          strip.placement = "outside",
          legend.position = "none")+
    ylab("PC site coefficient")+
    xlab(NULL)+
    ggtitle("")
  
}) %>% 
  cowplot::plot_grid(plotlist = ., 
                     labels = paste0("(", letters[1:4], ")", " ", "PC", 1:4))


ggsave(gg_fig4, 
       filename = path(fig_path_main, "Figure 4_lm.png"),
       # filename = "test.png",
       height = 15, width = 10)



# table 2 - climate correlation summaries ---------------------------------


dat_pca_clim_sgam

dat_pca_clim_sgam[, p.value.sym := symnum(p.value, corr = FALSE,
                                            cutpoints = c(0,  .001,.01,.05, .1, 1),
                                            symbols = c("***","**","*","'"," "))]
dat_pca_clim_sgam[, deviance.explained.perc := scales::percent(deviance.explained)]
dat_pca_clim_sgam[, lbl := sprintf("%s%-4s", deviance.explained.perc, p.value.sym)]


ft_table2 <-
dat_pca_clim_sgam %>% 
  dcast(variable + month_f ~ climvar_fct, value.var = "lbl") %>% 
  flextable() %>% 
  set_header_labels(variable = "PC number", 
                    month_f = "Month") %>% 
  merge_v(j = "variable") %>% 
  valign(valign = "top") %>% 
  fix_border_issues() %>% 
  fontsize(size = 9) %>% 
  autofit()

read_docx() %>% 
  body_add_flextable(ft_table2) %>% 
  print(target = path(fig_path_main, "Table 2_gam.docx"))





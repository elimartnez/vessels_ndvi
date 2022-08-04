# some PCA stuff



library(data.table)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(scico)
library(lubridate)
library(forcats)




# settings ----------------------------------------------------------------

min_month <- 1 #incl
max_month <- 12 #incl
# sites_to_remove <- c("ORI", "PAO")
agg_period <- 120
# corr_fig_height <- 9
path_out <- "fig/05-refined2/"



# prep data ---------------------------------------------------------------

dat_chrono <- readRDS("data/chrono-01.rds")
# dat_chrono <- dat_chrono[! site %in% sites_to_remove]
dat_chrono[, table(variable, site)]


dat_gimms <- readRDS("data/gimms-ndvi-02-moving.rds")
# dat_gimms <- dat_gimms[! site %in% sites_to_remove]
dat_gimms$site %>% table


dat_clim <- readRDS("data/climate-01.rds")
dat_gimms[, .(site, year)] %>% 
  unique %>% 
  merge(unique(dat_chrono[, .(site, year)])) %>% 
  merge(dat_clim) %>% 
  .[, lapply(.SD, mean), .(site, month)] -> dat_plot_clim

saveRDS(dat_plot_clim, file = paste0(path_out, "data-avg-climate.rds"))

dat_plot_clim <- dat_plot_clim[month >= min_month & month <= max_month]

dat_meta <- readRDS("data/chrono-00-meta.rds")


# loop over tree ring measures -----------------------------------------------------



for(eli_var in c("TRW", "Kh", "MVA")){
  
  
  dat_merge <- merge(dat_gimms[nn == (period + 1)], 
                     dat_chrono[variable == eli_var],
                     by = c("site", "year"))
  dat_merge[, nyear := length(unique(year)), site]
  
  # subset to 1 year
  dat_cor <- dat_merge[period < agg_period/15, 
                       cor.test(value, mean_ndvi)[c("estimate", "p.value")],
                       .(site, nyear, xx_half_month, period)]
  
  dat_cor[, xx_month := xx_half_month/2 + 0.5]
  dat_cor[, yy_period_days := period * 15]
  dat_cor[, yy_plot := yy_period_days + 15/2]
  
  # subset months
  dat_cor <- dat_cor[xx_month >= min_month & xx_month < (max_month + 1)]
  
  
  dat_cor[, .N, .(site)]
  
  # reorder sites
  dat_cor[, avg_estimate_jja := mean(estimate[xx_month >= 6 & xx_month < 9]), .(site)]
  dat_cor[, avg_estimate := mean(estimate), .(site)]
  dat_cor[, site_f := fct_reorder(site, avg_estimate)]
  
  
  gg_cor <- dat_cor %>% 
    ggplot(aes(xx_month, yy_plot, fill = estimate))+
    geom_raster()+
    # geom_point(data = dat_cor[p.value < 0.1], size = 0.1, colour = "grey20")+
    scale_fill_gradient2("cor", limits = c(-1,1))+
    coord_cartesian(expand = F)+
    # scale_x_continuous(breaks = 1:12 - 0.25,
    #                    labels = month.abb)+
    scale_x_continuous(breaks = c(min_month : max_month) - 0.25,
                       labels = month.abb[min_month : max_month])+
    scale_y_continuous(breaks = c(0, 30, 60, 90, 120))+
    # facet_wrap(~site + nyear, labeller = label_both)+
    facet_wrap(~site_f)+
    theme_bw()+
    # ggtitle(paste0("Correlation between different averages of NDVI and annual ", eli_var),
    #         "dots indicate p-value < 0.1; correlation is pearson")+
    ggtitle(paste0("Correlation between different averages of NDVI and annual ", eli_var))+
    xlab("End of aggregation period")+
    ylab("Number of days used to aggregate NDVI")
  # gg_cor
  
  
  ggsave(gg_cor,
         file = paste0(path_out, eli_var, "_corr-plot_order-avg-cor.pdf"), 
         width = 18, height = 9, units = "in")
  
  
  dat_cor[, site_f := fct_reorder(site, avg_estimate_jja)]
  
  ggsave(gg_cor,
         file = paste0(path_out, eli_var, "_corr-plot_order-avg-cor-summer.pdf"), 
         width = 18, height = 9, units = "in")
  
  
  # pca ---------------------------------------------------------------------
  
  
  
  dat_cor %>% 
    dcast(xx_month + yy_plot ~ site, value.var = "estimate") -> dat_wide
  
  dat_wide[, -c("xx_month", "yy_plot"), with = F] %>% 
    as.matrix -> mat_pca
  
  pca1 <- prcomp(mat_pca, scale. = T)
  
  
  
  # variance explained
  dat_var_exp <- data.table(PC = paste0("PC", 1:length(pca1$sdev)),
                            prop_variance = pca1$sdev^2 / sum(pca1$sdev^2))
  dat_var_exp[, prop_variance_cumulative := cumsum(prop_variance)]
  dat_var_exp[, prop_variance := round(prop_variance, 3)]
  dat_var_exp[, prop_variance_cumulative := round(prop_variance_cumulative, 3)]
  fwrite(dat_var_exp, file = paste0(path_out, eli_var, "_pca-variance-explained.csv"))
  
  
  # biplot
  # pdf(paste0(path_out, eli_var, "_pca-biplot.pdf"),
  #     width = 10, height = 10)
  # biplot(pca1)
  # dev.off()
  
  
  # pca - corr plot
  gg_pca_cor <-
  cbind(dat_wide[, .(xx_month, yy_plot)], pca1$x[, 1:3]) %>% 
    melt(id.vars = c("xx_month", "yy_plot")) %>% 
    .[, value_sc := scales::rescale(value, to = c(-1, 1))] %>% 
    ggplot(aes(xx_month, yy_plot, fill = value_sc))+
    geom_raster()+
    scale_fill_gradient2()+
    coord_cartesian(expand = F)+
    # scale_x_continuous(breaks = 1:12 - 0.25,
    #                    labels = month.abb)+
    scale_x_continuous(breaks = c(min_month : max_month) - 0.25,
                       labels = month.abb[min_month : max_month])+
    scale_y_continuous()+
    facet_wrap(~variable)+
    theme_bw()+
    ggtitle(paste0("PCs of the correlations: ", eli_var))+
    xlab("End of aggregation period")+
    ylab("Number of days used to aggregate NDVI")
  
  
  ggsave(gg_pca_cor,
         file = paste0(path_out, eli_var, "_pca-corr-plot.pdf"), 
         width = 12, height = 4, units = "in")
  
  
  
  
  
  # correlate PCs to climate
  dat_plot_pca <- data.table(pca1$rotation[, 1:3])
  dat_plot_pca[, site := rownames(pca1$rotation)]
  
  for(i_climvar in c("PPET", "tmean", "prec")){
    
    merge(dat_plot_clim,
          dat_plot_pca,
          by = "site") %>% 
      melt(id.vars = c("site", "month", i_climvar), measure.vars = paste0("PC", 1:3)) %>% 
      merge(dat_meta, by = "site") -> dat_pca_clim
    
    dat_pca_clim[, month_f := factor(month.abb[month], levels = month.abb)]
    
    
    gg_pca_clim <- dat_pca_clim[variable %in% paste0("PC", 1:3) & month %in% c(1:12)] %>% 
      ggplot(aes_string(i_climvar, "value"))+
      # geom_hline(yintercept = 0)+
      # geom_vline(xintercept = 0)+
      geom_smooth(method = lm, se = F, colour = "black", linetype = "dashed")+
      geom_point(aes(colour = species))+ 
      # geom_text(aes(label = site, colour = species), size = 3)+
      geom_text_repel(aes(label = site, colour = species), size = 2)+
      scale_color_brewer(palette = "Set1")+
      facet_grid(variable ~ month_f, scales = "free_x", space = "free_x")+
      theme_bw()+
      ylab("PC loading")
    # gg_pca_clim
    
    
    ggsave(gg_pca_clim,
           file = paste0(path_out, eli_var, "_pca-scatter-clim_", i_climvar, ".pdf"),
           width = 20, height = 6)
    
    
  }
  
  
  # scatter cor with pc's
  
  # for(i_pc in 1:4){
  #   
  #   cbind(dat_wide, PCx = pca1$x[, i_pc]) %>% 
  #     melt(id.vars = c("xx_month", "yy_plot", "PCx")) -> dat_pca_scat
  #   
  #   # dat_pca_scat[, table(xx_month, month)]
  #   dat_pca_scat[, month := floor(xx_month)]
  #   dat_pca_scat[, month_f := factor(month.abb[month], levels = month.abb)]
  #   
  #   gg_pca_raw_rot <- dat_pca_scat %>% 
  #     ggplot(aes(value, PCx, colour = month_f))+
  #     geom_point()+
  #     scale_color_scico_d("")+
  #     facet_wrap(~variable)+
  #     theme_bw()+
  #     ylab("Rotated data (from PC)")+
  #     xlab("Original data (correlation)")+
  #     ggtitle(paste0("Scatter plot of original data versus rotated data by PCA: PC", i_pc))
  #   # gg_pca_raw_rot  
  #   
  #   ggsave(gg_pca_raw_rot,
  #          file = paste0(path_out, eli_var, "_pca-scatter-orig-rot_PC", i_pc, ".pdf"),
  #          width = 14, height = 8)
  #   
  #   
  #   
  # }
  
  
  
}





# EOF ---------------------------------------------------------------------




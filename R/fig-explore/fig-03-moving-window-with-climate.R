# plot moving window stuff

library(data.table)
library(magrittr)
library(ggplot2)
library(purrr)
library(lubridate)
library(patchwork)

dat_chrono <- readRDS("data/chrono-01.rds")
dat_chrono[, table(variable, site)]

dat_climate <- readRDS("data/climate-01.rds")


# gimms ndvi --------------------------------------------------------------


dat_gimms <- readRDS("data/gimms-ndvi-02-moving.rds")
dat_gimms$site %>% table

dat_gimms[nn != (period + 1)]
dat_chrono[, .N, variable]

dat_gimms[, .(site, year)] %>% 
  unique %>% 
  merge(dat_climate) %>% 
  .[, lapply(.SD, mean), .(site, month)] -> dat_plot_climate

# plot all

for(eli_var in c("TRW", "Kh", "MVA")){
  dat_merge <- merge(dat_gimms[nn == (period + 1)], 
                     dat_chrono[variable == eli_var],
                     by = c("site", "year"))
  dat_merge[, nyear := length(unique(year)), site]
  
  # subset to 1 year
  dat_cor <- dat_merge[period < 365/15, 
                       cor.test(value, mean_ndvi)[c("estimate", "p.value")],
                       .(site, nyear, xx_half_month, period)]
  
  dat_cor[, xx_month := xx_half_month/2 + 0.5]
  dat_cor[, yy_period_days := period * 15]
  
  
  l_gg_patch <- foreach(
    i_site = sort(unique(dat_cor$site))
  ) %do% {
    
    dat_i_cor <- dat_cor[site == i_site]
    dat_i_clim <- dat_plot_climate[site == i_site]
    
    gg_cor <- ggplot(dat_i_cor, aes(xx_month, yy_period_days, fill = estimate))+
      geom_raster()+
      geom_point(data = dat_i_cor[p.value < 0.1], size = 0.1, colour = "grey20")+
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
    
    gg_ppet <- ggplot(dat_i_clim, aes(month, PPET))+
      geom_line()+
      geom_hline(yintercept = 0, linetype = "dotted")+
      theme_bw()+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())+
      scale_x_continuous(breaks = 1:12 - 0.5,
                         labels = month.abb,
                         limits = c(0.5, 12.5),
                         expand = c(0,0))+
      xlab(NULL)+ylab(NULL)
    # gg_ppet
    
    
    gg_climgraph <- 
      ggplot(dat_i_clim, aes(month))+
      geom_line(aes(y = tmean), col ="red")+
      geom_line(aes(y = prec/2), col = "blue")+
      # scale_y_continuous(name = expression("Temperature ("~degree~"C)"), 
      #                    sec.axis = sec_axis(~.*2, name = "Precipitation (mm)"))+
      scale_y_continuous(sec.axis = sec_axis(~.*2))+
      theme_bw()+
      # theme(panel.grid = element_blank())+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())+
      scale_x_continuous(breaks = 1:12 - 0.5,
                         labels = month.abb,
                         limits = c(0.5, 12.5),
                         expand = c(0,0))+
      xlab(NULL)+ylab(NULL)
    # gg_climgraph
    
    gg_spei <- ggplot(dat_i_clim, aes(month, SPEI))+
      geom_line()+
      geom_hline(yintercept = 0, linetype = "dotted")+
      theme_bw()+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())+
      ylim(-1, 1)+
      scale_x_continuous(breaks = 1:12 - 0.5,
                         labels = month.abb,
                         limits = c(0.5, 12.5),
                         expand = c(0,0))+
      xlab(NULL)+ylab(NULL)
    gg_spei
    
    
    gg_cor / gg_ppet / gg_climgraph / gg_spei +
      plot_layout(heights = c(3, 1, 1, 1))
    
    
  }
  
  gg_out <- wrap_plots(l_gg_patch, nrow = 3, ncol = 7)
  
  
  ggsave(gg_out,
         file = paste0("fig/02-with-climate/GIMMS-NDVI-", eli_var, ".pdf"), 
         width = 28, height = 12, units = "in")
  
}
       
  

# ndvi -------------------------------------------------------------


dat_modis <- readRDS("data/modis-ndvi-02-moving.rds")
dat_modis$site %>% table

dat_modis[nn != (period + 1)]
dat_chrono[, .N, variable]

dat_vars <- data.table(col_name = colnames(dat_modis) %>% grep("mean", ., value = T))
dat_vars[, short_name := toupper(gsub("mean_", "", col_name))]
dat_vars[, long_name := c("NDVI (original)",
                          "NDVI (QC and filled)")]
dat_vars[, path_out := rep(c("fig/02-with-climate/MODIS-original/", "fig/02-with-climate/MODIS-filled/"), 1)]


dat_modis[, .(site, year)] %>% 
  unique %>% 
  merge(dat_climate) %>% 
  .[, lapply(.SD, mean), .(site, month)] -> dat_plot_climate


# plot all
for(i in 1:nrow(dat_vars)){
  
  for(eli_var in c("TRW", "Kh", "MVA")){
    
    dat_merge <- merge(dat_modis[nn == (period + 1)], 
                       dat_chrono[variable == eli_var],
                       by = c("site", "year"))
    dat_merge[, nyear := length(unique(year)), site]
    
    setnames(dat_merge, dat_vars[i, col_name], "cor_var")
    
    # subset to 1 year
    # dat_cor <- dat_merge[period < 365/16, 
    #                      .(pcor = cor(value, cor_var, use = "p")),
    #                      .(site, nyear, doy, period)]
    
    dat_cor <- dat_merge[period < 365/16, 
                         cor.test(value, cor_var)[c("estimate", "p.value")],
                         .(site, nyear, doy, period)]
    
    dat_cor[, xx_date := ymd("1999-01-01")]
    yday(dat_cor$xx_date) <- dat_cor$doy
    
    dat_cor[, yy_period_days := period * 16]
    
    l_gg_patch <- foreach(
      i_site = sort(unique(dat_cor$site))
    ) %do% {
      
      dat_i_cor <- dat_cor[site == i_site]
      dat_i_clim <- dat_plot_climate[site == i_site]
      
  
      gg_cor <- ggplot(dat_i_cor, aes(xx_date, yy_period_days, fill = estimate))+
        geom_raster()+
        geom_point(data = dat_i_cor[p.value < 0.1], size = 0.1, colour = "grey20")+
        scale_fill_gradient2(limits = c(-1,1))+
        coord_cartesian(expand = F)+
        scale_x_date(date_breaks = "month",
                     date_labels = "%b")+
        scale_y_continuous()+
        facet_wrap(~site + nyear, labeller = label_both)+
        theme_bw()+
        theme(legend.position = "none")+
        # ggtitle(paste0("Correlation between different averages of ", dat_vars[i, long_name],
        #                " and annual ", eli_var),
        #         "dots indicate p-value < 0.1; correlation is pearson")+
        # xlab("End of aggregation period")+
        # ylab(paste0("Number of days used to aggregate ", dat_vars[i, long_name]))+
        xlab(NULL)+ylab(NULL)
      # gg_cor
      
      gg_ppet <- ggplot(dat_i_clim, aes(month, PPET))+
        geom_line()+
        geom_hline(yintercept = 0, linetype = "dotted")+
        theme_bw()+
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())+
        scale_x_continuous(breaks = 1:12 - 0.5,
                           labels = month.abb,
                           limits = c(0.5, 12.5),
                           expand = c(0,0))+
        xlab(NULL)+ylab(NULL)
      # gg_ppet
      
      
      gg_climgraph <- 
        ggplot(dat_i_clim, aes(month))+
        geom_line(aes(y = tmean), col ="red")+
        geom_line(aes(y = prec/2), col = "blue")+
        # scale_y_continuous(name = expression("Temperature ("~degree~"C)"), 
        #                    sec.axis = sec_axis(~.*2, name = "Precipitation (mm)"))+
        scale_y_continuous(sec.axis = sec_axis(~.*2))+
        theme_bw()+
        # theme(panel.grid = element_blank())+
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())+
        scale_x_continuous(breaks = 1:12 - 0.5,
                           labels = month.abb,
                           limits = c(0.5, 12.5),
                           expand = c(0,0))+
        xlab(NULL)+ylab(NULL)
      # gg_climgraph
      
      gg_spei <- ggplot(dat_i_clim, aes(month, SPEI))+
        geom_line()+
        geom_hline(yintercept = 0, linetype = "dotted")+
        theme_bw()+
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())+
        ylim(-1, 1)+
        scale_x_continuous(breaks = 1:12 - 0.5,
                           labels = month.abb,
                           limits = c(0.5, 12.5),
                           expand = c(0,0))+
        xlab(NULL)+ylab(NULL)
      # gg_spei
      
      
      gg_cor / gg_ppet / gg_climgraph / gg_spei +
        plot_layout(heights = c(3, 1, 1, 1))
      
      
    }
    
    gg_out <- wrap_plots(l_gg_patch, nrow = 3, ncol = 7)
    
    
    ggsave(gg_out,
           file = paste0(dat_vars[i, path_out], "MODIS-", dat_vars[i, short_name], "-", eli_var, ".pdf"), 
           width = 28, height = 12, units = "in")
    
  }
  
  
  
}

# gpp lai fpar et -------------------------------------------------------------


dat_modis <- readRDS("data/gpp-lai-fpar-et-02-moving.rds")
dat_modis$site %>% table

dat_modis[nn != (period + 1)]
dat_chrono[, .N, variable]

dat_vars <- data.table(col_name = colnames(dat_modis) %>% grep("mean", ., value = T))
dat_vars[, short_name := toupper(gsub("mean_", "", col_name))]
dat_vars[, long_name := c("LAI (original)",
                          "LAI (QC and filled)",
                          "FPAR (original)",
                          "FPAR (QC and filled)",
                          "ET (original)",
                          "ET (QC and filled)",
                          "GPP (original)",
                          "GPP (QC and filled)")]
dat_vars[, path_out := rep(c("fig/02-with-climate/MODIS-original/", "fig/02-with-climate/MODIS-filled/"), 4)]



dat_modis[, .(site, year)] %>% 
  unique %>% 
  merge(dat_climate) %>% 
  .[, lapply(.SD, mean), .(site, month)] -> dat_plot_climate




# plot all
for(i in 1:nrow(dat_vars)){
  
  for(eli_var in c("TRW", "Kh", "MVA")){
    
    dat_merge <- merge(dat_modis[nn == (period + 1)], 
                       dat_chrono[variable == eli_var],
                       by = c("site", "year"))
    dat_merge[, nyear := length(unique(year)), site]
    
    setnames(dat_merge, dat_vars[i, col_name], "cor_var")
    
    # subset to 1 year
    # dat_cor <- dat_merge[period < 365/8, 
    #                      .(pcor = cor(value, cor_var, use = "p")),
    #                      .(site, nyear, doy, period)]
    
    dat_merge[, n_obs := sum(!is.na(value) & !is.na(cor_var)), .(site, nyear, doy, period)]

    dat_cor <- dat_merge[period < 365/8 & n_obs > 2, 
                         cor.test(value, cor_var)[c("estimate", "p.value")],
                         .(site, nyear, doy, period)]
    
    
    dat_cor[, xx_date := ymd("1999-01-01")]
    yday(dat_cor$xx_date) <- dat_cor$doy
    
    dat_cor[, yy_period_days := period * 8]
    
    
    l_gg_patch <- foreach(
      i_site = sort(unique(dat_cor$site))
    ) %do% {
      
      dat_i_cor <- dat_cor[site == i_site]
      dat_i_clim <- dat_plot_climate[site == i_site]
      
      gg_cor <- ggplot(dat_i_cor, aes(xx_date, yy_period_days, fill = estimate))+
        geom_raster()+
        geom_point(data = dat_i_cor[p.value < 0.1], size = 0.1, colour = "grey20")+
        scale_fill_gradient2(limits = c(-1,1))+
        coord_cartesian(expand = F)+
        scale_x_date(date_breaks = "month",
                     date_labels = "%b")+
        scale_y_continuous()+
        facet_wrap(~site + nyear, labeller = label_both)+
        theme_bw()+
        theme(legend.position = "none")+
        # ggtitle(paste0("Correlation between different averages of ", dat_vars[i, long_name],
        #                " and annual ", eli_var),
        #         "dots indicate p-value < 0.1; correlation is pearson")+
        # xlab("End of aggregation period")+
        # ylab(paste0("Number of days used to aggregate ", dat_vars[i, long_name]))+
        xlab(NULL)+ylab(NULL)
      # gg_cor
      
      gg_ppet <- ggplot(dat_i_clim, aes(month, PPET))+
        geom_line()+
        geom_hline(yintercept = 0, linetype = "dotted")+
        theme_bw()+
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())+
        scale_x_continuous(breaks = 1:12 - 0.5,
                           labels = month.abb,
                           limits = c(0.5, 12.5),
                           expand = c(0,0))+
        xlab(NULL)+ylab(NULL)
      # gg_ppet
      
      
      gg_climgraph <- 
        ggplot(dat_i_clim, aes(month))+
        geom_line(aes(y = tmean), col ="red")+
        geom_line(aes(y = prec/2), col = "blue")+
        # scale_y_continuous(name = expression("Temperature ("~degree~"C)"), 
        #                    sec.axis = sec_axis(~.*2, name = "Precipitation (mm)"))+
        scale_y_continuous(sec.axis = sec_axis(~.*2))+
        theme_bw()+
        # theme(panel.grid = element_blank())+
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())+
        scale_x_continuous(breaks = 1:12 - 0.5,
                           labels = month.abb,
                           limits = c(0.5, 12.5),
                           expand = c(0,0))+
        xlab(NULL)+ylab(NULL)
      # gg_climgraph
      
      gg_spei <- ggplot(dat_i_clim, aes(month, SPEI))+
        geom_line()+
        geom_hline(yintercept = 0, linetype = "dotted")+
        theme_bw()+
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())+
        ylim(-1, 1)+
        scale_x_continuous(breaks = 1:12 - 0.5,
                           labels = month.abb,
                           limits = c(0.5, 12.5),
                           expand = c(0,0))+
        xlab(NULL)+ylab(NULL)
      # gg_spei
      
      
      gg_cor / gg_ppet / gg_climgraph / gg_spei +
        plot_layout(heights = c(3, 1, 1, 1))
      
      
    }
    
    gg_out <- wrap_plots(l_gg_patch, nrow = 3, ncol = 7)
    
    
    ggsave(gg_out,
           file = paste0(dat_vars[i, path_out], "MODIS-", dat_vars[i, short_name], "-", eli_var, ".pdf"), 
           width = 28, height = 12, units = "in")
    
    
  }
  
  
  
}




# EOF ---------------------------------------------------------------------







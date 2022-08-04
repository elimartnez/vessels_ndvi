# some PCA stuff



library(data.table)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(scico)
library(lubridate)
library(forcats)
library(ggrepel)
library(patchwork)
library(foreach)
library(cluster)

dat_meta <- readRDS("data/chrono-00-meta.rds")



# settings ----------------------------------------------------------------

min_month <- 3 #incl
max_month <- 10 #incl
sites_to_remove <- c("CHL", "EIC", "BOI")
agg_period <- 120
# corr_fig_height <- 9
path_out <- "fig/13-check-pca/"


# prep data ---------------------------------------------------------------

dat_chrono <- readRDS("data/chrono-01.rds")
dat_chrono <- dat_chrono[! site %in% sites_to_remove]
dat_chrono[, table(variable, site)]


dat_gimms <- readRDS("data/gimms-ndvi-02-moving.rds")
dat_gimms <- dat_gimms[! site %in% sites_to_remove]
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

# eli_var <- "MVA"

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
  set.seed(1234)
  pca1 <- prcomp(t(mat_pca), rank. = 10, center = F, scale. = F)
  
  
  
  # variance explained
  dat_var_exp <- data.table(PC = paste0("PC", 1:length(pca1$sdev)),
                            prop_variance = pca1$sdev^2 / sum(pca1$sdev^2))
  dat_var_exp[, prop_variance_cumulative := cumsum(prop_variance)]
  dat_var_exp[, prop_variance := round(prop_variance, 3)]
  dat_var_exp[, prop_variance_cumulative := round(prop_variance_cumulative, 3)]
  fwrite(dat_var_exp, file = paste0(path_out, eli_var, "_pcaCorr-variance-explained.csv"))
  
  
  # biplot
  pdf(paste0(path_out, eli_var, "_pcaCorr-biplot.pdf"),
      width = 10, height = 10)
  biplot(pca1)
  dev.off()
  
  
  # pca - corr plot
  gg_pca_cor <-
  cbind(dat_wide[, .(xx_month, yy_plot)], pca1$rotation[, 1:3]) %>% 
    melt(id.vars = c("xx_month", "yy_plot")) %>% 
    .[, value_sc := scales::rescale(value, to = c(-1, 1))] %>% 
    ggplot(aes(xx_month, yy_plot, fill = value))+
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
         file = paste0(path_out, eli_var, "_pcaCorr-corrplot.pdf"), 
         width = 12, height = 4, units = "in")
  
  
  
  
  # pca by site -------------------------------------------------------------
  
  
  dat_cor %>% 
    dcast(site ~ yy_period_days + xx_month, value.var = "estimate") -> dat_pca
  
  sites <- dat_pca$site
  dat_pca[, site := NULL]
  
  mat_pca_site <- as.matrix(dat_pca)
  set.seed(1234)
  pca_site <- prcomp(t(mat_pca_site), rank. = 10, center = F, scale. = F)
  
  # variance explained
  dat_var_exp <- data.table(PC = paste0("PC", 1:length(pca_site$sdev)),
                            prop_variance = pca_site$sdev^2 / sum(pca_site$sdev^2))
  dat_var_exp[, prop_variance_cumulative := cumsum(prop_variance)]
  dat_var_exp[, prop_variance := round(prop_variance, 3)]
  dat_var_exp[, prop_variance_cumulative := round(prop_variance_cumulative, 3)]
  fwrite(dat_var_exp, file = paste0(path_out, eli_var, "_pcaSite-variance-explained.csv"))
  
  
  
  cbind(dat_wide[, .(xx_month, yy_plot)], pca1$rotation[, 1:3]) %>% 
    melt(id.vars = c("xx_month", "yy_plot"),
         value.name = "pc_val") -> dat_cor_pca
  
  dat_cor_site_pca <- foreach(
    i_pc = sort(unique(dat_cor_pca$variable)),
    .final = rbindlist
  ) %do% {
    
    dat_cor_pca[variable == i_pc] %>% 
      merge(dat_cor) 
    
  }
  
  
  dat_plot_tile <- dat_cor_site_pca[, .(cor(pc_val, estimate)), .(site, variable)]
  dat_plot_tile[variable == "PC1", fct_reorder(site, V1)] %>% 
    levels -> site_levels_ord
  dat_plot_tile[, site2 := factor(site, levels = site_levels_ord)] 
  
  
  
  
  gg_loading <- data.table(site = sites, pca_site$rotation[, 1:3]) %>% 
    merge(dat_meta) %>% 
    melt(measure.vars = paste0("PC", 1:3)) %>% 
    ggplot(aes(lon, lat, colour = value))+
    borders()+
    geom_point(size = 3)+
    geom_text_repel(aes(label = site))+
    coord_quickmap(xlim = range(dat_meta$lon), ylim = range(dat_meta$lat))+
    theme_bw()+
    scale_color_gradient2()+
    facet_wrap(~variable, nrow = 1)+
    ggtitle("PC loading by site")
  
  gg_cor_pc <- 
    dat_plot_tile %>% 
    merge(dat_meta) %>% 
    ggplot(aes(lon, lat, colour = V1))+
    borders()+
    geom_point(size = 3)+
    geom_text_repel(aes(label = site))+
    coord_quickmap(xlim = range(dat_meta$lon), ylim = range(dat_meta$lat))+
    theme_bw()+
    scale_color_gradient2()+
    facet_wrap(~variable, nrow = 1)+
    ggtitle("Correlation between PC and site patterns")
  
  
  gg_out <- gg_loading / gg_cor_pc
  
  ggsave(paste0(path_out, eli_var, "pcaSite_map-comparison-PC.png"),
         gg_out,
         width = 12, height = 6)
  

# correlate PCs to climate ------------------------------------------------

  
  dat_plot_pca <- data.table(pca_site$rotation[, 1:3])
  dat_plot_pca[, site := sites]
  
  for(i_climvar in c("PPET", "tmean", "prec")){
    
    merge(dat_plot_clim,
          dat_plot_pca,
          by = "site") %>% 
      melt(id.vars = c("site", "month", i_climvar), measure.vars = paste0("PC", 1:3)) %>% 
      merge(dat_meta, by = "site") -> dat_pca_clim
    
    dat_pca_clim[, month_f := factor(month.abb[month], levels = month.abb)]
    
    # add sig and summary
    dat_pca_clim_summ <- dat_pca_clim[,
                 .(corr = cor(get(i_climvar), value),
                   corr_pval = cor.test(get(i_climvar), value)$p.value,
                   lm_rsq = summary(lm(value ~ get(i_climvar)))$r.squared,
                   lm_pval = broom::glance(lm(value ~ get(i_climvar)))$p.value,
                   xx = min(get(i_climvar)),
                   yy = max(value)
                 ),
                 .(month_f, variable)]
    
    yy_max <- max(dat_pca_clim$value)*1.3
    # dat_pca_clim_summ[, lbl := sprintf("cor = %.2f | pval = %.2f", corr, corr_pval)]
    dat_pca_clim_summ[, lbl := sprintf("rsq = %.2f\npval = %.2f", lm_rsq, lm_pval)]
    
    gg_pca_clim <-
    dat_pca_clim[variable %in% paste0("PC", 1:3) & month %in% c(1:12)] %>% 
      ggplot(aes_string(i_climvar, "value"))+
      # geom_hline(yintercept = 0)+
      # geom_vline(xintercept = 0)+
      geom_smooth(method = lm, formula = y ~ x, se = F, colour = "black", linetype = "dashed")+
      geom_point(aes(colour = species))+ 
      # geom_text(aes(label = site, colour = species), size = 3)+
      geom_text_repel(aes(label = site, colour = species), size = 2)+
      geom_text(inherit.aes = F, data = dat_pca_clim_summ,
                aes(x = xx, y = yy_max, label = lbl), 
                hjust = 0, vjust = 1, size = 3)+
      scale_color_brewer(palette = "Set1")+
      facet_grid(variable ~ month_f, scales = "free_x", space = "free_x")+
      theme_bw()+
      ylab("PC loading")
    # gg_pca_clim
    
    
    ggsave(gg_pca_clim,
           file = paste0(path_out, eli_var, "_pcaSite-scatter-clim_", i_climvar, ".pdf"),
           width = 20, height = 6)
    
    
    
  }
  

# kmeans ------------------------------------------------------------------


# ** chrono ------------------------------------------------------------------

  

  dat_chrono[variable == eli_var] %>% 
    dcast(year ~ site, value.var = "value") %>% 
    na.omit -> dat_chrono_km
  
  mat_sub <- as.matrix(dat_chrono_km)[, -1]
  mat_sub_sc <- scale(mat_sub)
  mat_clust_obs <- t(mat_sub_sc)
  
  set.seed(1234)
  dat_clust_obs <- foreach(
    kk = 2:5,
    .final = rbindlist
  ) %do% {
    
    km_fit <- kmeans(mat_clust_obs, kk)
    stn_clust <- km_fit$cluster
    
    km_sil <- silhouette(km_fit$cluster, dist(mat_clust_obs))
    
    data.table(kk = kk,
               site = names(km_fit$cluster), 
               cluster = km_sil[, 1],
               # cluster_neigh = km_sil[, 2],
               sil_width = km_sil[, 3])
    
  }

  gg_clust <- dat_clust_obs %>% 
    merge(dat_meta) %>% 
    ggplot(aes(lon, lat, colour = as.factor(cluster)))+
    borders()+
    geom_point(size = 3)+
    geom_text_repel(aes(label = site))+
    coord_quickmap(xlim = range(dat_meta$lon), ylim = range(dat_meta$lat))+
    theme_bw()+
    scale_color_brewer(palette = "Set1")+
    facet_grid(. ~ kk)+
    xlab(NULL)+ylab(NULL)
  
  gg_sil <- dat_clust_obs %>% 
    merge(dat_meta) %>% 
    ggplot(aes(lon, lat, colour = sil_width))+
    borders()+
    geom_point(size = 3)+
    geom_text_repel(aes(label = site))+
    coord_quickmap(xlim = range(dat_meta$lon), ylim = range(dat_meta$lat))+
    theme_bw()+
    scale_color_viridis_c(direction = -1)+
    facet_grid(. ~ kk)+
    xlab(NULL)+ylab(NULL)
  
  gg_out <- gg_clust / gg_sil
  
  ggsave(gg_out,
         filename = paste0(path_out, eli_var, "_kmeans_chrono.png"),
         width = 15, height = 6)
  
  
  # ** corr patterns ------------------------------------------------------------------
  

  mat_clust_corr <- t(mat_pca)
  
  set.seed(1234)
  dat_clust_corr <- foreach(
    kk = 2:5,
    .final = rbindlist
  ) %do% {
    
    km_fit <- kmeans(mat_clust_corr, kk)
    stn_clust <- km_fit$cluster
    
    km_sil <- silhouette(km_fit$cluster, dist(mat_clust_corr))
    
    data.table(kk = kk,
               site = names(km_fit$cluster), 
               cluster = km_sil[, 1],
               # cluster_neigh = km_sil[, 2],
               sil_width = km_sil[, 3])
    
  }
  
  gg_clust <- dat_clust_corr %>% 
    merge(dat_meta) %>% 
    ggplot(aes(lon, lat, colour = as.factor(cluster)))+
    borders()+
    geom_point(size = 3)+
    geom_text_repel(aes(label = site))+
    coord_quickmap(xlim = range(dat_meta$lon), ylim = range(dat_meta$lat))+
    theme_bw()+
    scale_color_brewer(palette = "Set1")+
    facet_grid(. ~ kk)+
    xlab(NULL)+ylab(NULL)
  
  gg_sil <- dat_clust_corr %>% 
    merge(dat_meta) %>% 
    ggplot(aes(lon, lat, colour = sil_width))+
    borders()+
    geom_point(size = 3)+
    geom_text_repel(aes(label = site))+
    coord_quickmap(xlim = range(dat_meta$lon), ylim = range(dat_meta$lat))+
    theme_bw()+
    scale_color_viridis_c(direction = -1)+
    facet_grid(. ~ kk)+
    xlab(NULL)+ylab(NULL)
  
  gg_out <- gg_clust / gg_sil
  
  ggsave(gg_out,
         filename = paste0(path_out, eli_var, "_kmeans_corr_patterns.png"),
         width = 15, height = 6)
  
  
  
  # ** site pca of corr patterns (maybe not) ------------------------------------------------------------------
  
  # 
  # mat_clust_corr_pca <- pca_site$rotation[, 1:5]
  # 
  # set.seed(1234)
  # dat_clust_corr_pca <- foreach(
  #   kk = 2:5,
  #   .final = rbindlist
  # ) %do% {
  #   
  #   km_fit <- kmeans(mat_clust_corr_pca, kk)
  #   stn_clust <- km_fit$cluster
  #   
  #   km_sil <- silhouette(km_fit$cluster, dist(mat_clust_corr_pca))
  #   
  #   data.table(kk = kk,
  #              site = sites, 
  #              cluster = km_sil[, 1],
  #              # cluster_neigh = km_sil[, 2],
  #              sil_width = km_sil[, 3])
  #   
  # }
  # 
  # gg_clust <- dat_clust_corr_pca %>% 
  #   merge(dat_meta) %>% 
  #   ggplot(aes(lon, lat, colour = as.factor(cluster)))+
  #   borders()+
  #   geom_point(size = 3)+
  #   geom_text_repel(aes(label = site))+
  #   coord_quickmap(xlim = range(dat_meta$lon), ylim = range(dat_meta$lat))+
  #   theme_bw()+
  #   scale_color_brewer(palette = "Set1")+
  #   facet_grid(. ~ kk)+
  #   xlab(NULL)+ylab(NULL)
  # 
  # gg_sil <- dat_clust_corr_pca %>% 
  #   merge(dat_meta) %>% 
  #   ggplot(aes(lon, lat, colour = sil_width))+
  #   borders()+
  #   geom_point(size = 3)+
  #   geom_text_repel(aes(label = site))+
  #   coord_quickmap(xlim = range(dat_meta$lon), ylim = range(dat_meta$lat))+
  #   theme_bw()+
  #   scale_color_viridis_c(direction = -1)+
  #   facet_grid(. ~ kk)+
  #   xlab(NULL)+ylab(NULL)
  # 
  # gg_out <- gg_clust / gg_sil
  # 
  # ggsave(gg_out,
  #        filename = paste0(path_out, eli_var, "_kmeans_site_pca_corr.png"),
  #        width = 15, height = 6)
  # 
  # 
  
  
}



# EOF ---------------------------------------------------------------------




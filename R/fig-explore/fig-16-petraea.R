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
sites_petraea <- dat_meta[species == "Quercus petraea", site]


# settings ----------------------------------------------------------------

min_month <- 3 #incl
max_month <- 10 #incl
sites_to_remove <- c("CHL", "EIC", "BOI")
agg_period <- 120
# corr_fig_height <- 9
path_out <- "fig/17-petraea/"


# prep data ---------------------------------------------------------------

dat_chrono <- readRDS("data/chrono-01.rds")
dat_chrono <- dat_chrono[! site %in% sites_to_remove & site %in% sites_petraea]
dat_chrono[, table(variable, site)]


dat_gimms <- readRDS("data/gimms-ndvi-04-detrend-moving.rds")
dat_gimms <- dat_gimms[! site %in% sites_to_remove & site %in% sites_petraea]
dat_gimms$site %>% table


dat_clim <- readRDS("data/climate-01.rds")
dat_gimms[, .(site, year)] %>% 
  unique %>% 
  merge(unique(dat_chrono[, .(site, year)])) %>% 
  merge(dat_clim) %>% 
  .[, lapply(.SD, mean), .(site, month)] -> dat_plot_clim

setnames(dat_plot_clim, "PPET", "PPET_diff")
dat_plot_clim[, PPET_ratio := prec / PET]

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
                       cor.test(value, mean_ndvi_detrend)[c("estimate", "p.value")],
                       .(site, nyear, xx_half_month, period)]
  
  dat_cor[, xx_month := xx_half_month/2 + 0.5]
  dat_cor[, yy_period_days := period * 15]
  dat_cor[, yy_plot := yy_period_days + 15/2]
  
  # subset months
  dat_cor <- dat_cor[xx_month >= min_month & xx_month < (max_month + 1)]
  
  dat_cor[, .N, .(site)]
  
  

  
  # reorder sites
  dat_cor[, site_f := factor(site)]
  dat_cor[, avg_estimate_jja := mean(estimate[xx_month >= 6 & xx_month < 9]), .(site)]
  dat_cor[, avg_estimate := mean(estimate), .(site)]
  
  
  
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
         file = paste0(path_out, eli_var, "_corr-plot_order-alphabet.pdf"), 
         width = 18, height = 9, units = "in")
  
  
  dat_cor[, site_f := fct_reorder(site, avg_estimate)]
  
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
  
  
  
  # pca data
  cbind(dat_wide[, .(xx_month, yy_plot)], pca1$rotation[, 1:3]) %>% 
    melt(id.vars = c("xx_month", "yy_plot")) %>% 
    .[, value_sc := scales::rescale(value, to = c(-1, 1))] -> dat_pca_cor
  
  
  dat_pca_site_wide <- data.table(site = colnames(mat_pca), pca1$x[, 1:3])
  
  dat_pca_site_wide %>% 
    melt(id.vars = c("site")) %>% 
    .[, value_sc := scales::rescale(value, to = c(-1, 1)), variable] -> dat_pca_site
  
  
  
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
   dat_pca_cor %>% 
    ggplot(aes(xx_month, yy_plot, fill = value))+
    geom_raster()+
    # scale_fill_gradient2()+
    scale_fill_scico(palette = "roma", rescaler = scales::rescale_mid)+
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
  
  
  # coefs scatter
  gg_pca_coef <- dat_pca_site %>% 
    ggplot(aes(value, fct_rev(site)))+
    geom_vline(xintercept = 0, linetype = "dashed")+
    geom_point()+
    facet_wrap(~variable)+
    theme_bw()+
    ylab("Coefficients per site and PC")+
    xlab(NULL)
  
  
  
  # coefs map
  gg_pca_coef_map <- dat_pca_site %>% 
    merge(dat_meta) %>% 
    ggplot(aes(lon, lat, colour = value))+
    borders()+
    geom_point(size = 3)+
    geom_text_repel(aes(label = site))+
    coord_quickmap(xlim = range(dat_meta$lon), ylim = range(dat_meta$lat))+
    theme_bw()+
    scale_color_scico(palette = "corkO", rescaler = scales::rescale_mid)+
    facet_wrap(~variable)+
    xlab(NULL)+
    ylab(NULL)+
    ggtitle("Coef displayed in map")
  
  
  
  gg_pca_comb <- gg_pca_cor / gg_pca_coef / gg_pca_coef_map
  
  
  ggsave(gg_pca_comb,
         file = paste0(path_out, eli_var, "_pcaCorr-corrplot-with-coef.pdf"), 
         width = 12, height = 12, units = "in")
  
  

# correlate PC coef to climate ------------------------------------------------


  for(i_climvar in c("PPET_diff", "PPET_ratio", "tmean", "prec")){
    
    merge(dat_plot_clim,
          dat_pca_site_wide,
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
                aes(x = xx, y = yy, label = lbl), 
                hjust = 0, vjust = 1, size = 3)+
      scale_color_brewer(palette = "Set1")+
      facet_grid(variable ~ month_f, scales = "free")+
      theme_bw()+
      ylab("PC site coefficient")
    # gg_pca_clim
    
    
    ggsave(gg_pca_clim,
           file = paste0(path_out, eli_var, "_pcaCorr-scatter-coef-clim_", i_climvar, ".pdf"),
           width = 20, height = 6)
    
    
    
  }
  

# decompose site pattern by PC --------------------------------------------

  dat_pca_site[, .(site, pc = variable, pc_coef = value)] %>% 
    split(by = c("site")) %>% 
    lapply(function(x) merge(x, dat_pca_cor[, .(xx_month, yy_plot, pc = variable, pc_cor = value)])) %>% 
    rbindlist() -> dat_pca_decompose
  
  dat_pca_decompose[, pc_cor_coef := pc_cor*pc_coef]
  
  dat_pca_decompose %>% 
    dcast(site + xx_month + yy_plot ~ pc, value.var = "pc_cor_coef") %>% 
    merge(dat_cor[, .(site, xx_month, yy_plot, obs_cor = estimate)]) -> dat_plot_decompose
  
  dat_plot_decompose[, rest := obs_cor - PC1 - PC2 - PC3]
  
  dat_plot_decompose %>% 
    melt(id.vars = c("site", "xx_month", "yy_plot")) -> dat_plot_decompose2
  
  dat_plot_decompose2[, variable_f := factor(variable,
                                             levels = c("obs_cor", "PC1", "PC2", "PC3", "rest"))]
  
  pdf(width = 14, height = 4,
      file = paste0(path_out, eli_var, "_pcaCorr-decompose.pdf"))
  for(i_site in sort(unique(dat_plot_decompose2$site))){
    
    gg <- dat_plot_decompose2[site == i_site] %>% 
      ggplot(aes(xx_month, yy_plot, fill = value))+
      geom_raster()+
      scale_fill_gradient2("cor", limits = c(-1,1))+
      coord_cartesian(expand = F)+
      scale_x_continuous(breaks = c(min_month : max_month) - 0.25,
                         labels = month.abb[min_month : max_month])+
      scale_y_continuous(breaks = c(0, 30, 60, 90, 120))+
      facet_grid(~ variable_f)+
      theme_bw()+
      ggtitle(i_site)+
      xlab("End of aggregation period")+
      ylab("Number of days used to aggregate NDVI")
    print(gg)
    
  }
  dev.off()
  
  
    
  

# # kmeans ------------------------------------------------------------------
# 
# 
# # ** chrono ------------------------------------------------------------------
# 
#   
# 
#   dat_chrono[variable == eli_var] %>% 
#     dcast(year ~ site, value.var = "value") %>% 
#     na.omit -> dat_chrono_km
#   
#   mat_sub <- as.matrix(dat_chrono_km)[, -1]
#   mat_sub_sc <- scale(mat_sub)
#   mat_clust_obs <- t(mat_sub_sc)
#   
#   set.seed(1234)
#   dat_clust_obs <- foreach(
#     kk = 2:5,
#     .final = rbindlist
#   ) %do% {
#     
#     km_fit <- kmeans(mat_clust_obs, kk)
#     stn_clust <- km_fit$cluster
#     
#     km_sil <- silhouette(km_fit$cluster, dist(mat_clust_obs))
#     
#     data.table(kk = kk,
#                site = names(km_fit$cluster), 
#                cluster = km_sil[, 1],
#                # cluster_neigh = km_sil[, 2],
#                sil_width = km_sil[, 3])
#     
#   }
# 
#   gg_clust <- dat_clust_obs %>% 
#     merge(dat_meta) %>% 
#     ggplot(aes(lon, lat, colour = as.factor(cluster)))+
#     borders()+
#     geom_point(size = 3)+
#     geom_text_repel(aes(label = site))+
#     coord_quickmap(xlim = range(dat_meta$lon), ylim = range(dat_meta$lat))+
#     theme_bw()+
#     scale_color_brewer(palette = "Set1")+
#     facet_grid(. ~ kk)+
#     xlab(NULL)+ylab(NULL)
#   
#   gg_sil <- dat_clust_obs %>% 
#     merge(dat_meta) %>% 
#     ggplot(aes(lon, lat, colour = sil_width))+
#     borders()+
#     geom_point(size = 3)+
#     geom_text_repel(aes(label = site))+
#     coord_quickmap(xlim = range(dat_meta$lon), ylim = range(dat_meta$lat))+
#     theme_bw()+
#     scale_color_viridis_c(direction = -1)+
#     facet_grid(. ~ kk)+
#     xlab(NULL)+ylab(NULL)
#   
#   gg_out <- gg_clust / gg_sil
#   
#   ggsave(gg_out,
#          filename = paste0(path_out, eli_var, "_kmeans_chrono.png"),
#          width = 15, height = 6)
#   
#   
#   # ** corr patterns ------------------------------------------------------------------
#   
# 
#   mat_clust_corr <- t(mat_pca)
#   
#   set.seed(1234)
#   dat_clust_corr <- foreach(
#     kk = 2:5,
#     .final = rbindlist
#   ) %do% {
#     
#     km_fit <- kmeans(mat_clust_corr, kk)
#     stn_clust <- km_fit$cluster
#     
#     km_sil <- silhouette(km_fit$cluster, dist(mat_clust_corr))
#     
#     data.table(kk = kk,
#                site = names(km_fit$cluster), 
#                cluster = km_sil[, 1],
#                # cluster_neigh = km_sil[, 2],
#                sil_width = km_sil[, 3])
#     
#   }
#   
#   gg_clust <- dat_clust_corr %>% 
#     merge(dat_meta) %>% 
#     ggplot(aes(lon, lat, colour = as.factor(cluster)))+
#     borders()+
#     geom_point(size = 3)+
#     geom_text_repel(aes(label = site))+
#     coord_quickmap(xlim = range(dat_meta$lon), ylim = range(dat_meta$lat))+
#     theme_bw()+
#     scale_color_brewer(palette = "Set1")+
#     facet_grid(. ~ kk)+
#     xlab(NULL)+ylab(NULL)
#   
#   gg_sil <- dat_clust_corr %>% 
#     merge(dat_meta) %>% 
#     ggplot(aes(lon, lat, colour = sil_width))+
#     borders()+
#     geom_point(size = 3)+
#     geom_text_repel(aes(label = site))+
#     coord_quickmap(xlim = range(dat_meta$lon), ylim = range(dat_meta$lat))+
#     theme_bw()+
#     scale_color_viridis_c(direction = -1)+
#     facet_grid(. ~ kk)+
#     xlab(NULL)+ylab(NULL)
#   
#   gg_out <- gg_clust / gg_sil
#   
#   ggsave(gg_out,
#          filename = paste0(path_out, eli_var, "_kmeans_corr_patterns.png"),
#          width = 15, height = 6)
#   
#   
#   
#   # plot corr plots kmeans by cluster
#   setorder(dat_clust_corr, kk, cluster)
#   
#   pdf(paste0(path_out, eli_var, "_kmeans_corr_patterns_grouping.pdf"),
#       width = 12, height = 12)
#   
#   for(i_kk in 2:5){
#     
#     l_gg <- dat_clust_corr[kk == i_kk] %>% split(by = "cluster") %>% 
#       lapply(function(dat){
#         
#         tit <- dat[1, paste0("Cluster: ", cluster)]
#         
#         dat %>% 
#           merge(dat_cor) %>% 
#           ggplot(aes(xx_month, yy_plot, fill = estimate))+
#           geom_raster()+
#           scale_fill_gradient2("cor", limits = c(-1,1), guide = F)+
#           coord_cartesian(expand = F)+
#           scale_x_continuous(breaks = c(min_month : max_month) - 0.25,
#                              labels = month.abb[min_month : max_month])+
#           scale_y_continuous(breaks = c(0, 30, 60, 90, 120))+
#           facet_wrap(~site, ncol = 5)+
#           theme_bw()+
#           ggtitle(tit)+
#           xlab(NULL)+
#           ylab(NULL)
#         
#       })
#     
#     dat_clust_corr[kk == i_kk, .N, cluster] %>% .[, N] -> n_clust
#     heights_clust <- ceiling(n_clust / 5)
#     
#     wrap_plots(l_gg, ncol = 1, heights = heights_clust / sum(heights_clust)) %>% print
#     
#   }
# 
#   dev.off()
#   
#   
}



# EOF ---------------------------------------------------------------------




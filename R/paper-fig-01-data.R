# paper figure: create data




library(data.table)
library(magrittr)
library(ggplot2)
library(lubridate)
library(forcats)
library(foreach)

# library(scico)
# library(ggrepel)
# library(patchwork)
# library(cluster)
# library(lemon)

dat_meta <- readRDS("data/chrono-00-meta.rds")



# settings ----------------------------------------------------------------

min_month <- 3 #incl
max_month <- 10 #incl
sites_to_remove <- c("CHL", "EIC", "BOI")
agg_period <- 120
# corr_fig_height <- 9
# path_out <- "fig/20-cleaning-up/"
n_pcs <- 4


# prep data ---------------------------------------------------------------

dat_chrono <- readRDS("data/chrono-02-kh1.rds")
dat_chrono <- dat_chrono[! site %in% sites_to_remove]
dat_chrono[, table(variable, site)]

# remove outlier CUG Kh, TRW, MVA
dat_chrono <- dat_chrono[!(site == "CUG" & variable == "Kh" & year == 2006)]
dat_chrono <- dat_chrono[!(site == "CUG" & variable == "MVA" & year == 1986)]



dat_gimms <- readRDS("data/gimms-ndvi-04-detrend-moving.rds")
dat_gimms <- dat_gimms[! site %in% sites_to_remove]
dat_gimms$site %>% table


dat_clim <- readRDS("data/climate-01.rds")
dat_gimms[, .(site, year)] %>% 
  unique %>% 
  merge(unique(dat_chrono[, .(site, year)])) %>% 
  merge(dat_clim) %>% 
  .[, lapply(.SD, mean), .(site, month)] -> dat_plot_clim

setnames(dat_plot_clim, "PPET", "PPET_diff")
dat_plot_clim[, PPET_ratio := prec / PET]

# saveRDS(dat_plot_clim, file = paste0(path_out, "data-avg-climate.rds"))

dat_plot_clim <- dat_plot_clim[month >= min_month & month <= max_month]

dat_meta <- readRDS("data/chrono-00-meta.rds")




# loop over tree ring measures -----------------------------------------------------

# eli_var <- "MVA"

l_cor <- list()

l_pca_varexp <- list()
l_pca_cor <- list()
l_pca_site <- list()

l_pca_clim <- list()
l_pca_clim_summ <- list()

l_decompose <- list()

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
  
  # save
  l_cor[[eli_var]] <- dat_cor
  
  
  
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

 
  # save
  l_pca_varexp[[eli_var]] <- dat_var_exp
    
  
  # fix sign (adopted from princomp)
  pca1_rot <- pca1$rotation[, 1:n_pcs]
  pca1_x <- pca1$x[, 1:n_pcs]
  
  # sign of first loading (coef)
  # signs <- ifelse(pca1_x[1, ] < 0, -1, 1) 
  
  # average of corr_plot values for first quarter (or so)
  signs <- ifelse(colMeans(pca1_rot[1:floor(nrow(pca1_rot)/8), ]) < 0, -1, 1)
  
  pca1_rot <- sweep(pca1_rot, 2L, signs, "*")
  pca1_x <- sweep(pca1_x, 2L, signs, "*")
  
  
  # pca data
  cbind(dat_wide[, .(xx_month, yy_plot)], pca1_rot) %>% 
    melt(id.vars = c("xx_month", "yy_plot")) %>% 
    .[, value_sc := scales::rescale(value, to = c(-1, 1))] -> dat_pca_cor
  
  
  dat_pca_site_wide <- data.table(site = colnames(mat_pca), pca1_x)
  
  dat_pca_site_wide %>% 
    melt(id.vars = c("site")) %>% 
    .[, value_sc := scales::rescale(value, to = c(-1, 1)), variable] -> dat_pca_site
  
  
  
  
  # save
  l_pca_cor[[eli_var]] <- dat_pca_cor
  l_pca_site[[eli_var]] <- dat_pca_site
  
  
  # correlate PC coef to climate ------------------------------------------------
  
  # i_climvar <- "prec"
  
  l_zz_clim <- list()
  l_zz_clim_summ <- list()
  
  for(i_climvar in c("PPET_ratio", "tmean", "prec")){
    
    merge(dat_plot_clim,
          dat_pca_site_wide,
          by = "site") %>% 
      melt(id.vars = c("site", "month", i_climvar), measure.vars = paste0("PC", 1:n_pcs)) %>% 
      merge(dat_meta, by = "site") -> dat_pca_clim
    
    
    dat_pca_clim[, month_f := factor(month.abb[month], levels = month.abb)]
    
    # add sig and summary
    dat_pca_clim_summ <- dat_pca_clim[,
                                      .(corr = cor(get(i_climvar), value, method = "p"),
                                        corr_pval = cor.test(get(i_climvar), value, method = "p")$p.value,
                                        lm_rsq = summary(lm(value ~ get(i_climvar)))$r.squared,
                                        lm_pval = broom::glance(lm(value ~ get(i_climvar)))$p.value,
                                        xx = min(get(i_climvar)),
                                        yy = max(value)
                                      ),
                                      .(month_f, variable)]
    
    
    # save
    setnames(dat_pca_clim, i_climvar, "climvar_value")
    l_zz_clim[[i_climvar]] <- dat_pca_clim
    l_zz_clim_summ[[i_climvar]] <- dat_pca_clim_summ
    
    
    
  }
  
  l_pca_clim[[eli_var]] <- rbindlist(l_zz_clim, idcol = "climvar")
  l_pca_clim_summ[[eli_var]] <- rbindlist(l_zz_clim_summ, idcol = "climvar")
  
  
  # decompose site pattern by PC --------------------------------------------
  
  dat_pca_site[, .(site, pc = variable, pc_coef = value)] %>% 
    split(by = c("site")) %>% 
    lapply(function(x) merge(x, dat_pca_cor[, .(xx_month, yy_plot, pc = variable, pc_cor = value)])) %>% 
    rbindlist() -> dat_pca_decompose
  
  dat_pca_decompose[, pc_cor_coef := pc_cor*pc_coef]
  
  dat_pca_decompose %>% 
    dcast(site + xx_month + yy_plot ~ pc, value.var = "pc_cor_coef") %>% 
    merge(dat_cor[, .(site, xx_month, yy_plot, obs_cor = estimate)]) -> dat_plot_decompose
  
  sum_all_pc <- rowSums(dat_plot_decompose[, paste0("PC", 1:n_pcs), with = F])
  dat_plot_decompose[, rest := obs_cor - sum_all_pc]
  
  dat_plot_decompose %>% 
    melt(id.vars = c("site", "xx_month", "yy_plot")) -> dat_plot_decompose2
  
  dat_plot_decompose2[, variable_f := factor(variable,
                                             levels = c("obs_cor", paste0("PC", 1:n_pcs), "rest"))]
  
  # save
  
  l_decompose[[eli_var]] <- dat_plot_decompose2
  
  
  
  
  
}



# combine and save


dat_cor <- rbindlist(l_cor, idcol = "tree_var")
dat_pca_varexp <- rbindlist(l_pca_varexp, idcol = "tree_var")
dat_pca_cor <- rbindlist(l_pca_cor, idcol = "tree_var")
dat_pca_site <- rbindlist(l_pca_site, idcol = "tree_var")
dat_pca_clim <- rbindlist(l_pca_clim, idcol = "tree_var")
dat_pca_clim_summ <- rbindlist(l_pca_clim_summ, idcol = "tree_var")
dat_decompose <- rbindlist(l_decompose, idcol = "tree_var")


save(dat_cor, 
     dat_pca_varexp, dat_pca_cor, dat_pca_site,
     dat_pca_clim, dat_pca_clim_summ,
     dat_decompose,
     min_month, max_month,
     file = "data/paper-fig-01.rda")



# EOF ---------------------------------------------------------------------







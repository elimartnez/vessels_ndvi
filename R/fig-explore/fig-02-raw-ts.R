# plot raw time series of MODIS data

library(data.table)
library(magrittr)
library(ggplot2)
library(purrr)




# gpp - check removed & filled values -----------------------------------------------


dat <- readRDS("data/gpp-lai-fpar-et-01.rds")

all_sites <- sort(unique(dat$site))

pdf(file = "fig/raw-ts/gpp-qc-fill.pdf",
    width = 14, height = 8)
walk(all_sites, function(i_site){
  
  pch_size <- 1.5
  
  gg <- dat[site == i_site] %>% 
    ggplot(aes(date))+
    geom_point(aes(y = gpp, colour = "1 original"), size = pch_size)+
    geom_point(aes(y = gpp_qc, colour = "2 after QC"), size = pch_size/4)+
    geom_point(data = dat[site == i_site & is.na(gpp_qc)],
               aes(y = gpp_f, colour = "3 filled"), size = pch_size)+
    facet_wrap(~year(date), scales = "free_x")+
    scale_x_date(date_breaks = "3 months", date_labels = "%b")+
    scale_shape_manual(values = c(19, 1))+
    theme_bw()+
    ggtitle(i_site)
  
  print(gg)
  
  
})
dev.off()




# lai - check removed & filled values -----------------------------------------------


dat <- readRDS("data/gpp-lai-fpar-et-01.rds")

all_sites <- sort(unique(dat$site))

pdf(file = "fig/raw-ts/lai-qc-fill.pdf",
    width = 14, height = 8)
walk(all_sites, function(i_site){
  
  pch_size <- 1.5
  
  gg <- dat[site == i_site] %>% 
    ggplot(aes(date))+
    geom_point(aes(y = lai, colour = "1 original"), size = pch_size)+
    geom_point(aes(y = lai_qc, colour = "2 after QC"), size = pch_size/4)+
    geom_point(data = dat[site == i_site & is.na(lai_qc)],
               aes(y = lai_f, colour = "3 filled"), size = pch_size)+
    facet_wrap(~year(date), scales = "free_x")+
    scale_x_date(date_breaks = "3 months", date_labels = "%b")+
    scale_shape_manual(values = c(19, 1))+
    theme_bw()+
    ggtitle(i_site)
  
  print(gg)
  
  
})
dev.off()


# fpar - check removed & filled values -----------------------------------------------


dat <- readRDS("data/gpp-lai-fpar-et-01.rds")

all_sites <- sort(unique(dat$site))

pdf(file = "fig/raw-ts/fpar-qc-fill.pdf",
    width = 14, height = 8)
walk(all_sites, function(i_site){
  
  pch_size <- 1.5
  
  gg <- dat[site == i_site] %>% 
    ggplot(aes(date))+
    geom_point(aes(y = fpar, colour = "1 original"), size = pch_size)+
    geom_point(aes(y = fpar_qc, colour = "2 after QC"), size = pch_size/4)+
    geom_point(data = dat[site == i_site & is.na(fpar_qc)],
               aes(y = fpar_f, colour = "3 filled"), size = pch_size)+
    facet_wrap(~year(date), scales = "free_x")+
    scale_x_date(date_breaks = "3 months", date_labels = "%b")+
    scale_shape_manual(values = c(19, 1))+
    theme_bw()+
    ggtitle(i_site)
  
  print(gg)
  
  
})
dev.off()


# et - check removed & filled values -----------------------------------------------


dat <- readRDS("data/gpp-lai-fpar-et-01.rds")

all_sites <- sort(unique(dat$site))

pdf(file = "fig/raw-ts/et-qc-fill.pdf",
    width = 14, height = 8)
walk(all_sites, function(i_site){
  
  pch_size <- 1.5
  
  gg <- dat[site == i_site] %>% 
    ggplot(aes(date))+
    geom_point(aes(y = et, colour = "1 original"), size = pch_size)+
    geom_point(aes(y = et_qc, colour = "2 after QC"), size = pch_size/4)+
    geom_point(data = dat[site == i_site & is.na(et_qc)],
               aes(y = et_f, colour = "3 filled"), size = pch_size)+
    facet_wrap(~year(date), scales = "free_x")+
    scale_x_date(date_breaks = "3 months", date_labels = "%b")+
    scale_shape_manual(values = c(19, 1))+
    theme_bw()+
    ggtitle(i_site)
  
  print(gg)
  
  
})
dev.off()



# ndvi - check removed & filled values -----------------------------------------------


dat <- readRDS("data/modis-ndvi-01.rds")

all_sites <- sort(unique(dat$site))

pdf(file = "fig/raw-ts/ndvi-qc-fill.pdf",
    width = 14, height = 8)
walk(all_sites, function(i_site){
  
  pch_size <- 1.5
  
  gg <- dat[site == i_site] %>% 
    ggplot(aes(date))+
    geom_point(aes(y = ndvi, colour = "1 original"), size = pch_size)+
    geom_point(aes(y = ndvi_qc, colour = "2 after QC"), size = pch_size/4)+
    geom_point(data = dat[site == i_site & is.na(ndvi_qc)],
               aes(y = ndvi_f, colour = "3 filled"), size = pch_size)+
    facet_wrap(~year(date), scales = "free_x")+
    scale_x_date(date_breaks = "3 months", date_labels = "%b")+
    scale_shape_manual(values = c(19, 1))+
    theme_bw()+
    ggtitle(i_site)
  
  print(gg)
  
  
})
dev.off()




# gimms ndvi --------------------------------------------------------------



dat <- readRDS("data/gimms-ndvi-01.rds")

all_sites <- sort(unique(dat$site))

pdf(file = "fig/raw-ts/gimms-ndvi.pdf",
    width = 14, height = 8)
walk(all_sites, function(i_site){
  
  gg <- dat[site == i_site] %>% 
    ggplot(aes(xx_half_month/2 + 0.5, ndvi))+
    geom_point()+
    facet_wrap(~year)+
    theme_bw()+
    ggtitle(i_site)+
    xlab("month")+
    scale_x_continuous(breaks = seq(1,12,by=2))
  
  print(gg)
  
  
})
dev.off()



# chronos time (after gimmms 1980) -----------------------------------------------------------------

dat_chrono <- readRDS("data/chrono-02-kh1.rds")
dat_chrono_sub <- dat_chrono[year >= 1980]


for(eli_var in c("TRW", "Kh", "MVA")){
  
  gg <- dat_chrono_sub[variable == eli_var] %>% 
    ggplot(aes(year, value))+
    geom_point()+
    facet_wrap(~site, scales = "free_y")+
    theme_bw()+
    xlab(NULL)+
    ylab(eli_var)
    
  ggsave(gg,
         filename = paste0("fig/raw-ts/chrono/", eli_var, ".pdf"),
         width = 14, height = 7)
}




# chronos time (all) -----------------------------------------------------------------

dat_chrono <- readRDS("data/chrono-02-kh1.rds")
# dat_chrono_sub <- dat_chrono[year >= 1980]


for(eli_var in c("TRW", "Kh", "MVA")){
  
  gg <- dat_chrono[variable == eli_var] %>% 
    ggplot(aes(year, value))+
    geom_vline(xintercept = 1981, linetype = "dashed")+
    geom_point()+
    facet_wrap(~site, scales = "free_y")+
    theme_bw()+
    xlab(NULL)+
    ylab(eli_var)
  
  ggsave(gg,
         filename = paste0("fig/raw-ts/chrono-full/", eli_var, ".pdf"),
         width = 14, height = 7)
}



# chrono scatters ---------------------------------------------------------


dat_chrono <- readRDS("data/chrono-02-kh2.rds")
dat_chrono_sub <- dat_chrono[year >= 1980]



dat_chrono_sub %>% 
  dcast(site + year ~ variable) %>% 
  ggplot(aes(MVA, Kh))+
  geom_point()+
  facet_wrap(~site, scales = "free")+
  theme_bw()

ggsave(filename = "fig/raw-ts/chrono-scatters/MVA-Kh2.pdf",
       width = 12,
       height = 10)

dat_chrono_sub %>% 
  dcast(site + year ~ variable) %>% 
  ggplot(aes(MVA, TRW))+
  geom_point()+
  facet_wrap(~site, scales = "free")+
  theme_bw()

ggsave(filename = "fig/raw-ts/chrono-scatters/MVA-TRW.pdf",
       width = 12,
       height = 10)

dat_chrono_sub %>% 
  dcast(site + year ~ variable) %>% 
  ggplot(aes(TRW, Kh))+
  geom_point()+
  facet_wrap(~site, scales = "free")+
  theme_bw()

ggsave(filename = "fig/raw-ts/chrono-scatters/TRW-Kh2.pdf",
       width = 12,
       height = 10)



# EOF ---------------------------------------------------------------------



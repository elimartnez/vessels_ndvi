# download gimms data

library(gimms)

ofl_v1 <- system.file("extdata", "inventory_ecv1.rds", package = "gimms")
readRDS(ofl_v1)

updateInventory()


downloadGimms(quiet = F, dsn = "/mnt/CEPH_PROJECTS/CLIRSNOW/gimms/")


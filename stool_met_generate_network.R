

library(meconetcomp)
library(microeco)
set.seed(123)
library(magrittr)


stool_met_network <- list()

stool_met_bray <- trans_network$new(dataset = stool_met, cor_method = "bray")
stool_met_bray$cal_network(COR_cut = 0.4)
stool_met_network$stool_met_bray <- stool_met_bray
rm(stool_met_bray)


# Pearson (R >= 0.4, P < 0.05)
stool_met_pear_r4 <- trans_network$new(dataset = stool_met, cor_method = "pearson")
stool_met_pear_r4$cal_network(COR_p_thres = 0.05, COR_cut = 0.4)
stool_met_network$stool_met_pear_r4 <- stool_met_pear_r4
rm(stool_met_pear_r4)
# Pearson (R >= 0.6, P < 0.05)
stool_met_pear_r6 <- trans_network$new(dataset = stool_met, cor_method = "pearson")
stool_met_pear_r6$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)
stool_met_network$stool_met_pear_r6 <- stool_met_pear_r6
rm(stool_met_pear_r6)
# Spearman (R >= 0.4, P < 0.05)
stool_met_spear_r4 <- trans_network$new(dataset = stool_met, cor_method = "spearman")
stool_met_spear_r4$cal_network(COR_p_thres = 0.05, COR_cut = 0.4)
stool_met_network$stool_met_spear_r4 <- stool_met_spear_r4
rm(stool_met_spear_r4)
# Spearman (R >= 0.5, P < 0.05)
stool_met_spear_r5 <- trans_network$new(dataset = stool_met, cor_method = "spearman")
stool_met_spear_r5$cal_network(COR_p_thres = 0.05, COR_cut = 0.5)
stool_met_network$stool_met_spear_r5 <- stool_met_spear_r5
rm(stool_met_spear_r5)
# Spearman (R >= 0.6, P < 0.05)
stool_met_spear_r6 <- trans_network$new(dataset = stool_met, cor_method = "spearman")
stool_met_spear_r6$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)
stool_met_network$stool_met_spear_r6 <- stool_met_spear_r6
rm(stool_met_spear_r6)

stool_met_network %<>% cal_module


save(stool_met_network, file = "stool_met_network.RData")





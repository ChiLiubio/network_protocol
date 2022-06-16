
library(microeco)
set.seed(123)

# create a list to store all the trans_network object
stool_met_network <- list()

# Bray
stool_met_bray <- trans_network$new(dataset = stool_met, cor_method = "bray")
stool_met_bray$cal_network(COR_cut = 0.4)
stool_met_network$stool_met_bray <- stool_met_bray
rm(stool_met_bray)

# pearson (R >= 0.4, P < 0.05)
stool_met_pear_r4 <- trans_network$new(dataset = stool_met, cor_method = "pearson")
stool_met_pear_r4$cal_network(COR_p_thres = 0.05, COR_cut = 0.4)
stool_met_network$stool_met_pear_r4 <- stool_met_pear_r4
rm(stool_met_pear_r4)

# Spearman (R >= 0.4, P < 0.05)
stool_met_spear_r4_p05 <- trans_network$new(dataset = stool_met, cor_method = "spearman")
stool_met_spear_r4_p05$cal_network(COR_p_thres = 0.05, COR_cut = 0.4)
stool_met_network$stool_met_spear_r4_p05 <- stool_met_spear_r4_p05
rm(stool_met_spear_r4_p05)


# Spearman (R >= 0.6, P < 0.05); result is same with (R >= 0.6, P < 0.001)
stool_met_spear_r6_p05 <- clone(stool_met_network$stool_met_spear_r4_p05)
stool_met_spear_r6_p05$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)
stool_met_network$stool_met_spear_r6_p05 <- stool_met_spear_r6_p05
rm(stool_met_spear_r6_p05)


# sparcc
stool_met_sparcc <- trans_network$new(dataset = stool_met, cor_method = "sparcc", use_sparcc_method = "SpiecEasi")
stool_met_sparcc$cal_network(COR_optimization = T, COR_p_thres = 0.05, COR_optimization_low_high = c(0, 0.8))
stool_met_network$stool_met_sparcc <- stool_met_sparcc
rm(stool_met_sparcc)


# SpiecEasi glasso
stool_met_spiec_glasso <- trans_network$new(dataset = stool_met, cor_method = NULL)
stool_met_spiec_glasso$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "glasso")
stool_met_network$stool_met_spiec_glasso <- stool_met_spiec_glasso
rm(stool_met_spiec_glasso)


# SpiecEasi mb
stool_met_spiec_mb <- trans_network$new(dataset = stool_met, cor_method = NULL)
stool_met_spiec_mb$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb")
stool_met_network$stool_met_spiec_mb <- stool_met_spiec_mb
rm(stool_met_spiec_mb)


# gcoda
stool_met_gcoda <- trans_network$new(dataset = stool_met, cor_method = NULL)
stool_met_gcoda$cal_network(network_method = "gcoda")
stool_met_network$stool_met_gcoda <- stool_met_gcoda
rm(stool_met_gcoda)


# FlashWeave
# heterogeneous option
stool_met_flash_het <- trans_network$new(dataset = stool_met, cor_method = NULL)
dir.create("FlashWeave_het_tempdir_het")
stool_met_flash_het$cal_network(network_method = "FlashWeave", FlashWeave_tempdir = "./FlashWeave_het_tempdir_het",
	FlashWeave_other_para = "alpha=0.01,sensitive=true,heterogeneous=true")
stool_met_network$stool_met_flash_het <- stool_met_flash_het
rm(stool_met_flash_het)

stool_met_flash_hom <- trans_network$new(dataset = stool_met, cor_method = NULL)
dir.create("FlashWeave_het_tempdir_hom")
stool_met_flash_hom$cal_network(network_method = "FlashWeave", FlashWeave_tempdir = "./FlashWeave_het_tempdir_hom",
	FlashWeave_other_para = "alpha=0.01,sensitive=true,heterogeneous=false")
stool_met_network$stool_met_flash_hom <- stool_met_flash_hom
rm(stool_met_flash_hom)

# beemStatic
stool_met_bs <- trans_network$new(dataset = stool_met, cor_method = NULL)
stool_met_bs$cal_network(network_method = "beemStatic")
stool_met_network$stool_met_bs <- stool_met_bs
rm(stool_met_bs)


save(stool_met_network, file = "stool_met_network.RData")





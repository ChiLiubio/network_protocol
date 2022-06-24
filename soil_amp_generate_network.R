

library(magrittr)
library(dplyr)
library(file2meco)
library(microeco)
set.seed(123)


soil_amp_network <- list()


# Bray

soil_amp_bray <- trans_network$new(dataset = soil_amp, cor_method = "bray")
soil_amp_bray$cal_network(COR_cut = 0.5)
soil_amp_network$soil_amp_bray <- soil_amp_bray
rm(soil_amp_bray)

# pearson (R >= 0.6, P < 0.05)
soil_amp_pear_r6 <- trans_network$new(dataset = soil_amp, cor_method = "pearson", use_WGCNA_pearson_spearman = TRUE)
soil_amp_pear_r6$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)
soil_amp_network$soil_amp_pear_r6 <- soil_amp_pear_r6

rm(soil_amp_pear_r6)


# Spearman (R >= 0.5, p < 0.01)
soil_amp_spear_r5 <- trans_network$new(dataset = soil_amp, cor_method = "spearman", use_WGCNA_pearson_spearman = TRUE)
soil_amp_spear_r5$cal_network(COR_p_thres = 0.01, COR_cut = 0.5)
soil_amp_network$soil_amp_spear_r5 <- soil_amp_spear_r5

# Spearman (R >= 0.6, p < 0.01)
soil_amp_spear_r6 <- clone(soil_amp_spear_r5)
soil_amp_spear_r6$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
soil_amp_network$soil_amp_spear_r6 <- soil_amp_spear_r6

# Spearman (R >= 0.7, p < 0.01)
soil_amp_spear_r7 <- clone(soil_amp_spear_r5)
soil_amp_spear_r7$cal_network(COR_p_thres = 0.01, COR_cut = 0.7)
soil_amp_network$soil_amp_spear_r7 <- soil_amp_spear_r7


rm(soil_amp_spear_r5)
rm(soil_amp_spear_r6)
rm(soil_amp_spear_r7)


# SparCC
soil_amp_sparcc <- trans_network$new(dataset = soil_amp, cor_method = "sparcc", use_sparcc_method  = "NetCoMi")
soil_amp_sparcc$cal_network(COR_optimization = TRUE)
soil_amp_network$soil_amp_sparcc <- soil_amp_sparcc
rm(soil_amp_sparcc)



# SpiecEasi mb
soil_amp_spiec_mb <- trans_network$new(dataset = soil_amp, cor_method = NULL)
soil_amp_spiec_mb$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb")
soil_amp_network$soil_amp_spiec_mb <- soil_amp_spiec_mb

rm(soil_amp_spiec_mb)


# FlashWeave with metadata
soil_amp_flash_meta <- trans_network$new(dataset = soil_amp, cor_method = NULL, env_cols = 8:15)
dir.create("FlashWeave_het_tempdir_meta")
soil_amp_flash_meta$cal_network(network_method = "FlashWeave", FlashWeave_meta_data = TRUE, FlashWeave_tempdir = "./FlashWeave_het_tempdir_meta",
	FlashWeave_other_para = "alpha=0.01,sensitive=true,heterogeneous=true")
soil_amp_flash_meta$res_network <- soil_amp_flash_meta$subset_network(node = rownames(soil_amp$otu_table))
soil_amp_network$soil_amp_flash_meta <- soil_amp_flash_meta

# FlashWeave without metadata
soil_amp_flash_nom <- trans_network$new(dataset = soil_amp, cor_method = NULL)
dir.create("FlashWeave_het_tempdir_nometa")
soil_amp_flash_nom$cal_network(network_method = "FlashWeave", FlashWeave_meta_data = FALSE, FlashWeave_tempdir = "./FlashWeave_het_tempdir_nometa",
	FlashWeave_other_para = "alpha=0.01,sensitive=true,heterogeneous=false")
soil_amp_network$soil_amp_flash_nom <- soil_amp_flash_nom

rm(soil_amp_flash_meta)
rm(soil_amp_flash_nom)


# beemStatic network
soil_amp_bs <- trans_network$new(dataset = soil_amp, cor_method = NULL)
soil_amp_bs$cal_network(network_method = "beemStatic")
soil_amp_network$soil_amp_bs <- soil_amp_bs
rm(soil_amp_bs)

for(i in names(soil_amp_network)){
	soil_amp_network[[i]]$res_cor_p <- NULL
}

save(soil_amp_network, file = "soil_amp_network.RData")







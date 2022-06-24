

library(pacman)
p_load(microeco, magrittr, igraph, dplyr, ggplot2, grid, tibble)
set.seed(123)
theme_set(theme_bw())


data("stool_met.RData")
data("stool_met_network.RData")
data("prok_func_NJC19_list")

plot_edge_metabolite <- function(label = NULL, substance = NULL, ylab_name = NULL){
	long_name <- c("Indole", "Isobutyrate (2-Methylpropanoic acid)", "Isovalerate (3-Methylbutanoic acid)", "Niacin (Vitamin B3, Nicotinic acid, Nicotinate, Nicotinamide, Vitamin B3)", 
	"Pyruvate", "Sulfate (Sulfuric acid, Sulfite)", "Menaquinone (Vitamin K2)", "Methylamine (Monomethylamine)", "Propanoate (Propionate)", "Urea", 
	"Butyrate", "Chenodeoxycholic acid (Chenodeoxycholate)", "Cholic acid (Cholate)", 
	"Methanol", "Riboflavin (Vitamin B2)", "1,2-propanediol (Propene diol, Propylene glycol, [R]-1,2-propanediol, [R]-propane-1,2-diol, [S]-1,2-propanediol, [S]-propane-1,2-diol)", 
	"[R,R]-2,3-Butanediol ([R,R]-Butanediol, [S,S]-2,3-Butanediol, [S,S]-Butanediol)", "Ethanol", "H2 (Hydrogen)", "Succinate", "Formate", "Thiamine (Vitamin B1, Thiamin)", 
	"NH3 (Ammonia, NH4+, Ammonium)", "CO2", "Folic acid (Folate, Vitamin B9)", "Biotin (Vitamin B7)", "Pyridoxal (Vitamin B6, Pyridoxine, Pyridoxamine, Vitamin B6)", 
	"Adenosylcobalamin (Vitamin B12, Cobamide coenzyme, Cob[II]alamin, Vitamin B12r, Cob[I]alamin, Vitamin B12s, Aquacobalamin, Cobalamin [III])", 
	"L-Lactate ([S]-Lactate, Lactate, D-Lactate, [R]-Lactate)", "Acetate")
	short_name <- c("Indole", "Isobutyrate (2-Methylpropanoic acid)", "Isovalerate (3-Methylbutanoic acid)", 
	"Niacin (Vitamin B3, Nicotinate, Nicotinamide)", "Pyruvate", "Sulfate (Sulfuric acid, Sulfite)", "Menaquinone (Vitamin K2)", 
	"Methylamine (Monomethylamine)", "Propanoate (Propionate)", "Urea", "Butyrate", "Chenodeoxycholic acid (Chenodeoxycholate)", "Cholic acid (Cholate)", 
	"Methanol", "Riboflavin (Vitamin B2)", "1,2-propanediol (Propylene glycol)", "[R,R]-2,3-Butanediol ([S,S]-2,3-Butanediol)", 
	"Ethanol", "H2 (Hydrogen)", "Succinate", "Formate", "Thiamine (Vitamin B1, Thiamin)", "NH3 (Ammonia, NH4+, Ammonium)", "CO2", "Folic acid (Folate, Vitamin B9)", 
	"Biotin (Vitamin B7)", "Pyridoxal (Vitamin B6, Pyridoxine, Pyridoxamine)", "Adenosylcobalamin (Vitamin B12)", 
	"L-Lactate ([S]-Lactate, D-Lactate, [R]-Lactate)", "Acetate")

	# species_edge_metabolism to store the match results
	species_edge_metabolism <- data.frame()
	for(i in names(stool_met_network)){
		tmp <- stool_met_network[[i]]$res_edge_table
		tmp1 <- tmp
		tmp1[, 1] <- tmp[, 2]
		tmp1[, 2] <- tmp[, 1]
		tmp <- rbind(tmp, tmp1)
		tmp$net <- i
		species_edge_metabolism <- rbind(species_edge_metabolism, tmp)
	}
	species_edge_metabolism$for_match <- paste0(species_edge_metabolism[, 1], "-", species_edge_metabolism[, 2])
	species_edge_metabolism %<>% .[.$label == label, ]

	species_edge_metabolism$PS <- NA
	species_edge_metabolism$SS <- NA
	species_edge_metabolism$PP <- NA
	species_edge_metabolism %<>% .[.[, 1] %in% names(prok_func_NJC19_list), ]
	species_edge_metabolism %<>% .[.[, 2] %in% names(prok_func_NJC19_list), ]

	for(i in seq_len(nrow(species_edge_metabolism))){
		tmp1 <- prok_func_NJC19_list[[species_edge_metabolism[i, 1]]]
		tmp2 <- prok_func_NJC19_list[[species_edge_metabolism[i, 2]]]
		species_edge_metabolism[i, "PS"] <- intersect(tmp1$Production, tmp2$Consumption) %>% paste0(collapse = "|")
		species_edge_metabolism[i, "SS"] <- intersect(tmp1$Consumption, tmp2$Consumption) %>% paste0(collapse = "|")
		species_edge_metabolism[i, "PP"] <- intersect(tmp1$Production, tmp2$Production) %>% paste0(collapse = "|")
	}
	# parse required
	tmp <- species_edge_metabolism %>% .[, substance] %>% lapply(., function(x){unlist(strsplit(x, "\\|"))}) %>% unlist %>% table %>% names
	tmp <- data.frame(name = tmp)
	for(i in names(stool_met_network)){
		if(i %in% species_edge_metabolism$net){
			tmp1 <- species_edge_metabolism %>% .[.$net == i, substance] %>% unique
			if(identical(tmp1, "")){
				next
			}else{
				tmp1 %<>% lapply(., function(x){unlist(strsplit(x, "\\|"))}) %>% unlist %>% table %>% as.data.frame
				colnames(tmp1) <- c("name", i)
				tmp <- left_join(tmp, tmp1, by = c("name" = "name"))
			}
		}else{
			next
		}
	}
	rownames(tmp) <- tmp$name
	tmp <- tmp[, -1]
	tmp <- tmp[long_name, ]
	rownames(tmp) <- short_name
	tmp[is.na(tmp)] <- 0
	tmp <- tmp[apply(tmp, 1, sum) != 0, ]

	for(i in names(stool_met_network)){
		tmp3 <- stool_met_network[[i]]$res_edge_table
		tmp3 %<>% .[.$label == label, ]
		if(i %in% colnames(tmp)){
			tmp[, i] %<>% {. / nrow(tmp3)}
		}else{
			next
		}
	}

	tmp1 <- reshape2::melt(rownames_to_column(tmp, var = "name"), id.vars = "name")
	colnames(tmp1)[3] <- "Abundance"
	tmp1$Abundance %<>% {. * 100}
	tmp2 <- reshape2::dcast(tmp1, name~variable, value.var = "Abundance") %>% `row.names<-`(.[,1]) %>% .[, -1, drop = FALSE]
	tmp2[is.na(tmp2)] <- 0
	x1 <- hclust(dist(tmp2)) 
	x1 %<>% {.$labels[.$order]}
	x2 <- names(stool_met_network)

	tmp2 <- tmp1
	tmp2$Abundance %<>% round(., 1)

	g <- ggplot(tmp1, aes_string(x = "variable", y = "name", label = formatC("Abundance", format = "f", digits = 2))) +
		geom_tile(aes(fill = Abundance), colour = "white", size = 0.5) +
		theme(axis.text.y = element_text(size = 11)) + theme(plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")) + 
		geom_text(data = tmp2, size = 4, colour = "grey10") + 
		scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")[-11]), trans = "log10", na.value = RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")[11]) + 
		labs(x = "", y = ylab_name, fill = "Percentage") +
		theme(axis.text.x = element_text(angle = 25, colour = "black", vjust = 1, hjust = 1, size = 11)) + 
		theme(axis.title.y = element_text(size = 15)) +
		scale_y_discrete(limits = x1) + scale_x_discrete(limits = x2)

	g
}

g1 <- plot_edge_metabolite(label = "+", substance = "PS", ylab_name = "Metabolites for exchange in positive edges")
g2 <- plot_edge_metabolite(label = "+", substance = "SS", ylab_name = "Substrates for competition in positive edges")
g3 <- plot_edge_metabolite(label = "-", substance = "PS", ylab_name = "Metabolites for exchange in negative edges")
g4 <- plot_edge_metabolite(label = "-", substance = "SS", ylab_name = "Substrates for competition in negative edges")


